#include <stdio.h>
#include <math.h>

#include <include/sim.h>
#include <include/utils.h>
#include <math/linalg.h>
#include <include/poisson.h>

#include <include/fem/functions.h>
#include <include/sim/material.h>

#include <math/pyvisual.h>

#include <include/fem/mesh.h>

double bulkNetCharge(Bulk bulk, Dopant* dopants, size_t dopant_count, const double x, const double fermi_lvl, const double Ec, SimParameters simParams)
{
    double charge = -ELECTRON_CHARGE * bulkElectronConcSingle(bulk, fermi_lvl, Ec, simParams.temperature, simParams.tol_conc);

    charge += ELECTRON_CHARGE * bulkHoleConcSingle(bulk, fermi_lvl, Ec - bulk.bandgap, simParams.temperature, simParams.tol_conc);

    for(size_t i = 0; i < dopant_count; i++)
    {
        charge += dopantIonizedChargeSingle(dopants[i], fermi_lvl, Ec, x, simParams);
    }
    return charge;
}

void bulkNetChargeVec(Bulk bulk, Dopant* dopants, size_t dopant_count, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* net_charge)
{
    Vec charge = vecInitZerosA(net_charge->len);

    bulkElectronConc(bulk, fermi_lvl, Ec, simParams, net_charge);
    vecScale(-ELECTRON_CHARGE, *net_charge, net_charge);

    bulkHoleConc(bulk, fermi_lvl, Ec, simParams, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge, charge, net_charge);

    for(size_t i = 0; i < dopant_count; i++)
    {
        dopantIonizedCharge(dopants[i], fermi_lvl, Ec, simParams, &charge);
        vecAdd(*net_charge, charge, net_charge);
    }

    freeVec(&charge);
}

void bulkNetChargeVecDerivative(Bulk bulk, Dopant* dopants, size_t dopant_count, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* net_charge_d)
{
    Vec charge = vecInitZerosA(net_charge_d->len);

    bulkElectronConcDerivative(bulk, fermi_lvl, Ec, simParams, net_charge_d);
    vecScale(-ELECTRON_CHARGE, *net_charge_d, net_charge_d);

    bulkHoleConcDerivative(bulk, fermi_lvl, Ec, simParams, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge_d, charge, net_charge_d);

    for(size_t i = 0; i < dopant_count; i++)
    {
        dopantIonizedChargeDerivative(dopants[i], fermi_lvl, Ec, simParams, &charge);
        vecAdd(*net_charge_d, charge, net_charge_d);
    }

    freeVec(&charge);
}

double bulkMinConc(Bulk bulk, Dopant* dopants, size_t dopant_count, const double x, const double fermi_lvl, const double Ec, SimParameters simParams)
{
    double minConc = bulkElectronConcSingle(bulk, fermi_lvl, Ec, simParams.temperature, simParams.tol_conc);

    minConc = min_d(bulkHoleConcSingle(bulk, fermi_lvl, Ec - bulk.bandgap, simParams.temperature, simParams.tol_conc), minConc);

    for(size_t i = 0; i < dopant_count; i++)
    {
        minConc = min_d(dopantIonizedConcSingle(dopants[i], fermi_lvl, Ec, x, simParams), minConc);
    }
    return minConc;
}

// solve for fermi-level at x = 0, assuming Ec = 0 at x = 0(ref)
double solveFermiLvl(Bulk bulk, Dopant* dopants, size_t dopant_count, SimParameters simParams)
{
    double low_u = 0;
    double high_u = -bulk.bandgap;

    double low_value = bulkNetCharge(bulk, dopants, dopant_count, 0.0, low_u, 0.0, simParams);
    double high_value = bulkNetCharge(bulk, dopants, dopant_count, 0.0, high_u, 0.0, simParams);

    double min_conc = bulkMinConc(bulk, dopants, dopant_count, 0.0, low_u, 0.0, simParams);
    min_conc = min_d(bulkMinConc(bulk, dopants, dopant_count, 0.0, high_u, 0.0, simParams), min_conc);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_u - low_u) > simParams.tol_potential*fabs(high_u + low_u))
    {
        double new_u = (low_u + high_u) / 2;

        double error = bulkNetCharge(bulk, dopants, dopant_count, 0.0, new_u, 0.0, simParams);
        min_conc = bulkMinConc(bulk, dopants, dopant_count, 0.0, new_u, 0.0, simParams);

        if(error*low_value < 0)
        {
            high_u = new_u;
            high_value = error;
        }
        else
        {
            low_u = new_u;
            low_value = error;
        }
    }

    return (high_u + low_u) / 2;
}

// solve for Ec at x = L
double solveEc(Bulk bulk, Dopant* dopants, size_t dopant_count, SimParameters simParams, double fermi_lvl)
{
    double low_Ec = fermi_lvl;
    double high_Ec = fermi_lvl + bulk.bandgap;

    double low_value = bulkNetCharge(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, low_Ec, simParams);
    double high_value = bulkNetCharge(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, high_Ec, simParams);

    double min_conc = bulkMinConc(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, low_Ec, simParams);
    min_conc = min_d(bulkMinConc(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, high_Ec, simParams), min_conc);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_Ec - low_Ec) > simParams.tol_potential*fabs(high_Ec + low_Ec))
    {
        double new_Ec = (low_Ec + high_Ec) / 2;

        double error = bulkNetCharge(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, new_Ec, simParams);
        min_conc = bulkMinConc(bulk, dopants, dopant_count, simParams.end_x, fermi_lvl, new_Ec, simParams);

        if(error*low_value < 0)
        {
            high_Ec = new_Ec;
            high_value = error;
        }
        else
        {
            low_Ec = new_Ec;
            low_value = error;
        }
    }

    return (high_Ec + low_Ec) / 2;
}

int main()
{
    size_t N = 1 << 8;

    // make a simulation setup
    SimParameters simParams = simConstructA(N, 0.0, 0.2e-5, 300, 1e-4, 1e-6);

    Bulk silicon;

    silicon.bandgap = 1.1;
    silicon.effMassElectron = 0.98 * ELECTRON_MASS;
    silicon.effMassHoles = 0.48 * ELECTRON_MASS;
    silicon.epsilon = 11.68 * VACCUME_PERMITTIVITY; // convert to C/(eV m)

    double n_doping = 1e21;
    Vec conc = vecConstruct((double[]){0.0, n_doping}, 2);
    Vec    x = vecConstruct((double[]){0.0, simParams.end_x}, 2);
    Fx1D phosphorous_conc = fxConstruct1D(conc, x, INTERP_NEAREST);
    Dopant phosphorous = dopantConstructDonor1DA(phosphorous_conc, 0.045, 2, silicon, simParams);

    double p_doping = 1e21;
    conc = vecConstruct((double[]){p_doping, 0.0}, 2);
    x    = vecConstruct((double[]){0.0 , simParams.end_x}, 2);
    Fx1D boron_conc = fxConstruct1D(conc, x, INTERP_NEAREST);
    Dopant boron = dopantConstructAcceptor1DA(boron_conc, 0.045, 4, silicon, simParams);

    Dopant dopants[] = {phosphorous, boron};

    double fermi_lvl = solveFermiLvl(silicon, dopants, 2, simParams);
    double Ec = solveEc(silicon, dopants, 2, simParams, fermi_lvl);

    Vec fermi_lvl_vec = vecInitA(fermi_lvl, simParams.sample_count);

    printf("Fermi-Level: %lf eV\n", fermi_lvl);
    printf("Built-in Potential: %lf V\n", -Ec);

    Vec rho = vecInitZerosA(N);
    Vec rho_tmp = vecInitZerosA(N);
    Vec rho_d = vecInitZerosA(N);
    Vec V = vecInitZerosA(N);
    Vec Ec_vec = vecInitZerosA(N);
    Vec res = vecInitZerosA(N);

    Vec subdiag_jacobian = vecInitA(-1, N);
    Vec ref_diag_jacobian = vecInitA(2.0, N);
    Vec diag_jacobian = vecInitA(0.0, N);

    Vec scratch = vecInitZerosA(N);

    PyVi pyvi = pyviInitA("v.pyvi");

    PyViBase pos = pyviCreateParameter(&pyvi, "position", simParams.x_samples);
    PyViSec pyvi_pot = pyviCreateSection(&pyvi, "potential", pos);
    PyViSec pyvi_ch = pyviCreateSection(&pyvi, "charge", pos);
    PyViSec pyvi_chd = pyviCreateSection(&pyvi, "charge derivative", pos);
    PyViSec pyvi_res = pyviCreateSection(&pyvi, "residual", pos);
    PyViSec pyvi_J = pyviCreateSection(&pyvi, "Jacobian", pos);
    PyViSec pyvi_F = pyviCreateSection(&pyvi, "Fermi lvl", pos);

    // setup initial guess

    double xn = (simParams.start_x + simParams.end_x)/2.0 + sqrt(2 * silicon.epsilon / ELECTRON_CHARGE * (p_doping / n_doping) /( p_doping + n_doping ) * (-Ec));
    double xp = (simParams.start_x + simParams.end_x)/2.0 - sqrt(2 * silicon.epsilon / ELECTRON_CHARGE * (n_doping / p_doping) /( p_doping + n_doping ) * (-Ec));

    printf("x_n: %le\n", xn - (simParams.start_x + simParams.end_x)/2.0);
    printf("x_p: %le\n", xp - (simParams.start_x + simParams.end_x)/2.0);

    Fx1D guessRho = fxConstruct1D(vecConstruct((double[]){  0.0, -ELECTRON_CHARGE * p_doping, ELECTRON_CHARGE * n_doping, 0.0 }, 4), 
                                  vecConstruct((double[]){ xp - simParams.spacing,        xp,                         xn, xn + simParams.spacing }, 4), INTERP_NEAREST);
    
    // solve for potential
    fxInterpolateSample1D(guessRho, simParams.x_samples, &rho);
    //poisson1DSolve(&V, rho, silicon.epsilon, 0.0, -Ec, simParams.spacing);

    Mesh mesh = meshInitPieceUniformA(
        vecConstruct((double[]){0.0, xp, xn, 0.2e-5}, 4),
        (size_t[]){N / 8,  3 * (N / 4), N / 8},
        3
    );

    meshSetDirichletBC(&mesh, 0.0, -Ec);

    PyViBase pos_mesh = pyviCreateParameter(&pyvi, "position_nu", mesh.x);
    PyViSec pyvi_M = pyviCreateSection(&pyvi, "Mesh Test", pos_mesh);
    
    Vec accum = vecInitZerosA(N);
    Vec mesh_test = vecInitZerosA(N);
    Fx1D fxPotential = fxConstruct1D(V, simParams.x_samples, INTERP_LINEAR);
    //poissonEvaluate(V, &accum, 0.0, -Ec, silicon.epsilon, simParams.spacing);
    //vecScale(silicon.epsilon / (simParams.spacing*simParams.spacing), accum, &accum);

    pyviSectionPush(pyvi_pot, V);
    pyviSectionPush(pyvi_ch, rho);
    pyviSectionPush(pyvi_chd, rho_d);
    pyviSectionPush(pyvi_J, diag_jacobian);
    pyviSectionPush(pyvi_F, fermi_lvl_vec);
    pyviSectionPush(pyvi_res, accum);
    pyviSectionPush(pyvi_M, mesh_test);

    vecPrint(mesh.dx);
    printf(":len = %zu\n", mesh.len);

    for(size_t i = 0; i < 10; i++)
    {
        Vec ch_scaled = vecInitZerosA(N);

        vecScale(-1, V, &Ec_vec);
        poissonEvaluate(V, &accum, 0.0, -Ec, silicon.epsilon, simParams.spacing);
        bulkNetChargeVec(silicon, dopants, 2, fermi_lvl_vec, Ec_vec, simParams, &rho);
        vecScale(simParams.spacing*simParams.spacing / silicon.epsilon, rho, &ch_scaled);

        vecSub(accum, ch_scaled, &accum);

        bulkNetChargeVecDerivative(silicon, dopants, 2, fermi_lvl_vec, Ec_vec, simParams, &rho_d);
        vecScale(simParams.spacing*simParams.spacing / silicon.epsilon, rho_d, &ch_scaled);
        vecAdd(ref_diag_jacobian, ch_scaled, &diag_jacobian);

        solveTridiagonalSymm(accum, subdiag_jacobian, diag_jacobian, scratch);
        vecSub(V, accum, &accum);

        //vecCopy(V, &mesh.potential);
        fxInterpolateSample1D(fxPotential, mesh.x, &mesh.potential);
        meshPoissonEvaluate(mesh, silicon.epsilon, &mesh_test);

        pyviSectionPush(pyvi_pot, V);
        pyviSectionPush(pyvi_ch, rho);
        pyviSectionPush(pyvi_chd, rho_d);
        pyviSectionPush(pyvi_J, diag_jacobian);
        pyviSectionPush(pyvi_F, fermi_lvl_vec);
        pyviSectionPush(pyvi_res, accum);
        pyviSectionPush(pyvi_M, mesh_test);

        // set V to new value
        vecCopy(accum, &V);
    }

    pyviWrite(pyvi);

    freePyVi(&pyvi);

    freeMesh(&mesh);

    freeDopant(&phosphorous);
    freeDopant(&boron);

    freeSim(&simParams);
}
