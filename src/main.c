#include <stdio.h>
#include <math.h>

#include <include/sim/env.h>
#include <include/utils.h>
#include <math/linalg.h>
#include <include/poisson.h>

#include <include/fem/functions.h>
#include <include/sim/bulk.h>
#include <include/sim/dopant.h>

#include <math/pyvisual.h>

#include <include/fem/mesh.h>

double bulkNetCharge(Bulk bulk, Dopant* dopants, size_t dopant_count, const double x, const double fermi_lvl, const double Ec, Environment env)
{
    double charge = -ELECTRON_CHARGE * bulkElectrons(bulk, fermi_lvl, Ec, env);

    charge += ELECTRON_CHARGE * bulkHoles(bulk, fermi_lvl, Ec - bulk.bandgap, env);

    for(size_t i = 0; i < dopant_count; i++)
    {
        charge += dopantIonized(dopants[i], fermi_lvl, Ec, x, env);
    }
    return charge;
}

void bulkNetChargeVec(Bulk bulk, Dopant* dopants, size_t dopant_count, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* net_charge)
{
    Vec charge = vecInitZerosA(net_charge->len);

    bulkElectronsV(bulk, fermi_lvl, Ec, env, net_charge);
    vecScale(-ELECTRON_CHARGE, *net_charge, net_charge);

    bulkHolesV(bulk, fermi_lvl, Ec, env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge, charge, net_charge);

    for(size_t i = 0; i < dopant_count; i++)
    {
        dopantIonizedV(dopants[i], fermi_lvl, Ec, env, &charge);
        vecAdd(*net_charge, charge, net_charge);
    }

    freeVec(&charge);
}

void bulkNetChargeVecDerivative(Bulk bulk, Dopant* dopants, size_t dopant_count, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* net_charge_d)
{
    Vec charge = vecInitZerosA(net_charge_d->len);

    bulkElectronsDV(bulk, fermi_lvl, Ec, env, net_charge_d);
    vecScale(-ELECTRON_CHARGE, *net_charge_d, net_charge_d);

    bulkHolesDV(bulk, fermi_lvl, Ec, env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge_d, charge, net_charge_d);

    for(size_t i = 0; i < dopant_count; i++)
    {
        dopantIonizedDV(dopants[i], fermi_lvl, Ec, env, &charge);
        vecAdd(*net_charge_d, charge, net_charge_d);
    }

    freeVec(&charge);
}

double bulkMinConc(Bulk bulk, Dopant* dopants, size_t dopant_count, const double x, const double fermi_lvl, const double Ec, Environment env)
{
    double minConc = bulkElectrons(bulk, fermi_lvl, Ec, env);

    minConc = min_d(bulkHoles(bulk, fermi_lvl, Ec - bulk.bandgap, env), minConc);

    for(size_t i = 0; i < dopant_count; i++)
    {
        minConc = min_d(dopantIonizedConc(dopants[i], fermi_lvl, Ec, x, env), minConc);
    }
    return minConc;
}

// solve for fermi-level at x = 0, assuming Ec = 0 at x = 0(ref)
double solveFermiLvl(Bulk bulk, Dopant* dopants, size_t dopant_count, Environment env)
{
    double low_u = 0;
    double high_u = -bulk.bandgap;

    double low_value = bulkNetCharge(bulk, dopants, dopant_count, 0.0, low_u, 0.0, env);
    double high_value = bulkNetCharge(bulk, dopants, dopant_count, 0.0, high_u, 0.0, env);

    double min_conc = bulkMinConc(bulk, dopants, dopant_count, 0.0, low_u, 0.0, env);
    min_conc = min_d(bulkMinConc(bulk, dopants, dopant_count, 0.0, high_u, 0.0, env), min_conc);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_u - low_u) > env.tol_potential*fabs(high_u + low_u))
    {
        double new_u = (low_u + high_u) / 2;

        double error = bulkNetCharge(bulk, dopants, dopant_count, 0.0, new_u, 0.0, env);
        min_conc = bulkMinConc(bulk, dopants, dopant_count, 0.0, new_u, 0.0, env);

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
double solveEc(Bulk bulk, Dopant* dopants, size_t dopant_count, Environment env, double fermi_lvl)
{
    double low_Ec = fermi_lvl;
    double high_Ec = fermi_lvl + bulk.bandgap;

    double low_value = bulkNetCharge(bulk, dopants, dopant_count, env.end_x, fermi_lvl, low_Ec, env);
    double high_value = bulkNetCharge(bulk, dopants, dopant_count, env.end_x, fermi_lvl, high_Ec, env);

    double min_conc = bulkMinConc(bulk, dopants, dopant_count, env.end_x, fermi_lvl, low_Ec, env);
    min_conc = min_d(bulkMinConc(bulk, dopants, dopant_count, env.end_x, fermi_lvl, high_Ec, env), min_conc);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_Ec - low_Ec) > env.tol_potential*fabs(high_Ec + low_Ec))
    {
        double new_Ec = (low_Ec + high_Ec) / 2;

        double error = bulkNetCharge(bulk, dopants, dopant_count, env.end_x, fermi_lvl, new_Ec, env);
        min_conc = bulkMinConc(bulk, dopants, dopant_count, env.end_x, fermi_lvl, new_Ec, env);

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
    Environment env;
    env.temperature = 300;
    env.tol_potential = 1e-4;
    env.tol_conc = 1e-6;
    env = envPopulate(env);

    Bulk silicon;
    silicon.bandgap = 1.1;
    silicon.effElectronMass = 0.98 * ELECTRON_MASS;
    silicon.effHoleMass = 0.48 * ELECTRON_MASS;
    silicon.epsilon = 11.68 * VACCUME_PERMITTIVITY; // convert to C/(eV m)

    double n_doping = 1e21;
    Vec conc = vecConstruct((double[]){0.0, n_doping}, 2);
    Vec    x = vecConstruct((double[]){0.0, env.end_x}, 2);
    Fx1D phosphorous_conc = fxConstruct1D(conc, x, INTERP_NEAREST);
    Dopant phosphorous = dopantInitDonor1DA(phosphorous_conc, 0.045, 2, silicon, env);

    double p_doping = 1e21;
    conc = vecConstruct((double[]){p_doping, 0.0}, 2);
    x    = vecConstruct((double[]){0.0 , env.end_x}, 2);
    Fx1D boron_conc = fxConstruct1D(conc, x, INTERP_NEAREST);
    Dopant boron = dopantInitAcceptor1DA(boron_conc, 0.045, 4, silicon, env);

    Dopant dopants[] = {phosphorous, boron};

    double fermi_lvl = solveFermiLvl(silicon, dopants, 2, env);
    double Ec = solveEc(silicon, dopants, 2, env, fermi_lvl);

    Vec fermi_lvl_vec = vecInitA(fermi_lvl, env.sample_count);

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

    PyViBase pos = pyviCreateParameter(&pyvi, "position", env.x_samples);
    PyViSec pyvi_pot = pyviCreateSection(&pyvi, "potential", pos);
    PyViSec pyvi_ch = pyviCreateSection(&pyvi, "charge", pos);
    PyViSec pyvi_chd = pyviCreateSection(&pyvi, "charge derivative", pos);
    PyViSec pyvi_res = pyviCreateSection(&pyvi, "residual", pos);
    PyViSec pyvi_J = pyviCreateSection(&pyvi, "Jacobian", pos);
    PyViSec pyvi_F = pyviCreateSection(&pyvi, "Fermi lvl", pos);

    // setup initial guess

    double xn = (env.start_x + env.end_x)/2.0 + sqrt(2 * silicon.epsilon / ELECTRON_CHARGE * (p_doping / n_doping) /( p_doping + n_doping ) * (-Ec));
    double xp = (env.start_x + env.end_x)/2.0 - sqrt(2 * silicon.epsilon / ELECTRON_CHARGE * (n_doping / p_doping) /( p_doping + n_doping ) * (-Ec));

    printf("x_n: %le\n", xn - (env.start_x + env.end_x)/2.0);
    printf("x_p: %le\n", xp - (env.start_x + env.end_x)/2.0);

    Fx1D guessRho = fxConstruct1D(vecConstruct((double[]){  0.0, -ELECTRON_CHARGE * p_doping, ELECTRON_CHARGE * n_doping, 0.0 }, 4), 
                                  vecConstruct((double[]){ xp - env.spacing,        xp,                         xn, xn + env.spacing }, 4), INTERP_NEAREST);
    
    // solve for potential
    fxInterpolateSample1D(guessRho, env.x_samples, &rho);
    //poisson1DSolve(&V, rho, silicon.epsilon, 0.0, -Ec, env.spacing);

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
    Fx1D fxPotential = fxConstruct1D(V, env.x_samples, INTERP_LINEAR);
    //poissonEvaluate(V, &accum, 0.0, -Ec, silicon.epsilon, env.spacing);
    //vecScale(silicon.epsilon / (env.spacing*env.spacing), accum, &accum);

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
        poissonEvaluate(V, &accum, 0.0, -Ec, silicon.epsilon, env.spacing);
        bulkNetChargeVec(silicon, dopants, 2, fermi_lvl_vec, Ec_vec, env, &rho);
        vecScale(env.spacing*env.spacing / silicon.epsilon, rho, &ch_scaled);

        vecSub(accum, ch_scaled, &accum);

        bulkNetChargeVecDerivative(silicon, dopants, 2, fermi_lvl_vec, Ec_vec, env, &rho_d);
        vecScale(env.spacing*env.spacing / silicon.epsilon, rho_d, &ch_scaled);
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

    freeSim(&env);
}
