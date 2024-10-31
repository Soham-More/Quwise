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

#include <include/sim/semiconductor.h>

int main()
{
    size_t N = 1 << 8;

    double length = 2e-5;
    double n_doping = 1e21;
    double p_doping = 1e21;

    // make a simulation setup
    Environment env;
    env.temperature = 300;
    env.tol_potential = 1e-4;
    env.tol_conc = 1e-6;
    env = envPopulate(env);

    BulkInfo siliconInfo;
    siliconInfo.bandgap = 1.1;
    siliconInfo.relativeEffElectronMass = 0.98;
    siliconInfo.relativeEffHoleMass = 0.48;
    siliconInfo.relativePermittivity = 11.68;

    DopantInfo phosphorousInfo;
    phosphorousInfo.sampled_conc = vecConstruct((double[]){0.0, n_doping}, 2);
    phosphorousInfo.sampled_x    = vecConstruct((double[]){0.0, length}, 2);
    phosphorousInfo.delE = 0.045;
    phosphorousInfo.degeneracy = 2;
    phosphorousInfo.samplingMode = INTERP_NEAREST;

    DopantInfo boronInfo;
    boronInfo.sampled_conc = vecConstruct((double[]){p_doping, 0.0}, 2);
    boronInfo.sampled_x    = vecConstruct((double[]){0.0 , length}, 2);
    boronInfo.delE = 0.045;
    boronInfo.degeneracy = 2;
    boronInfo.samplingMode = INTERP_NEAREST;

    SemiConductor silicon = scInit(siliconInfo, env);

    double v_estimate = env.thermal_potential * (10 * log(10));

    double xn = length/2.0 + sqrt(2 * silicon.bulk.epsilon / ELECTRON_CHARGE * (p_doping / n_doping) /( p_doping + n_doping ) * v_estimate);
    double xp = length/2.0 - sqrt(2 * silicon.bulk.epsilon / ELECTRON_CHARGE * (n_doping / p_doping) /( p_doping + n_doping ) * v_estimate);

    printf("x_n: %le\n", xn - length/2.0);
    printf("x_p: %le\n", xp - length/2.0);

    Mesh mesh = meshInitPieceUniformA(
        vecConstruct((double[]){0.0, xp, xn, 0.2e-5}, 4),
        (size_t[]){N / 8,  3 * (N / 4), N / 8},
        3
    );

    scDopeAcceptor(&silicon, mesh, phosphorousInfo);
    scDopeDonor(&silicon, mesh, boronInfo);

    // solve for actual fermi and conduction level
    double boundary_ef = scSolveBoundary(silicon, 0.0, 0.0, SC_SOLVE_FERMI);
    double boundary_ec = scSolveBoundary(silicon, 0.0, boundary_ef, SC_SOLVE_EC);

    printf("Fermi-Level: %lf eV\n", boundary_ef);
    printf("Built-in Potential: %lf V\n", -boundary_ec);

    meshSetDirichletBC(&mesh, 0.0, -boundary_ec);

    Vec fermi_lvl = vecInitA(boundary_ef, mesh.len);    

    Vec rho = vecInitZerosA(N);
    Vec rho_tmp = vecInitZerosA(N);
    Vec rho_d = vecInitZerosA(N);
    Vec Ec_vec = vecInitZerosA(N);
    Vec res = vecInitZerosA(N);

    Vec subdiag_jacobian = vecInitA(-1, N);
    Vec ref_diag_jacobian = vecInitA(2.0, N);
    Vec diag_jacobian = vecInitA(0.0, N);

    Vec scratch = vecInitZerosA(N);

    PyVi pyvi = pyviInitA("v.pyvi");

    PyViBase pos = pyviCreateParameter(&pyvi, "position", mesh.x);
    PyViSec pyvi_pot = pyviCreateSection(&pyvi, "potential", pos);
    PyViSec pyvi_ch = pyviCreateSection(&pyvi, "charge", pos);
    PyViSec pyvi_chd = pyviCreateSection(&pyvi, "charge derivative", pos);
    PyViSec pyvi_res = pyviCreateSection(&pyvi, "residual", pos);
    PyViSec pyvi_J = pyviCreateSection(&pyvi, "Jacobian", pos);
    PyViSec pyvi_F = pyviCreateSection(&pyvi, "Fermi lvl", pos);

    // setup initial guess
    Vec accum = vecInitZerosA(N);

    pyviSectionPush(pyvi_pot, mesh.potential);
    pyviSectionPush(pyvi_ch, rho);
    pyviSectionPush(pyvi_chd, rho_d);
    pyviSectionPush(pyvi_J, diag_jacobian);
    pyviSectionPush(pyvi_res, accum);

    vecPrint(mesh.dx);
    printf(":len = %zu\n", mesh.len);

    for(size_t i = 0; i < 10; i++)
    {
        Vec ch_scaled = vecInitZerosA(N);

        vecScale(-1, mesh.potential, &Ec_vec);
        meshPoissonEvaluate(mesh, silicon.bulk.epsilon, &accum);

        scTotalQV(silicon, fermi_lvl, Ec_vec, env, &rho);
        vecSub(accum, rho, &accum);

        bulkNetChargeVecDerivative(siliconInfo, dopants, 2, fermi_lvl_vec, Ec_vec, env, &rho_d);
        vecScale(env.spacing*env.spacing / siliconInfo.epsilon, rho_d, &ch_scaled);
        vecAdd(ref_diag_jacobian, ch_scaled, &diag_jacobian);

        solveTridiagonalSymm(accum, subdiag_jacobian, diag_jacobian, scratch);
        vecSub(V, accum, &accum);

        //vecCopy(V, &mesh.potential);
        fxInterpolateSample1D(fxPotential, mesh.x, &mesh.potential);
        meshPoissonEvaluate(mesh, siliconInfo.epsilon, &mesh_test);

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
