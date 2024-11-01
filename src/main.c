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
    /*
    Mat2d scratch = mat2DInitZerosA(5, 6);
    Mat2d matrix = mat2DConstruct((double[]){
         2.0,-1.0, 0.0, 0.0, 0.0,
        -3.0, 2.0,-1.0, 0.0, 0.0,
         0.0,-1.0, 2.0,-1.0, 0.0,
         0.0, 0.0,-1.0, 2.0,-1.0,
         0.0, 0.0, 0.0,-1.0, 2.0,
    }, 5, 5);

    size_t arr[] = {0.0, 0.0, 0.0, 0.0, 0.0};

    Vec ones = vecInitOnesA(5);
    Vec ans  = vecInitOnesA(5);

    mat2DSqSolve(matrix, ones, &scratch, arr, &ans);

    vecPrint(ans);
    printf("\n");

    freeMat2D(&scratch);
    */

    size_t N = 1 << 6;

    double length = 0.2e-5;
    double n_doping = 1e21;
    double p_doping = 1e21;

    // make a simulation setup
    Environment env;
    env.temperature = 300;
    env.tol_potential = 1e-6;
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
    boronInfo.degeneracy = 4;
    boronInfo.samplingMode = INTERP_NEAREST;

    SemiConductor silicon = scInit(siliconInfo, env);

    double v_estimate = env.thermal_potential * (10 * log(10));

    double xn = length/2.0 + sqrt(2 * silicon.bulk.epsilon / ELECTRON_CHARGE * (p_doping / n_doping) /( p_doping + n_doping ) * v_estimate);
    double xp = length/2.0 - sqrt(2 * silicon.bulk.epsilon / ELECTRON_CHARGE * (n_doping / p_doping) /( p_doping + n_doping ) * v_estimate);

    printf("x_n: %le\n", xn - length/2.0);
    printf("x_p: %le\n", xp - length/2.0);


    Mesh mesh = meshInitPieceUniformA(
        vecConstruct((double[]){0.0, length / 2 - 0.1e-5, length / 2 + 0.1e-5, length}, 4),
        (size_t[]){N / 4, N / 2 + 1, N / 4 - 1},
        3
    );

/*
    Mesh mesh = meshInitPieceUniformA(
        vecConstruct((double[]){0.0, length}, 2),
        (size_t[]){N},
        1
    );
*/
    scDopeAcceptor(&silicon, mesh, boronInfo);
    scDopeDonor(&silicon, mesh, phosphorousInfo);

    // solve for actual fermi and conduction level
    double boundary_ef = scSolveBoundary(silicon, 0.0, 0.0, SC_SOLVE_FERMI);
    double boundary_ec = scSolveBoundary(silicon, length, boundary_ef, SC_SOLVE_EC);

    printf("Fermi-Level: %lf eV\n", boundary_ef);
    printf("Built-in Potential: %lf V\n", -boundary_ec);

    meshSetDirichletBC(&mesh, 0.0, -boundary_ec);

    Vec fermi_lvl = vecInitA(boundary_ef, mesh.len);    

    Vec rho = vecInitZerosA(N);
    Vec rho_tmp = vecInitZerosA(N);
    Vec rho_d = vecInitZerosA(N);
    Vec Ec_vec = vecInitZerosA(N);
    Vec res = vecInitZerosA(N);

    MatTriDiag jacobian = triDiagInitZeroA(N);
    Vec accum = vecInitZerosA(N);

    Vec scratch = vecInitZerosA(N);

    PyVi pyvi = pyviInitA("v.pyvi");

    PyViBase pos = pyviCreateParameter(&pyvi, "position", mesh.x);
    PyViSec pyvi_pot = pyviCreateSection(&pyvi, "potential", pos);
    PyViSec pyvi_ch = pyviCreateSection(&pyvi, "charge", pos);
    PyViSec pyvi_chd = pyviCreateSection(&pyvi, "charge derivative", pos);
    PyViSec pyvi_res = pyviCreateSection(&pyvi, "residual", pos);
    PyViSec pyvi_J = pyviCreateSection(&pyvi, "Jacobian", pos);
    //PyViSec pyvi_F = pyviCreateSection(&pyvi, "Fermi lvl", pos);

    pyviSectionPush(pyvi_pot, mesh.potential);
    pyviSectionPush(pyvi_ch, rho);
    pyviSectionPush(pyvi_chd, rho_d);
    pyviSectionPush(pyvi_J, jacobian.superdiagonal);
    pyviSectionPush(pyvi_res, accum);

    vecPrint(mesh.dx);
    printf(":len = %zu\n", mesh.len);

    for(size_t i = 0; i < 20; i++)
    {
        Vec ch_scaled = vecInitZerosA(N);

        vecScale(-1, mesh.potential, &Ec_vec);
        meshPoissonEvaluate(mesh, silicon.bulk.epsilon, &accum);

        scVecTotalQ(silicon, fermi_lvl, Ec_vec, &rho);
        vecSub(accum, rho, &accum);

        meshPoissonJacobian(mesh, silicon.bulk.epsilon, &jacobian);
        scVecTotalQD(silicon, fermi_lvl, Ec_vec, &rho_d);
        triDiagAddDiagonalSelf(&jacobian, rho_d);

        pyviSectionPush(pyvi_J, mesh.dx);

        triDiagSolveDestructive(&jacobian, &accum);
        vecSub(mesh.potential, accum, &accum);

        pyviSectionPush(pyvi_pot, mesh.potential);
        pyviSectionPush(pyvi_ch, rho);
        pyviSectionPush(pyvi_chd, rho_d);
        //pyviSectionPush(pyvi_J, diag_jacobian);
        pyviSectionPush(pyvi_res, accum);

        // set V to new value
        vecCopy(accum, &mesh.potential);
    }

    pyviWrite(pyvi);

    freePyVi(&pyvi);
    freeMesh(&mesh);
}
