#include <include/fem/mesh.h>

// initialize the grid with piecewise uniform sampling
Mesh meshInitPieceUniformA(Vec ranges, size_t sampleCount[], size_t count)
{
    Mesh mesh;
    mesh.len = 0;

    for(size_t i = 0; i < count; i++) mesh.len += sampleCount[i];

    mesh.x = vecInitZerosA(mesh.len);
    mesh.dx = vecInitZerosA(mesh.len);
    mesh.potential = vecInitZerosA(mesh.len);
    mesh.min = VEC_INDEX(ranges, 0);
    mesh.max = VEC_INDEX(ranges, count);

    // populate the mesh
    size_t point_count = 0;
    for(size_t i = 0; i < count - 1; i++)
    {
        double interval = (VEC_INDEX(ranges, i + 1) - VEC_INDEX(ranges, i));
        double spacing = interval / (double)(sampleCount[i]);

        for(size_t j = 0; j < sampleCount[i]; j++)
        {
            VEC_INDEX(mesh.x, point_count) = VEC_INDEX(ranges, i) + interval * ((double)(j + 1) / (double)(sampleCount[i]));
            if(point_count > 0) VEC_INDEX(mesh.dx, point_count - 1) = spacing;
            point_count++;
        }
    }
    double interval = (VEC_INDEX(ranges, count) - VEC_INDEX(ranges, count - 1));
    double spacing = interval / (double)(sampleCount[count - 1] + 1);
    double samples = (double)(sampleCount[count - 1] + 1);

    for(size_t i = 0; i < sampleCount[count - 1]; i++)
    {
        VEC_INDEX(mesh.x, point_count) = VEC_INDEX(ranges, count - 1) + interval * ((double)(i + 1) / samples);
        if(point_count > 0) VEC_INDEX(mesh.dx, point_count - 1) = spacing;
        point_count++;
    }
    if(point_count > 0) VEC_INDEX(mesh.dx, point_count - 1) = spacing;

    /*
    for(size_t i = 0; i < mesh.len - 1; i++)
    {
        VEC_INDEX(mesh.dx, i) = VEC_INDEX(mesh.x, i + 1) - VEC_INDEX(mesh.x, i);
    }
    VEC_INDEX(mesh.dx, mesh.len - 1) = mesh.max - VEC_INDEX(mesh.x, mesh.len - 1);
    */

    mesh.jacobianD2 = triDiagInitZeroA(mesh.len);

    // Calculate the jacobian of poission equation
    double h = (VEC_INDEX(mesh.dx, 0) - mesh.min) + VEC_INDEX(mesh.dx, 0);
    double h_eff = (VEC_INDEX(mesh.dx, 0) - mesh.min) * VEC_INDEX(mesh.dx, 0);
    VEC_INDEX(mesh.jacobianD2.diagonal, 0) = -2.0 / h_eff;
    VEC_INDEX(mesh.jacobianD2.superdiagonal, 0) = 2.0 / ( VEC_INDEX(mesh.dx, 0) * h );
    
    for(size_t i = 1; i < mesh.len - 1; i++)
    {
        h = VEC_INDEX(mesh.dx, i - 1) + VEC_INDEX(mesh.dx, i);
        h_eff = VEC_INDEX(mesh.dx, i - 1) * VEC_INDEX(mesh.dx, i);
        VEC_INDEX(mesh.jacobianD2.subdiagonal, i - 1) = 2.0 / ( VEC_INDEX(mesh.dx, i - 1) * h );
        VEC_INDEX(mesh.jacobianD2.diagonal, i) = -2.0 /  h_eff;
        VEC_INDEX(mesh.jacobianD2.superdiagonal, i) = 2.0 / ( VEC_INDEX(mesh.dx, i) * h );
    }

    h = VEC_INDEX(mesh.dx, mesh.len - 2) + VEC_INDEX(mesh.dx, mesh.len - 1);
    h_eff = VEC_INDEX(mesh.dx, mesh.len - 2) * VEC_INDEX(mesh.dx, mesh.len - 1);
    VEC_INDEX(mesh.jacobianD2.subdiagonal, mesh.len - 2) = 2.0 / ( VEC_INDEX(mesh.dx, mesh.len - 2) * h );
    VEC_INDEX(mesh.jacobianD2.diagonal, mesh.len - 1) = -2.0 /  h_eff;

    return mesh;
}

void meshSetDirichletBC(Mesh* mesh, double V0, double V1)
{
    mesh->v0 = V0;
    mesh->v1 = V1;
}

// calculate the first Derivative wrt position for a vector on this mesh
void meshFirstDerivative(Mesh mesh, Vec var, double boundaries[2], Vec* out)
{
    double r_b = (VEC_INDEX(mesh.x, 0) - mesh.min) / VEC_INDEX(mesh.dx, 0);
    double r_f = VEC_INDEX(mesh.dx, 0) / (VEC_INDEX(mesh.x, 0) - mesh.min);
    double h = VEC_INDEX(mesh.dx, 0) + (VEC_INDEX(mesh.x, 0) - mesh.min);

    VEC_INDEX(*out, 0) = (r_b * VEC_INDEX(var, 1) - r_f * boundaries[0] + (r_b - r_f) * VEC_INDEX(var, 0)) / h;

    for(size_t i = 1; i < mesh.len - 1; i++)
    {
        r_b = VEC_INDEX(mesh.dx, i - 1) / VEC_INDEX(mesh.dx, i);
        r_f = VEC_INDEX(mesh.dx, i) / VEC_INDEX(mesh.dx, i - 1);
        h = VEC_INDEX(mesh.dx, i) + VEC_INDEX(mesh.dx, i - 1);

        VEC_INDEX(*out, 0) = r_b * VEC_INDEX(var, i + 1) - r_f * VEC_INDEX(var, i - 1) + (r_b - r_f) * VEC_INDEX(var, i);
    }

    r_b = VEC_INDEX(mesh.dx, mesh.len - 2) / VEC_INDEX(mesh.dx, mesh.len - 1);
    r_f = VEC_INDEX(mesh.dx, mesh.len - 1) / VEC_INDEX(mesh.dx, mesh.len - 2);
    h = VEC_INDEX(mesh.dx, mesh.len - 2) + VEC_INDEX(mesh.dx, mesh.len - 1);

    VEC_INDEX(*out, mesh.len - 1) = (r_b * boundaries[1] - r_f * VEC_INDEX(var, mesh.len - 2) + (r_b - r_f) * VEC_INDEX(var, mesh.len - 1)) / h;
}

// calculate the second Derivative wrt position for a vector on this mesh
void meshSecondDerivative(Mesh mesh, Vec var, double boundaries[2], Vec* out)
{
    double h = VEC_INDEX(mesh.dx, 0) + (VEC_INDEX(mesh.x, 0) - mesh.min);
    double r_b = 2 * (VEC_INDEX(mesh.x, 0) - mesh.min) / h;
    double r_f = 2 * VEC_INDEX(mesh.dx, 0) / h;
    double h_eff = VEC_INDEX(mesh.dx, 0) * (VEC_INDEX(mesh.x, 0) - mesh.min);

    VEC_INDEX(*out, 0) = (r_b * VEC_INDEX(var, 1) + r_f * boundaries[0] - 2.0 * VEC_INDEX(var, 0)) / h_eff;

    for(size_t i = 1; i < mesh.len - 1; i++)
    {
        h = VEC_INDEX(mesh.dx, i) + VEC_INDEX(mesh.dx, i - 1);
        r_b = 2 * VEC_INDEX(mesh.dx, i - 1) / h;
        r_f = 2 * VEC_INDEX(mesh.dx, i) / h;
        h_eff = VEC_INDEX(mesh.dx, i) * VEC_INDEX(mesh.dx, i - 1);

        VEC_INDEX(*out, i) = (r_b * VEC_INDEX(var, i + 1) + r_f * VEC_INDEX(var, i - 1) - 2.0 * VEC_INDEX(var, i)) / h_eff;
    }

    h = VEC_INDEX(mesh.dx, mesh.len - 1) + VEC_INDEX(mesh.dx, mesh.len - 2);
    r_b = 2 * VEC_INDEX(mesh.dx, mesh.len - 2) / h;
    r_f = 2 * VEC_INDEX(mesh.dx, mesh.len - 1) / h;
    h_eff = VEC_INDEX(mesh.dx, mesh.len - 1) * VEC_INDEX(mesh.dx, mesh.len - 2);

    VEC_INDEX(*out, mesh.len - 1) = (r_b * boundaries[1] + r_f * VEC_INDEX(var, mesh.len - 2) - 2.0 * VEC_INDEX(var, mesh.len - 1)) / h_eff;
}

// evaluate the poisson equation
void meshPoissonEvaluate(Mesh mesh, double epsilon, Vec* rho)
{
    meshSecondDerivative(mesh, mesh.potential, (double[2]){mesh.v0, mesh.v1}, rho);
    vecScale(-epsilon, *rho, rho);
}

// get the jacobian of poission equation
void meshPoissonJacobian(Mesh mesh, double epsilon, MatTriDiag* jacobian)
{
    triDiagScale(-epsilon, mesh.jacobianD2, jacobian);
}

void freeMesh(Mesh* mesh)
{
    freeVec(&mesh->x);
    freeVec(&mesh->dx);
    freeVec(&mesh->potential);
}
