#pragma once

#include <math/linalg.h>
#include <math/stack.h>

typedef struct Mesh
{
    // stores the grid points
    Vec x;
    // stores the step sizes
    // x[i] + dx[i] = x[i + 1]
    Vec dx;
    // potential at each point
    Vec potential;
    // the jacobian of poission equation
    // wrt potential
    MatTriDiag jacobianD2;

    // size of mesh
    size_t len;

    // bounds of the mesh
    double min;
    double max;

    // the value of potential at boundaries
    double v0;
    double v1;
} Mesh;

// initialize the grid with piecewise uniform sampling
// assumes ranges are in increasing order
// Eg: ranges[] = { 0.0, 0.4, 0.5, 1.0 }, sampleCount[] = { 5, 10, 5 }
// (0.0, 0.4] => 5 samples
// (0.4, 0.5] => 10 samples
// (0.5, 1.0) => 5 samples**
// count = no of elements in sampleCount
Mesh meshInitPieceUniformA(Vec ranges, size_t sampleCount[], size_t count);

void meshSetDirichletBC(Mesh* mesh, double V0, double V1);

// calculate the first(central) Derivative wrt position for a vector on this mesh
void meshFirstDerivative(Mesh mesh, Vec var, double boundaries[2], Vec* out);

// calculate the second(central) Derivative wrt position for a vector on this mesh
void meshSecondDerivative(Mesh mesh, Vec var, double boundaries[2], Vec* out);

// evaluate the poisson equation
void meshPoissonEvaluate(Mesh mesh, double epsilon, Vec* rho);

// get the jacobian of poission equation
void meshPoissonJacobian(Mesh mesh, double epsilon, MatTriDiag* jacobian);

void freeMesh(Mesh* mesh);
