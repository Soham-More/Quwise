#include <include/linalg.h>

// Solves poisson equation(exactly) for uniform sampling
// assumes x and rho_i have same dimension
void poisson1DSolve(Vec* x, Vec rho_i, const double epsilon, const double left_boundary, const double right_boundary, const double x_step);

// assumes V and rho_i have same dimension
void poissonEvaluate(Vec V, Vec* rho, const double left_boundary, const double right_boundary, const double epsilon, const double x_step);


void solveTridiagonalSymm(Vec x, Vec a, Vec b, Vec scratch);
void toPoissonChargeVec(Vec rho, const double epsilon, const double left_boundary, const double right_boundary, const double x_step, Vec* poisson_charge_vec);
