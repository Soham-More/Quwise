#include <math/linalg.h>
#include <include/poisson.h>

// Solves poisson equation(exactly) for uniform sampling
// assumes x and rho_i have same dimension
void poisson1DSolve(Vec* x, Vec rho_i, const double epsilon, const double left_boundary, const double right_boundary, const double x_step)
{
    double csum = left_boundary;

    for(size_t i = 0; i < rho_i.len; i++)
    {
        double c_i = 1.0 - 1.0 / (i + 2);
        if (i == rho_i.len - 1) csum += (double)rho_i.len * right_boundary;

        csum += (double)(i + 1) * VEC_INDEX(rho_i, i) * (x_step * x_step) / epsilon;
        VEC_INDEX(*x, i) = csum * c_i / ( (double)(i + 1)*(double)(i + 1) );
    }
    //printf("Dbg:\n");
    //vecPrint(*x);
    //printf("\n");
    csum = 0;
    for(size_t i = rho_i.len - 1; 1; i--)
    {
        csum += VEC_INDEX(*x, i);
        VEC_INDEX(*x, i) = csum * (double)(i + 1);
        if(i == 0) break;
    }
}

// assumes V and rho_i have same dimension
void poissonEvaluate(Vec V, Vec* rho, const double left_boundary, const double right_boundary, const double epsilon, const double x_step)
{
    for(size_t i = 0; i < V.len; i++)
    {
        VEC_INDEX(*rho, i) = 2.0 * VEC_INDEX(V, i);

        if(i > 0) VEC_INDEX(*rho, i) -= VEC_INDEX(V, i - 1);
        else VEC_INDEX(*rho, i) -= left_boundary;

        if(i < V.len - 1) VEC_INDEX(*rho, i) -= VEC_INDEX(V, i + 1);
        else VEC_INDEX(*rho, i) -= right_boundary;
    }
}

void solveTridiagonalSymm(Vec x, Vec a, Vec b, Vec scratch)
{
    VEC_INDEX(scratch, 0) = VEC_INDEX(a, 0) / VEC_INDEX(b, 0);
    VEC_INDEX(x, 0) = VEC_INDEX(x, 0) / VEC_INDEX(b, 0);

    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < b.len; ix++)
    {
        if (ix < b.len-1)
        {
            VEC_INDEX(scratch, ix) = VEC_INDEX(a, ix) / (VEC_INDEX(b, ix) - VEC_INDEX(a, ix) * VEC_INDEX(scratch, ix - 1));
        }
        VEC_INDEX(x, ix) = (VEC_INDEX(x, ix) - VEC_INDEX(a, ix) * VEC_INDEX(x, ix - 1)) / (VEC_INDEX(b, ix) - VEC_INDEX(a, ix) * VEC_INDEX(scratch, ix-1));
    }

    /* loop from X - 2 to 0 inclusive */
    for (size_t ix = b.len - 2; ix > 0; ix--)
    {
        VEC_INDEX(x, ix) -= VEC_INDEX(scratch, ix) * VEC_INDEX(x, ix + 1);
    }
    VEC_INDEX(x, 0) -= VEC_INDEX(scratch, 0) * VEC_INDEX(x, 1);
}

void toPoissonChargeVec(Vec rho, const double epsilon, const double left_boundary, const double right_boundary, const double x_step, Vec* poisson_charge_vec)
{
    vecScale((x_step*x_step) / epsilon, rho, poisson_charge_vec);
    VEC_INDEX(*poisson_charge_vec, 0) += left_boundary;
    VEC_INDEX(*poisson_charge_vec, poisson_charge_vec->len - 1) += right_boundary;
}

