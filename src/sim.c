#include "include/sim.h"
#include <include/utils.h>

SimParameters simConstructA(size_t sample_count, double start_x, double end_x, double temp, double tol_poisson, double tol_conc)
{
    SimParameters simParams;

    simParams.sample_count = sample_count;
    simParams.start_x = start_x;
    simParams.end_x = end_x;
    simParams.temperature = temp;
    simParams.thermal_potential = temp * BOLTZMANN_CONSTANT_EV;
    simParams.tol_conc = tol_conc;
    simParams.tol_potential = tol_poisson;

    simParams.x_samples = vecInitZerosA(sample_count);
    simParams.spacing = (end_x - start_x) / (double)(sample_count + 2);

    for(size_t i = 0; i < sample_count; i++)
    {
        VEC_INDEX(simParams.x_samples, i) = (end_x - start_x) * ((double)(i + 1) / (double)(sample_count + 2));
    }

    return simParams;
}
void freeSim(SimParameters* simParams)
{
    freeVec(&simParams->x_samples);
}

