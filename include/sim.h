#pragma once

#include <math/linalg.h>

typedef struct SimParameters
{
    size_t sample_count;
    double start_x;
    double end_x;

    // samples on the open interval, (start_x, end_x)
    Vec x_samples;
    double spacing;

    double tol_potential;
    double tol_conc;

    double temperature;
    double thermal_potential;
} SimParameters;

SimParameters simConstructA(size_t sample_count, double start_x, double end_x, double temp, double tol_poisson, double tol_conc);
void freeSim(SimParameters* simParams);


