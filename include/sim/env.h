#pragma once

#include <math/linalg.h>

typedef struct Environment
{
    double tol_potential;
    double tol_conc;
    double temperature;
    double thermal_potential;
} Environment;

// calculate the rest of env members
Environment envPopulate(Environment env);


