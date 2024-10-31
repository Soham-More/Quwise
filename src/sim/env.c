#include "include/sim/env.h"
#include <include/utils.h>

// calculate the rest of env members
Environment envPopulate(Environment env)
{
    env.thermal_potential = env.temperature * BOLTZMANN_CONSTANT_EV;
    return env;
}