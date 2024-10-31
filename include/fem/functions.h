#pragma once

#include <math/linalg.h>

#include <stdint.h>

#define INTERP_NONE 0x0
#define INTERP_NEAREST 0x1
#define INTERP_LINEAR 0x2
#define INTERP_LAGRANGE 0x4

// represents a sampled function(1D)
typedef struct Fx1D
{
    Vec values;
    Vec x;
    uint8_t interp_flags;
} Fx1D;

// this keeps the pointer to values and x
Fx1D fxConstruct1D(Vec values, Vec x, uint8_t flags);

// evaluate function at x, using interpolation specified in Fx1D
double fxInterpolateSingle1D(Fx1D fx, double x);

// evaluate function at all x, using interpolation specified in Fx1D
void fxInterpolateSample1D(Fx1D fx, Vec x, Vec* nvalues);

