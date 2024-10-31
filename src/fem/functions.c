#include <include/fem/functions.h>
#include <math.h>

double interp_linear(Vec values, Vec x, double x_eval)
{
    size_t nearest_low = 0;
    double nearest_low_x = -INFINITY;

    size_t nearest_high = 0;
    double nearest_high_x = INFINITY;

    for(size_t i = 0; i < x.len; i++)
    {
        double x_sample = VEC_INDEX(x, i);

        if(nearest_low_x < x_sample && x_sample < x_eval)
        {
            nearest_low_x = x_sample;
            nearest_low = i;
        }
        else if(nearest_high_x > x_sample && x_sample > x_eval)
        {
            nearest_high_x = x_sample;
            nearest_high = i;
        }
    }

    double m = (VEC_INDEX(values, nearest_high) - VEC_INDEX(values, nearest_low)) / (nearest_high_x - nearest_low_x);

    return VEC_INDEX(values, nearest_low) + m * ( x_eval - nearest_low_x );
}

double interp_nearest(Vec values, Vec x, double x_eval)
{
    double smallest_dist = INFINITY;
    size_t smallest_index = 0;

    for(size_t i = 0; i < x.len; i++)
    {
        double sample_dist = fabs(VEC_INDEX(x, i) - x_eval);

        if(sample_dist < smallest_dist)
        {
            smallest_dist = sample_dist;
            smallest_index = i;
        }
    }

    return VEC_INDEX(values, smallest_index);
}

Fx1D fxConstruct1D(Vec values, Vec x, uint8_t flags)
{
    Fx1D fx1D;
    fx1D.values = values;
    fx1D.x = x;
    fx1D.interp_flags = flags;

    return fx1D;
}

// evaluate function at x, using interpolation specified in Fx1D
double fxInterpolateSingle1D(Fx1D fx, double x)
{
    switch (fx.interp_flags)
    {
        case INTERP_NEAREST:
            return interp_nearest(fx.values, fx.x, x);
            break;
        case INTERP_LINEAR:
            return interp_linear(fx.values, fx.x, x);
            break;
        default:
            break;
    }
    printf("Warning: no valid interpolation type selected!\n");
    return NAN;
}

// evaluate function at all x, using interpolation specified in Fx1D
void fxInterpolateSample1D(Fx1D fx, Vec x, Vec* nvalues)
{
    for(size_t i = 0; i < x.len; i++)
    {
        VEC_INDEX(*nvalues, i) = fxInterpolateSingle1D(fx, VEC_INDEX(x, i));
    }
}


