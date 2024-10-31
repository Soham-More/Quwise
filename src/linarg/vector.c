#include "include/linalg.h"

#include <stdlib.h>
#include <memory.h>
#include <math.h>

// gets value at index from vector by reference(dereferenced)
// allows for syntax like: LA_VIDX(a, 2) = 5;
// DOES NOT CHECK FOR OUT OF BOUNDS ACCESS
#define LA_VIDX(vector, index) *(vector.x + vector.offset * index)

// gets value at index from vector(ptr) by reference(dereferenced)
// allows for syntax like: LA_VIDX(a, 2) = 5;
// DOES NOT CHECK FOR OUT OF BOUNDS ACCESS
#define LA_VIDX_PTR(vector, index) *(vector->x + vector->offset * index)

// initialzie the vector on the heap with some initial value
Vec vecInitA(double value, size_t len)
{
    if(len == 0)
    {
        LINALG_REPORT_ERROR("invalid zero length vector requested!");
        return (Vec){ NULL, 0, 0 };
    }
    Vec x = { (double*)calloc(len, sizeof(double)), len, 1 };
    LINALG_ASSERT_ERROR(!x.x, x, "unkown error occured when allocation memory!");
    for(size_t i = 0; i < x.len; i++) x.x[i] = value;
    return x;
}
// initialize the vector on the heap to zeros
Vec vecInitZerosA(size_t len)
{
    return vecInitA(0.0, len);
}
// initialize the vector on the heap to ones
Vec vecInitOnesA(size_t len)
{
    return vecInitA(1.0, len);
}

// make a copy of a vector on heap
Vec vecCopyA(Vec vector)
{
    if(vector.len == 0)
    {
        LINALG_REPORT_ERROR("invalid zero length vector requested!");
        return (Vec){ NULL, 0, 0 };
    }
    else if(vector.x == NULL)
    {
        LINALG_REPORT_ERROR("invalid source pointer(null)!");
        return (Vec){ NULL, 0, 0 };
    }
    Vec x = { (double*)calloc(vector.len, sizeof(double)), vector.len, 1 };
    LINALG_ASSERT_ERROR(!x.x, x, "unkown error occured when allocation memory!");
    for(size_t i = 0; i < x.len; i++) x.x[i] = LA_VIDX(vector, i);
    return x;
}
// make a copy of a vector on another vector
int vecCopy(Vec src, Vec* dst)
{
    LINALG_ASSERT_ERROR(src.len == 0, LINALG_ERROR, "source length is zero!");
    LINALG_ASSERT_ERROR(src.len != dst->len, LINALG_ERROR, "attempt to copy vectors with unqeual dimensions %zu to %zu!", src.len, dst->len);
    for(size_t i = 0; i < src.len; i++) LA_VIDX_PTR(dst, i) = LA_VIDX(src, i);
    return LINALG_OK;
}

// construct a vector from a pointer(does not allocate)
Vec vecConstruct(double* ptr, size_t len)
{
    Vec x = { ptr, len, 1 };
    return x;
}

// pretty print a vector
void vecPrint(Vec a)
{
    printf("[");
    for(size_t i = 0; i < a.len - 1; i++)
    {
        printf("%le, ", LA_VIDX(a, i));
    }
    printf("%le]", LA_VIDX(a, a.len - 1));
}

// gets the nth value in a vector
// handles buffer offsets
// checks for out-of-bounds
double vecGet(Vec a, size_t n)
{
    LINALG_ASSERT_ERROR(n >= a.len, NAN, "out of bounds vector access!");
    return *(a.x + a.offset * n);
}

// gets the nth value in a vector
// handles buffer offsets
// checks for out-of-bounds
double* vecRef(Vec a, size_t n)
{
    LINALG_ASSERT_ERROR(!a.x, NULL, "null vector!");
    LINALG_ASSERT_ERROR(n >= a.len, NULL, "out of bounds vector access!");
    return (a.x + a.offset * n);
}

// add 2 vectors and get result into another vector
int vecAdd(Vec a, Vec b, Vec* result)
{
    LINALG_ASSERT_ERROR(a.len != b.len, LINALG_ERROR, "attempt to add vectors with dimension %zu and %zu!", a.len, b.len);
    LINALG_ASSERT_ERROR(!result || !result->x, LINALG_ERROR, "result vector is null!");
    LINALG_ASSERT_ERROR(!a.x || !b.x, LINALG_ERROR, "input vector/s is/are null!");
    LINALG_ASSERT_ERROR(b.len < result->len, LINALG_ERROR, "output vector not big enough to store result!");
    LINALG_ASSERT_ERROR(b.len > result->len, LINALG_ERROR, "output dimension larger than input dimension!");

    for(size_t i = 0; i < a.len; i++)
    {
        LA_VIDX_PTR(result, i) = LA_VIDX(a, i) + LA_VIDX(b, i);
    }

    return LINALG_OK;
}
// subtract 2 vectors(a - b) and get result into another vector, prints error if input is invalid
int vecSub(Vec a, Vec b, Vec* result)
{
    LINALG_ASSERT_ERROR(a.len != b.len, LINALG_ERROR, "attempt to add vectors with dimension %zu and %zu!", a.len, b.len);
    LINALG_ASSERT_ERROR(!result || !result->x, LINALG_ERROR, "result vector is null!");
    LINALG_ASSERT_ERROR(!a.x || !b.x, LINALG_ERROR, "input vector/s is/are null!");
    LINALG_ASSERT_ERROR(b.len < result->len, LINALG_ERROR, "output vector not big enough to store result!");
    LINALG_ASSERT_ERROR(b.len > result->len, LINALG_ERROR, "output dimension larger than input dimension!");

    for(size_t i = 0; i < a.len; i++)
    {
        LA_VIDX_PTR(result, i) = LA_VIDX(a, i) - LA_VIDX(b, i);
    }

    return LINALG_OK;
}
// multiply scalar value to vectors and get result into another vector
int vecScale(double a, Vec b, Vec* result)
{
    LINALG_ASSERT_ERROR(!result || !result->x, LINALG_ERROR, "resultant vector is null!");
    LINALG_ASSERT_ERROR(!b.x, LINALG_ERROR, "input vector/s is/are null!");
    LINALG_ASSERT_ERROR(b.len < result->len, LINALG_ERROR, "output vector not big enough to store result!");
    LINALG_ASSERT_ERROR(b.len > result->len, LINALG_ERROR, "output dimension larger than input dimension!");

    for(size_t i = 0; i < b.len; i++)
    {
        LA_VIDX_PTR(result, i) = a * LA_VIDX(b, i);
    }

    return LINALG_OK;
}
// divide a scalar value to vector and get result into another vector, prints error if input is invalid
int vecRScale(double a, Vec b, Vec* result)
{
    LINALG_ASSERT_ERROR(!result || !result->x, LINALG_ERROR, "resultant vector is null!");
    LINALG_ASSERT_ERROR(!b.x, LINALG_ERROR, "input vector/s is/are null!");
    LINALG_ASSERT_ERROR(b.len < result->len, LINALG_ERROR, "output vector not big enough to store result!");
    LINALG_ASSERT_ERROR(b.len > result->len, LINALG_ERROR, "output dimension larger than input dimension!");

    for(size_t i = 0; i < b.len; i++)
    {
        LA_VIDX_PTR(result, i) = LA_VIDX(b, i) / a;
    }

    return LINALG_OK;
}
// calculate exp of every component in vector and get result into another vector, prints error if input is invalid
int vecExp(Vec b, Vec* result)
{
    LINALG_ASSERT_ERROR(!result || !result->x, LINALG_ERROR, "resultant vector is null!");
    LINALG_ASSERT_ERROR(!b.x, LINALG_ERROR, "input vector/s is/are null!");
    LINALG_ASSERT_ERROR(b.len < result->len, LINALG_ERROR, "output vector not big enough to store result!");
    LINALG_ASSERT_ERROR(b.len > result->len, LINALG_ERROR, "output dimension larger than input dimension!");

    for(size_t i = 0; i < b.len; i++)
    {
        LA_VIDX_PTR(result, i) = exp(LA_VIDX(b, i));
    }

    return LINALG_OK;
}
// unit vector of the norm
int vecNormalize(Vec a, Vec* result)
{
    return vecScale(1 / vecMagnitude(a), a, result);
}
// get the dot product between 2 variables
double vecDot(Vec a, Vec b)
{
    LINALG_ASSERT_ERROR(a.len != b.len, NAN, "attempt to take dot product of vectors with dimension %zu and %zu!", a.len, b.len);
    LINALG_ASSERT_ERROR(!a.x || !b.x, NAN, "input vector/s is/are null!");

    double result = 0.0;

    for(size_t i = 0; i < a.len; i++)
    {
        result += LA_VIDX(a, i) * LA_VIDX(b, i);
    }

    return result;
}
// get the L2 norm of vector
double vecMagnitude(Vec a)
{
    return sqrt(vecDot(a, a));
}
// get the L_p norm of vector, prints warning if p < 1
// for p = inf use vecMax 
double vecNorm(Vec a, double p)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_ERROR(p < 1, NAN, "L_p is not a valid norm for p = %f", p);

    double result = 0.0;

    for(size_t i = 0; i < a.len; i++)
    {
        result +=  pow(fabs(LA_VIDX(a, i)), p);
    }

    return pow(result, 1/p);
}
// maximum value in the vector
double vecMax(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, INFINITY, "max of a zero dimension vector");

    double result = -INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        result = result > LA_VIDX(a, i) ? result : LA_VIDX(a, i);
    }

    return result;
}
// maximum value(abs) in the vector
double vecMaxAbs(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, INFINITY, "max of a zero dimension vector");

    double result = -INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        result = result > fabs(LA_VIDX(a, i)) ? result : fabs(LA_VIDX(a, i));
    }

    return result;
}
// minimum value in the vector
double vecMin(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, -INFINITY, "min of a zero dimension vector");

    double result = INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        result = result < LA_VIDX(a, i) ? result : LA_VIDX(a, i);
    }

    return result;
}
// minimum value(abs) in the vector
double vecMinAbs(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, -INFINITY, "min of a zero dimension vector");

    double result = INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        result = result < fabs(LA_VIDX(a, i)) ? result : fabs(LA_VIDX(a, i));
    }

    return result;
}
// sum all values in a vector
double vecSum(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, 0, "sum of a zero dimension vector");

    double result = 0;
    for(size_t i = 0; i < a.len; i++) result += LA_VIDX(a, i);

    return result;
}
// return the product of all values in a vector
double vecProd(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, 1, "product of a zero dimension vector");

    double result = 1;
    for(size_t i = 0; i < a.len; i++) result *= LA_VIDX(a, i);

    return result;
}
// get the range of vector, i.e max - min
double vecRange(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, 0, "input is a zero dimension vector");

    double max_val = -INFINITY;
    double min_val = INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        max_val = max_val > LA_VIDX(a, i) ? max_val : LA_VIDX(a, i);
        min_val = min_val < LA_VIDX(a, i) ? min_val : LA_VIDX(a, i);
    }

    return max_val - min_val;
}
// get the (relative) range of vector, i.e (max - min) / min( |max|, |min| )
double vecRangeRelative(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, 0, "input is a zero dimension vector");

    double max_val = -INFINITY;
    double min_val = INFINITY;
    double abs_min = INFINITY;

    for(size_t i = 0; i < a.len; i++)
    {
        max_val = max_val > LA_VIDX(a, i) ? max_val : LA_VIDX(a, i);
        min_val = min_val < LA_VIDX(a, i) ? min_val : LA_VIDX(a, i);
        abs_min = abs_min < fabs(LA_VIDX(a, i)) ? min_val : fabs(LA_VIDX(a, i));
    }

    return (max_val - min_val) / abs_min;
}
// get the standard deviation of the vector
double vecStandardDeviation(Vec a)
{
    LINALG_ASSERT_ERROR(!a.x, NAN, "input vector is null!");
    LINALG_ASSERT_WARN(a.len == 0, 0, "checking a zero dimension vector");

    double mean = vecSum(a)/(double)a.len;
    double dev = 0;
    for(size_t i = 0; i < a.len; i++) dev += fabs(LA_VIDX(a, i) - mean);

    return dev / (double)a.len;
}

//swap two vectors
void swapVec(Vec* a, Vec* b)
{
    Vec imm = *a;
    *a = *b;
    *b = imm;
}

// free the vector on the heap
void freeVec(Vec* vec)
{
    // skip if null
    if(!vec->x) return;

    free(vec->x);
    vec->x = NULL;
    vec->len = 0;
}