#pragma once

// this library helps with better visualizing C vectors in python

#include <stdio.h>
#include <include/linalg.h>
#include <include/stack.h>

typedef struct PyViBase
{
    const char* name;
    Vec axis;
} PyViBase;

typedef struct PyVi
{
    FILE* file;
    DynStack/*PyViSection*/ sections;
    DynStack/*PyViParameter*/ parameters;
} PyVi;

typedef struct PyViSec
{
    size_t id;
    PyVi* pyvi;
} PyViSec;

PyVi pyviInitA(const char* filename);

PyViSec pyviCreateSection(PyVi* pyvi, const char* section_name, PyViBase p);

// Copies p by reference. DO NOT FREE p BEFORE PYVI is freed
PyViBase pyviCreateParameter(PyVi* pyvi, const char* param_name, Vec p);

// push a vector fx varying with parameter x, copies the vector
void pyviSectionPush(PyViSec section, Vec fx);

// writes all the data to file
void pyviWrite(PyVi pyvi);

// free the PyVi struct
// all related pyvi* data are invalid from 
// this point onwards
void freePyVi(PyVi* pyvi);

