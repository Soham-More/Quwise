#pragma once

// this library helps with better visualizing C vectors in python

#include <stdio.h>
#include <include/linalg.h>
#include <include/stack.h>

typedef struct PyViParameter
{
    // section name
    const char* name;
    Vec axis;
} PyViParameter;

typedef struct PyViSection
{
    // section name
    const char* name;
    // all the Vec data, in this section
    DynStack/*Vec*/ data;
    PyViParameter parameter;
} PyViSection;

typedef struct PyVi
{
    FILE* file;
    DynStack/*PyViSection*/ sections;
    DynStack/*PyViParameter*/ parameters;
} PyVi;

PyVi pyviInitA(const char* filename);

PyViSection* pyviCreateSection(PyVi* pyvi, const char* section_name, PyViParameter p);

// Copies p by reference. DO NOT FREE p BEFORE PYVI is freed
PyViParameter pyviCreateParameter(PyVi* pyvi, const char* param_name, Vec p);

// push a vector fx varying with parameter x, copies the vector
void pyviSectionPush(PyViSection* section, Vec fx);

// writes all the data to file
void pyviWrite(PyVi pyvi);

// free the PyVi struct
// all related pyvi* data are invalid from 
// this point onwards
void freePyVi(PyVi* pyvi);

