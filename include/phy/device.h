#pragma once

#include <include/utils.h>
#include <include/stack.h>
#include <include/linalg.h>

#include <include/functions.h>
#include <include/material.h>

#include <include/mesh.h>

typedef struct DeviceState
{
    // potential
    Vec V;

    // currents
    Vec Jp;
    Vec Jn;
    Vec Jnet;

    // bandgaps
    Vec Ec;
    Vec Ef;
    Vec Ev;

    // carrier concentrations
    Vec p;
    Vec n;
} DeviceState;

typedef struct Device
{
    Bulk bulk;
    DynStack dopants;
    Mesh mesh;
    // reference potential
    double vref;
    // fermi-level at x = 0
    double fermi_x0;
} Device;

;


