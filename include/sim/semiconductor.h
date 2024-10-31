#pragma once

#include <stdint.h>
#include <math/stack.h>

#include "bulk.h"
#include <include/fem/mesh.h>

#define SC_SOLVE_FERMI 0x1
#define SC_SOLVE_EC 0x2

typedef struct SemiConductor
{
    Bulk bulk;
    DynStack dopants;
    Environment env;
} SemiConductor;

typedef struct BulkInfo
{
    double bandgap;
    double relativePermittivity;
    double relativeEffElectronMass;
    double relativeEffHoleMass;
} BulkInfo;

typedef struct DopantInfo
{
    Vec sampled_conc;
    Vec sampled_x;
    uint8_t samplingMode;

    // energy difference from relavent band
    double delE;
    double degeneracy;
} DopantInfo;

SemiConductor scInit(const BulkInfo bulkInfo, const Environment env);

void scDopeAcceptor(SemiConductor* sc, Mesh mesh, const DopantInfo dopantInfo);
void scDopeDonor(SemiConductor* sc, Mesh mesh, const DopantInfo dopantInfo);

double scTotalQ(SemiConductor sc, const double x, const double fermi_lvl, const double Ec);

void scVecTotalQ(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Vec* net_charge);
void scVecTotalQD(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Vec* net_charge_d);

double scSolveBoundary(SemiConductor sc, const double x, const double E, uint8_t flags);

void freeSemiConductor(SemiConductor* sc);

