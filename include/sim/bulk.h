#pragma once
#include <stdint.h>

#include <math/linalg.h>

#include "env.h"
#include <include/fem/functions.h>

// Bulk Semiconductor properties
typedef struct Bulk
{
    // in SI
    double epsilon;
    // in eV
    double bandgap;
     // in SI
    double effElectronMass;
     // in SI
    double effHoleMass;
} Bulk;

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkElectrons(Bulk bulk, const double fermi_lvl, const double Ec, const Environment env);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkHoles(Bulk bulk, const double fermi_lvl, const double Ev, const Environment env);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
// derivative wrt Ec
double bulkElectronsD(Bulk bulk, const double fermi_lvl, const double Ec, const Environment env);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
// derivative wrt Ec
double bulkHolesD(Bulk bulk, const double fermi_lvl, const double Ev, const Environment env);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronsV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* nconc);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHolesV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* hconc);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronsDV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* nconc);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHolesDV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* hconc);
