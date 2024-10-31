#pragma once
#include <stdint.h>

#include <include/sim.h>
#include <math/linalg.h>
#include <include/functions.h>

// Bulk Semiconductor properties
typedef struct Bulk
{
    double epsilon; // in SI
    double bandgap; // in eV
    double effMassElectron; // in SI
    double effMassHoles; // in SI
} Bulk;

// doping properties, multiple dopants supported
typedef struct Dopant
{
    double delE; // in eV, wrt Ev/Ec depending on dopant type
    double degeneracy; // the degeneracy of the dopant
    uint8_t type; // type of dopant
    Vec interp_conc; // interpolated conc data
    Bulk bulk; // the SC bulk, this dopant is part of
    Fx1D conc;
} Dopant;

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkElectronConcSingle(Bulk bulk, const double fermi_lvl, const double Ec, const double T, const double tol);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkHoleConcSingle(Bulk bulk, const double fermi_lvl, const double Ev, const double T, const double tol);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkElectronConcSingleDerivative(Bulk bulk, const double fermi_lvl, const double Ec, const double T, const double tol);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkHoleConcSingleDerivative(Bulk bulk, const double fermi_lvl, const double Ev, const double T, const double tol);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronConc(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* nconc);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHoleConc(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* hconc);

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronConcDerivative(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* nconc);
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHoleConcDerivative(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* hconc);

// allocates memory
Dopant dopantConstructAcceptor1DA(Fx1D conc, const double Ed_Ev, const double degeneracy, Bulk bulk, SimParameters simParams);
// allocates memory
Dopant dopantConstructDonor1DA(Fx1D conc, const double Ec_Ed, const double degeneracy, Bulk bulk, SimParameters simParams);

// calculate the amount of dopant ionized
double dopantIonizedConcSingle(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, SimParameters simParams);
// calculate the the charge conc due to dopants
double dopantIonizedChargeSingle(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, SimParameters simParams);

// calculate the amount of dopant ionized
void dopantIonizedConc(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* ionized_conc);
// calculate the the charge conc due to dopants
void dopantIonizedCharge(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* charge);

// calculate the the charge conc due to dopants
void dopantIonizedChargeDerivative(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* charge);

void freeDopant(Dopant* dopant);

