#pragma once
#include <stdint.h>

#include <math/linalg.h>

#include "env.h"
#include "bulk.h"

#include <include/fem/functions.h>

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

// allocates memory
Dopant dopantInitAcceptor1DA(Fx1D conc, Vec mesh, const double Ed_Ev, const double degeneracy, Bulk bulk);
// allocates memory
Dopant dopantInitDonor1DA(Fx1D conc, Vec mesh, const double Ec_Ed, const double degeneracy, Bulk bulk);

// calculate the amount of dopant ionized
double dopantIonizedConc(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, Environment env);
// calculate the the charge conc due to dopants
double dopantIonized(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, Environment env);

// calculate the amount of dopant ionized
void dopantIonizedConcV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* ionized_conc);
// calculate the the charge conc due to dopants
void dopantIonizedV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* charge);

// calculate the the charge conc due to dopants
void dopantIonizedDV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* charge);

void freeDopant(Dopant* dopant);


