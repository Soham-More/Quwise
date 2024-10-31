#include "include/sim/material.h"

#include <math.h>
#include <include/utils.h>

#define DOPANT_HOLE 0x1
#define DOPANT_ELECTRON 0x2

// calculates Li_3/2(-e^(-x))
double calcPolyLogSeries(const double x, const double rel_error)
{
    double sum = 0;
    for(size_t k = 1;1;k++)
    {
        double t_k = exp(-(double)k * x) / pow(k, 1.5);

        if((k & 1) == 1) sum -= t_k;
        else           sum += t_k;

        // if t_k is within tolerence, return the value
        if(t_k < fabs(rel_error*sum) || t_k == 0.0) return sum;
    }
}
// calculate F_1/2(-x)
double calcFermiDiracHalfSeries(const double x, const double rel_error)
{
    const long double n_eta[] = 
    {
        0.7651470246254079453672687586034781795124679693458281781499491988,
        0.6048986434216303702472659142359554997597625451302473803785466480,
        0.3801048126096840167775421565518083625709371693918072509339302678,
        0.1186808707198402120435985572491929878560124016317782635845912186,
        -0.087841120721362842395232450051556648962734925662990349735537893,
        -0.096047604045123188661553125936681010842979342156641037163678475,
        0.1368213093058865486562233965254776791088062638036662212318929393,
        0.2391213204154196273351526507385141003549175466055537151570551835,
        -0.494471340576897219641496234666616819016701072128773235962114906,
        -1.180249705900827552847557670546601800135155723223697848022692353,
        3.1931333205421318063108292877510736947030664083035091442877553289,
        9.6556654171083717463012808156872975204555006128218758266716547726,
        -32.27147290195405328880232733747521541322771036312031676606381601,
        -118.1315244047967516255043546327349631179153334486332217922818883,
        470.03005997705269625724332628731639911620311427214666225221108048,
        2019.8049128538492730952445503994277187355777737010939113805322758,
    };

    size_t running_fact = 1;

    double sum = n_eta[0];

    for(size_t k = 1; k < sizeof(n_eta)/sizeof(long double); k++)
    {
        running_fact *= k;

        double t_k = (n_eta[k] / running_fact) * pow(x, k);

        if((k & 1) == 1) sum -= t_k;
        else           sum += t_k;

        // if t_k is within tolerence, return the value
        if(fabs(t_k) < rel_error*sum) return sum;
    }

    printf("Warning: concentration might be inaccurate: %le\n", x);

    return sum;
}

// calculates Li_3/2(-e^(-x))
double calcPolyLogSeriesDerivative(const double x, const double rel_error)
{
    double sum = 0;
    for(size_t k = 1;1;k++)
    {
        double t_k = -(double)k * exp(-(double)k * x) / pow(k, 1.5);

        if((k & 1) == 1) sum -= t_k;
        else           sum += t_k;

        // if t_k is within tolerence, return the value
        if(fabs(t_k) < fabs(rel_error*sum) || t_k == 0.0) return sum;
    }
}
// calculate F_1/2(-x)
double calcFermiDiracHalfSeriesDerivative(const double x, const double rel_error)
{
    const long double n_eta[] = 
    {
        0.7651470246254079453672687586034781795124679693458281781499491988,
        0.6048986434216303702472659142359554997597625451302473803785466480,
        0.3801048126096840167775421565518083625709371693918072509339302678,
        0.1186808707198402120435985572491929878560124016317782635845912186,
        -0.087841120721362842395232450051556648962734925662990349735537893,
        -0.096047604045123188661553125936681010842979342156641037163678475,
        0.1368213093058865486562233965254776791088062638036662212318929393,
        0.2391213204154196273351526507385141003549175466055537151570551835,
        -0.494471340576897219641496234666616819016701072128773235962114906,
        -1.180249705900827552847557670546601800135155723223697848022692353,
        3.1931333205421318063108292877510736947030664083035091442877553289,
        9.6556654171083717463012808156872975204555006128218758266716547726,
        -32.27147290195405328880232733747521541322771036312031676606381601,
        -118.1315244047967516255043546327349631179153334486332217922818883,
        470.03005997705269625724332628731639911620311427214666225221108048,
        2019.8049128538492730952445503994277187355777737010939113805322758,
        -9322.395392144822830535309362206116099640307039492393916925564803,
    };

    size_t running_fact = 1;

    double sum = n_eta[1];

    for(size_t k = 1; k < sizeof(n_eta)/sizeof(long double) - 1; k++)
    {
        running_fact *= k;

        double t_k = (n_eta[k + 1] / running_fact) * pow(x, k);

        if((k & 1) == 1) sum += t_k;
        else           sum -= t_k;

        // if t_k is within tolerence, return the value
        if(fabs(t_k) < rel_error*sum) return sum;
    }

    printf("Warning: concentration might be inaccurate: %le\n", x);

    return sum;
}

// doping = Na - Nd, doping > 0: p-type, doping < 0: n-type
// fermi-level and Ec should be in eV
// minimim tolerance supported is 1e-6 or 1e-4 % rel. error
double bulkElectronConcSingle(Bulk bulk, const double fermi_lvl, double Ec, const double T, const double tol)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    const double thermal_potential = BOLTZMANN_CONSTANT_EV * T;

    double density_of_states = 2 * pow(2*PI*bulk.effMassElectron*BOLTZMANN_CONSTANT*T, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (Ec - fermi_lvl) / thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ec-Ef/kbT = %le\n", r);

    if (r > 0.9)
    {
        return -density_of_states * calcPolyLogSeries(r, tol);
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeries(r, tol);
    }
}
double bulkHoleConcSingle(Bulk bulk, const double fermi_lvl, const double Ev, const double T, const double tol)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    const double thermal_potential = BOLTZMANN_CONSTANT_EV * T;

    double density_of_states = 2 * pow(2*PI*bulk.effMassHoles*BOLTZMANN_CONSTANT*T, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (fermi_lvl - Ev) / thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ef-Ev/kbT = %le\n", r);

    if (r > 1)
    {
        return -density_of_states * calcPolyLogSeries(r, tol);
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeries(r, tol);
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkElectronConcSingleDerivative(Bulk bulk, const double fermi_lvl, const double Ec, const double T, const double tol)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    const double thermal_potential = BOLTZMANN_CONSTANT_EV * T;

    double density_of_states = 2 * pow(2*PI*bulk.effMassElectron*BOLTZMANN_CONSTANT*T, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (Ec - fermi_lvl) / thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ec-Ef/kbT = %le\n", r);

    if (r > 0.9)
    {
        return -density_of_states * calcPolyLogSeriesDerivative(r, tol) / thermal_potential;
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeriesDerivative(r, tol) / thermal_potential;
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkHoleConcSingleDerivative(Bulk bulk, const double fermi_lvl, const double Ev, const double T, const double tol)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    const double thermal_potential = BOLTZMANN_CONSTANT_EV * T;

    double density_of_states = 2 * pow(2*PI*bulk.effMassHoles*BOLTZMANN_CONSTANT*T, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (fermi_lvl - Ev) / thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ef-Ev/kbT = %le\n", r);

    if (r > 1)
    {
        return density_of_states * calcPolyLogSeriesDerivative(r, tol) / thermal_potential;
    }
    else
    {
        return -density_of_states * calcFermiDiracHalfSeriesDerivative(r, tol) / thermal_potential;
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronConc(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* nconc)
{
    for(size_t i = 0; i < simParams.sample_count; i++)
    {
        VEC_INDEX(*nconc, i) = bulkElectronConcSingle(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i), simParams.temperature, simParams.tol_conc);
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHoleConc(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* hconc)
{
    for(size_t i = 0; i < simParams.sample_count; i++)
    {
        VEC_INDEX(*hconc, i) = bulkHoleConcSingle(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i) - bulk.bandgap, simParams.temperature, simParams.tol_conc);
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronConcDerivative(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* nconc)
{
    for(size_t i = 0; i < simParams.sample_count; i++)
    {
        VEC_INDEX(*nconc, i) = bulkElectronConcSingleDerivative(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i), simParams.temperature, simParams.tol_conc);
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHoleConcDerivative(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const SimParameters simParams, Vec* hconc)
{
    for(size_t i = 0; i < simParams.sample_count; i++)
    {
        VEC_INDEX(*hconc, i) = bulkHoleConcSingleDerivative(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i) - bulk.bandgap, simParams.temperature, simParams.tol_conc);
    }
}

Dopant dopantConstructAcceptor1DA(Fx1D conc, const double Ed_Ev, const double degeneracy, const Bulk bulk, SimParameters simParams)
{
    Dopant dopant;
    dopant.delE = Ed_Ev;
    dopant.type = DOPANT_HOLE;
    dopant.degeneracy = degeneracy;
    dopant.bulk = bulk;
    dopant.conc = conc;
    
    dopant.interp_conc = vecInitZerosA(simParams.sample_count);
    fxInterpolateSample1D(conc, simParams.x_samples, &dopant.interp_conc);

    return dopant;
}
Dopant dopantConstructDonor1DA(Fx1D conc, const double Ec_Ed, const double degeneracy, const Bulk bulk, SimParameters simParams)
{
    Dopant dopant;
    dopant.delE = Ec_Ed;
    dopant.type = DOPANT_ELECTRON;
    dopant.degeneracy = degeneracy;
    dopant.bulk = bulk;
    dopant.conc = conc;
    
    dopant.interp_conc = vecInitZerosA(simParams.sample_count);
    fxInterpolateSample1D(conc, simParams.x_samples, &dopant.interp_conc);

    return dopant;
}

// calculate the amount of dopant ionized
double dopantIonizedConcSingle(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, SimParameters simParams)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        double Ed = Ec - dopant.delE;
        double r = (fermi_lvl - Ed) / simParams.thermal_potential;

        return fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        double Ed = Ec - dopant.bulk.bandgap + dopant.delE;
        double r = (Ed - fermi_lvl) / simParams.thermal_potential;

        return fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
    return NAN;
}
// calculate the amount of dopant ionized
double dopantIonizedChargeSingle(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, SimParameters simParams)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        double Ed = Ec - dopant.delE;
        double r = (fermi_lvl - Ed) / simParams.thermal_potential;

        return ELECTRON_CHARGE * fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        double Ed = Ec - dopant.bulk.bandgap + dopant.delE;
        double r = (Ed - fermi_lvl) / simParams.thermal_potential;

        return -ELECTRON_CHARGE * fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
    return NAN;
}

// calculate the amount of dopant ionized, 
void dopantIonizedConc(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* ionized_conc)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / simParams.thermal_potential;

            VEC_INDEX(*ionized_conc, i) = VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / simParams.thermal_potential;

            VEC_INDEX(*ionized_conc, i) = VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}
// calculate the the charge conc due to dopants
void dopantIonizedCharge(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* charge)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / simParams.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / simParams.thermal_potential;

            VEC_INDEX(*charge, i) = -ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}

void dopantIonizedChargeDerivative(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, SimParameters simParams, Vec* charge)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / simParams.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) * dopant.degeneracy * exp(r) / (simParams.thermal_potential * (1 + dopant.degeneracy * exp(r))*(1 + dopant.degeneracy * exp(r)));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < simParams.sample_count; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / simParams.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) * dopant.degeneracy * exp(r) / (simParams.thermal_potential * (1 + dopant.degeneracy * exp(r))*(1 + dopant.degeneracy * exp(r)));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}

void freeDopant(Dopant* dopant)
{
    freeVec(&dopant->interp_conc);
}

