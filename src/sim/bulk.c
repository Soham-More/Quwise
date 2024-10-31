#include "include/sim/bulk.h"

#include <math.h>
#include <include/utils.h>

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
double bulkElectrons(Bulk bulk, const double fermi_lvl, double Ec, const Environment env)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    double density_of_states = 2 * pow(2*PI*bulk.effMassElectron*BOLTZMANN_CONSTANT*env.temperature, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (Ec - fermi_lvl) / env.thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ec-Ef/kbT = %le\n", r);

    if (r > 0.9)
    {
        return -density_of_states * calcPolyLogSeries(r, env.tol_conc);
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeries(r, env.tol_conc);
    }
}
double bulkHoles(Bulk bulk, const double fermi_lvl, const double Ev, const Environment env)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    double density_of_states = 2 * pow(2*PI*bulk.effMassHoles*BOLTZMANN_CONSTANT*env.temperature, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (fermi_lvl - Ev) / env.thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ef-Ev/kbT = %le\n", r);

    if (r > 1)
    {
        return -density_of_states * calcPolyLogSeries(r, env.tol_conc);
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeries(r, env.tol_conc);
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkElectronsD(Bulk bulk, const double fermi_lvl, const double Ec, const Environment env)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    double density_of_states = 2 * pow(2*PI*bulk.effMassElectron*BOLTZMANN_CONSTANT*env.temperature, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (Ec - fermi_lvl) / env.thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ec-Ef/kbT = %le\n", r);

    if (r > 0.9)
    {
        return -density_of_states * calcPolyLogSeriesDerivative(r, env.tol_conc) / env.thermal_potential;
    }
    else
    {
        return density_of_states * calcFermiDiracHalfSeriesDerivative(r, env.tol_conc) / env.thermal_potential;
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
double bulkHolesD(Bulk bulk, const double fermi_lvl, const double Ev, const Environment env)
{
    // 2 strategies: if Ec - Ef > thermal_potential: use polylogarithm power series
    // if Ec - Ef < thermal_potential: use fermi-dirac function power series

    double density_of_states = 2 * pow(2*PI*bulk.effMassHoles*BOLTZMANN_CONSTANT*env.temperature, 1.5) / pow(PLACK_CONSTANT, 3);

    double r = (fermi_lvl - Ev) / env.thermal_potential;

    if (r < 0) printf("Warning: Calculating conc in untested range Ef-Ev/kbT = %le\n", r);

    if (r > 1)
    {
        return density_of_states * calcPolyLogSeriesDerivative(r, env.tol_conc) / env.thermal_potential;
    }
    else
    {
        return -density_of_states * calcFermiDiracHalfSeriesDerivative(r, env.tol_conc) / env.thermal_potential;
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronsV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* nconc)
{
    for(size_t i = 0; i < fermi_lvl.len; i++)
    {
        VEC_INDEX(*nconc, i) = bulkElectrons(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i), env);
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHolesV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* hconc)
{
    for(size_t i = 0; i < fermi_lvl.len; i++)
    {
        VEC_INDEX(*hconc, i) = bulkHoles(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i) - bulk.bandgap, env);
    }
}

// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkElectronsDV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* nconc)
{
    for(size_t i = 0; i < fermi_lvl.len; i++)
    {
        VEC_INDEX(*nconc, i) = bulkElectronsD(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i), env);
    }
}
// fermi-level, and Ec should be in eV, tolerance is (max) relative error permitted
// assumes Ec - Ef > 0
void bulkHolesDV(Bulk bulk, const Vec fermi_lvl, const Vec Ec, const Environment env, Vec* hconc)
{
    for(size_t i = 0; i < fermi_lvl.len; i++)
    {
        VEC_INDEX(*hconc, i) = bulkHolesD(bulk, VEC_INDEX(fermi_lvl, i), VEC_INDEX(Ec, i) - bulk.bandgap, env);
    }
}
