#include <include/sim/dopant.h>

#include <math.h>
#include <include/utils.h>

#define DOPANT_HOLE 0x1
#define DOPANT_ELECTRON 0x2

Dopant dopantInitAcceptor1DA(Fx1D conc, Vec mesh, const double Ed_Ev, const double degeneracy, const Bulk bulk, Environment env)
{
    Dopant dopant;
    dopant.delE = Ed_Ev;
    dopant.type = DOPANT_HOLE;
    dopant.degeneracy = degeneracy;
    dopant.bulk = bulk;
    dopant.conc = conc;
    
    dopant.interp_conc = vecInitZerosA(mesh.len);
    fxInterpolateSample1D(conc, mesh, &dopant.interp_conc);

    return dopant;
}
Dopant dopantInitDonor1DA(Fx1D conc, Vec mesh, const double Ec_Ed, const double degeneracy, const Bulk bulk, Environment env)
{
    Dopant dopant;
    dopant.delE = Ec_Ed;
    dopant.type = DOPANT_ELECTRON;
    dopant.degeneracy = degeneracy;
    dopant.bulk = bulk;
    dopant.conc = conc;
    
    dopant.interp_conc = vecInitZerosA(mesh.len);
    fxInterpolateSample1D(conc, mesh, &dopant.interp_conc);

    return dopant;
}

// calculate the amount of dopant ionized
double dopantIonizedConc(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, Environment env)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        double Ed = Ec - dopant.delE;
        double r = (fermi_lvl - Ed) / env.thermal_potential;

        return fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        double Ed = Ec - dopant.bulk.bandgap + dopant.delE;
        double r = (Ed - fermi_lvl) / env.thermal_potential;

        return fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
    return NAN;
}
// calculate the amount of dopant ionized
double dopantIonized(const Dopant dopant, const double fermi_lvl, const double Ec, const double x, Environment env)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        double Ed = Ec - dopant.delE;
        double r = (fermi_lvl - Ed) / env.thermal_potential;

        return ELECTRON_CHARGE * fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        double Ed = Ec - dopant.bulk.bandgap + dopant.delE;
        double r = (Ed - fermi_lvl) / env.thermal_potential;

        return -ELECTRON_CHARGE * fxInterpolateSingle1D(dopant.conc, x) / (1 + dopant.degeneracy * exp(r));
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
    return NAN;
}

// calculate the amount of dopant ionized, 
void dopantIonizedConcV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* ionized_conc)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / env.thermal_potential;

            VEC_INDEX(*ionized_conc, i) = VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / env.thermal_potential;

            VEC_INDEX(*ionized_conc, i) = VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}
// calculate the the charge conc due to dopants
void dopantIonizedV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* charge)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / env.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / env.thermal_potential;

            VEC_INDEX(*charge, i) = -ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) / (1 + dopant.degeneracy * exp(r));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}

void dopantIonizedDV(const Dopant dopant, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* charge)
{
    if(dopant.type & DOPANT_ELECTRON)
    {
        // TODO: vectorize this loop
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.delE;
            double r = (VEC_INDEX(fermi_lvl, i) - Ed) / env.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) * dopant.degeneracy * exp(r) / (env.thermal_potential * (1 + dopant.degeneracy * exp(r))*(1 + dopant.degeneracy * exp(r)));
        }
        return;
    }
    else if(dopant.type & DOPANT_HOLE)
    {
        for(size_t i = 0; i < fermi_lvl.len; i++)
        {
            double Ed = VEC_INDEX(Ec, i) - dopant.bulk.bandgap + dopant.delE;
            double r = (Ed - VEC_INDEX(fermi_lvl, i)) / env.thermal_potential;

            VEC_INDEX(*charge, i) = ELECTRON_CHARGE * VEC_INDEX(dopant.interp_conc, i) * dopant.degeneracy * exp(r) / (env.thermal_potential * (1 + dopant.degeneracy * exp(r))*(1 + dopant.degeneracy * exp(r)));
        }
        return;
    }

    printf("Warning: invalid dopant type 0x%x\n", dopant.type);
}

void freeDopant(Dopant* dopant)
{
    freeVec(&dopant->interp_conc);
}

