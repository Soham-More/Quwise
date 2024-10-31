#include <include/sim/semiconductor.h>

#include <math.h>

#include <include/utils.h>
#include <include/sim/dopant.h>

SemiConductor scInit(const BulkInfo bulkInfo, const Environment env)
{
    SemiConductor sc;

    sc.bulk.bandgap = bulkInfo.bandgap;
    sc.bulk.effElectronMass = bulkInfo.relativeEffElectronMass * ELECTRON_MASS;
    sc.bulk.effHoleMass = bulkInfo.relativeEffHoleMass * ELECTRON_MASS;
    sc.bulk.epsilon = bulkInfo.relativePermittivity * VACCUME_PERMITTIVITY;
    sc.env = env;
    sc.dopants = dynStackInit(sizeof(Dopant));

    return sc;
}

void scDopeAcceptor(SemiConductor* sc, Mesh mesh, const DopantInfo dopantInfo)
{
    Fx1D conc = fxConstruct1D(dopantInfo.sampled_conc, dopantInfo.sampled_x, dopantInfo.samplingMode);

    Dopant acceptor = dopantInitAcceptor1DA(conc, mesh.x, dopantInfo.delE, dopantInfo.degeneracy, sc->bulk, sc->env);

    dynStackPush(&sc->dopants, &acceptor);
}
void scDopeDonor(SemiConductor* sc, Mesh mesh, const DopantInfo dopantInfo)
{
    Fx1D conc = fxConstruct1D(dopantInfo.sampled_conc, dopantInfo.sampled_x, dopantInfo.samplingMode);

    Dopant acceptor = dopantInitDonor1DA(conc, mesh.x, dopantInfo.delE, dopantInfo.degeneracy, sc->bulk, sc->env);

    dynStackPush(&sc->dopants, &acceptor);
}

double scTotalQ(SemiConductor sc, const double x, const double fermi_lvl, const double Ec)
{
    double charge = -ELECTRON_CHARGE * bulkElectrons(sc.bulk, fermi_lvl, Ec, sc.env);

    charge += ELECTRON_CHARGE * bulkHoles(sc.bulk, fermi_lvl, Ec - sc.bulk.bandgap, sc.env);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        charge += dopantIonized(dopant, fermi_lvl, Ec, x, sc.env);
    }
    return charge;
}

void scVecTotalQ(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Vec* net_charge)
{
    // TODO: Optimize this
    Vec charge = vecInitZerosA(net_charge->len);

    bulkElectronsV(sc.bulk, fermi_lvl, Ec, sc.env, net_charge);
    vecScale(-ELECTRON_CHARGE, *net_charge, net_charge);

    bulkHolesV(sc.bulk, fermi_lvl, Ec, sc.env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge, charge, net_charge);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        dopantIonizedV(dopant, fermi_lvl, Ec, sc.env, &charge);
        vecAdd(*net_charge, charge, net_charge);
    }

    freeVec(&charge);
}
void scVecTotalQD(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Vec* net_charge_d)
{
    // TODO: Optimize this
    Vec charge = vecInitZerosA(net_charge_d->len);

    bulkElectronsDV(sc.bulk, fermi_lvl, Ec, sc.env, net_charge_d);
    vecScale(-ELECTRON_CHARGE, *net_charge_d, net_charge_d);

    bulkHolesDV(sc.bulk, fermi_lvl, Ec, sc.env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge_d, charge, net_charge_d);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        dopantIonizedDV(dopant, fermi_lvl, Ec, sc.env, &charge);
        vecAdd(*net_charge_d, charge, net_charge_d);
    }

    freeVec(&charge);
}

double scSolveBEC(SemiConductor sc, const double x, const double Ef)
{
    double low_Ec = Ef;
    double high_Ec = Ef + sc.bulk.bandgap;

    double low_value = scTotalQ(sc, x, Ef, low_Ec);
    double high_value = scTotalQ(sc, x, Ef, high_Ec);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_Ec - low_Ec) > sc.env.tol_potential*fabs(high_Ec + low_Ec))
    {
        double new_Ec = (low_Ec + high_Ec) / 2;
        double error = scTotalQ(sc, x, Ef, new_Ec);

        if(error*low_value < 0)
        {
            high_Ec = new_Ec;
            high_value = error;
        }
        else
        {
            low_Ec = new_Ec;
            low_value = error;
        }
    }

    return (high_Ec + low_Ec) / 2;
}
double scSolveBEF(SemiConductor sc, const double x, const double Ec)
{
    double low_u = Ec;
    double high_u = Ec - sc.bulk.bandgap;

    double low_value = scTotalQ(sc, x, low_u, Ec);
    double high_value = scTotalQ(sc, x, high_u, Ec);

    if(low_value*high_value > 0)
    {
        printf("Error: Both bounds are less than zero!\n");
        return NAN;
    }

    while(fabs(high_u - low_u) > sc.env.tol_potential*fabs(high_u + low_u))
    {
        double new_u = (low_u + high_u) / 2;

        double error = scTotalQ(sc, x, new_u, Ec);

        if(error*low_value < 0)
        {
            high_u = new_u;
            high_value = error;
        }
        else
        {
            low_u = new_u;
            low_value = error;
        }
    }

    return (high_u + low_u) / 2;
}

double scSolveBoundary(SemiConductor sc, const double x, const double E, uint8_t flags)
{
    if(flags & SC_SOLVE_EC) return scSolveBEC(sc, x, E);
    if(flags & SC_SOLVE_FERMI) return scSolveBEF(sc, x, E);
    
    printf("[Semiconductor] Error: scSolveBoundary flags is invalid %zu", flags);
    return NAN;
}

void freeSemiConductor(SemiConductor* sc)
{
    for(size_t i = 0; i < sc->dopants.len; i++)
    {
        freeDopant(dynStackGet(sc->dopants, i));
    }

    freeDynStack(&sc->dopants);
}



