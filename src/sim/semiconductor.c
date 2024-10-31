#include <include/sim/semiconductor.h>

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

double scTotalQ(SemiConductor sc, const double x, const double fermi_lvl, const double Ec, Environment env)
{
    double charge = -ELECTRON_CHARGE * bulkElectrons(sc.bulk, fermi_lvl, Ec, env);

    charge += ELECTRON_CHARGE * bulkHoles(sc.bulk, fermi_lvl, Ec - sc.bulk.bandgap, env);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        charge += dopantIonized(dopant, fermi_lvl, Ec, x, env);
    }
    return charge;
}

void scVecTotalQ(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* net_charge)
{
    Vec charge = vecInitZerosA(net_charge->len);

    bulkElectronsV(sc.bulk, fermi_lvl, Ec, env, net_charge);
    vecScale(-ELECTRON_CHARGE, *net_charge, net_charge);

    bulkHolesV(sc.bulk, fermi_lvl, Ec, env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge, charge, net_charge);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        dopantIonizedV(dopant, fermi_lvl, Ec, env, &charge);
        vecAdd(*net_charge, charge, net_charge);
    }

    freeVec(&charge);
}
void scVecTotalQD(SemiConductor sc, const Vec fermi_lvl, const Vec Ec, Environment env, Vec* net_charge_d)
{
    Vec charge = vecInitZerosA(net_charge_d->len);

    bulkElectronsDV(sc.bulk, fermi_lvl, Ec, env, net_charge_d);
    vecScale(-ELECTRON_CHARGE, *net_charge_d, net_charge_d);

    bulkHolesDV(sc.bulk, fermi_lvl, Ec, env, &charge);
    vecScale(ELECTRON_CHARGE, charge, &charge);
    vecAdd(*net_charge_d, charge, net_charge_d);

    for(size_t i = 0; i < sc.dopants.len; i++)
    {
        Dopant dopant = *(Dopant*)dynStackGet(sc.dopants, i);
        dopantIonizedDV(dopant, fermi_lvl, Ec, env, &charge);
        vecAdd(*net_charge_d, charge, net_charge_d);
    }

    freeVec(&charge);
}

void freeSemiConductor(SemiConductor* sc)
{
    for(size_t i = 0; i < sc->dopants.len; i++)
    {
        freeDopant(dynStackGet(sc->dopants, i));
    }

    freeDynStack(&sc->dopants);
}



