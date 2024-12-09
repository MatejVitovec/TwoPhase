#include <cmath>
#include <iostream>
#include <fstream>

#include "Iapws95Thermo.hpp"


Vars<3> Iapws95Thermo::updateThermo(const Compressible& data, const ThermoVar& thermoData) const
{
    double rho = data.density();
    //double TGuess = ((1.29 - 1.0)*(data.totalEnergy() - 0.5*data.absVelocity2()))/specGasConst; //ideal gas guess
    double T = tFromRhoE(rho, data.internalEnergy(), thermoData.temperature());

    return Vars<3>({T, p(rho, T), a(rho, T)});
}

Compressible Iapws95Thermo::primitiveToConservative(const Vars<5>& primitive) const
{
    double rho = primitive[0];
    double p = primitive[4];
    double T = tFromRhoP(rho, p, p/(specGasConst*rho)); //odhad T pomoci idealniho plynu

    return Compressible({rho,
                         rho*primitive[1],
                         rho*primitive[2],
                         rho*primitive[3],
                         rho*(e(rho, T) + 0.5*(primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3]))});
}

Compressible Iapws95Thermo::stagnationState(double TTot, double pTot) const
{
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return Compressible({rhoTot,
                         0.0,
                         0.0,
                         0.0,
                         rhoTot*(e(rhoTot, TTot))});
}

double Iapws95Thermo::calcM2is(double p0, double T0, double p2) const
{
    double rho0 = rhoFromTP(T0, p0, p0/(specGasConst*T0));
    double s0 = s(rho0, T0);

    std::pair<double, double> state2 = RhoTFromSP(s0, p2, rho0, T0*1.2);
    double h0 = h(rho0, T0);
    double h2 = h(state2.first, state2.second);

    return std::sqrt((2*(h0 - h2))/std::pow(a(state2.first, state2.second), 2));
}

Compressible Iapws95Thermo::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const
{
    thermoIn = updateThermo(stateIn, thermoIn);

    double pIn = std::min(thermoIn.pressure(), pTot);

    double guessRho = stateIn.density();
    double guessT = thermoIn.temperature();

    std::pair<double, double> result = RhoTFromSP(sTot, pIn, guessRho, guessT);

    double rho = result.first;
    double T = result.second;

    double absU2 = std::max(2.0*hTot - 2.0*h(rho, T), 0.0);
    double absU = std::sqrt(absU2);

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         rho*(0.5*absU2 + e(rho, T))});
}

std::array<double, 3> Iapws95Thermo::initPressureTemperatureInlet(double pTot, double TTot) const
{
    double rhoTot = rhoFromTP(TTot, pTot, pTot/(specGasConst*TTot));

    return std::array<double, 3>({rhoTot, s(rhoTot, TTot), h(rhoTot, TTot)});
}