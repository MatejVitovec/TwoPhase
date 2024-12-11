#include <cmath>

#include "IdealGasThermo.hpp"

double IdealGasThermo::p(double T, double rho) const
{
    return rho*T*R;
}

double IdealGasThermo::rho(double p, double T) const
{
    return p/(T*R);
}

double IdealGasThermo::T(double rho, double p) const
{
    return p/(rho*R);
}

double IdealGasThermo::a(double rho, double p) const
{
    return std::sqrt((gamma*std::max(0.0, p))/rho);
}

double IdealGasThermo::e(double p, double T) const
{
    return p/((gamma - 1.0)*rho(p, T)) + energyShift;
}
        
double IdealGasThermo::pFromRho_e(double rho, double e) const
{
    return (gamma - 1.0)*rho*(e - energyShift);
}


Vars<3> IdealGasThermo::updateThermo(const Compressible& data, const ThermoVar& thermoData) const
{
    //return Vars<3>({temperature(data), pressure(data), soundSpeed(data)});

    const double density = data.density();
    const double pressure = pFromRho_e(density, data.internalEnergy());
    return Vars<3>({T(density, pressure), pressure, a(density, pressure)});
}

Vars<3> IdealGasThermo::updateThermo(const Primitive& data, const PrimitiveThermoVar& thermoData) const
{
    const double pressure = data.pressure();
    const double temperature = data.temperature();
    const double density = rho(pressure, temperature);
    return Vars<3>({density, e(pressure, temperature), a(density, pressure)});
}

void IdealGasThermo::updateThermo(ComponentThermoVar& thermoData) const
{
    const double density = thermoData.density();
    const double pressure = pFromRho_e(density, thermoData.internalEnergy());
    thermoData[ComponentThermoVar::T] = T(density, pressure);
    thermoData[ComponentThermoVar::P] = pressure;
    thermoData[ComponentThermoVar::A] = a(density, pressure);
}


Compressible IdealGasThermo::primitiveToConservative(const Vars<5>& primitive) const
{
    double velocity2 = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

    Compressible out = Compressible({primitive[0],
                                     primitive[0]*primitive[1],
                                     primitive[0]*primitive[2],
                                     primitive[0]*primitive[3],
                                     0.5*primitive[0]*velocity2 + (primitive[4])/(gamma - 1.0) + primitive[0]*energyShift});

    return out;
}

Compressible IdealGasThermo::stagnationState(double TTot, double pTot) const
{
    Compressible out = Compressible({pTot/(R*TTot),
                                     0.0,
                                     0.0,
                                     0.0,
                                     pTot/(gamma - 1.0) + pTot/(R*TTot)*energyShift});

    return out;
}


Compressible IdealGasThermo::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const
{
    double p = std::min(pFromRho_e(stateIn.density(), stateIn.internalEnergy()), pTot);
    double M2 = (2.0/(gamma - 1.0))*(std::pow((pTot/p), ((gamma - 1.0)/gamma)) - 1.0);
    double T = TTot/(1.0 + ((gamma - 1.0)/2)*M2);
    double rho = p/(R*T);
    double a = std::sqrt((gamma*p)/rho);
    double absU= std::sqrt(M2)*a;

    return Compressible({rho,
                         rho*absU*velocityDirection[0],
                         rho*absU*velocityDirection[1],
                         rho*absU*velocityDirection[2],
                         0.5*rho*absU*absU + p/(gamma - 1.0) + rho*energyShift});
}

std::array<double, 3> IdealGasThermo::initPressureTemperatureInlet(double pTot, double TTot) const
{
    return std::array<double, 3>({0.0, 0.0, 0.0});
}