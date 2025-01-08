#include <cmath>
#include "Fluid.hpp"

double Fluid::density() const
{
    return thermoData[RHO];
}

double Fluid::internalEnergy() const
{
    return thermoData[INT_E];
}

double Fluid::totalEnergy() const
{
    return thermoData[INT_E] + 0.5*absVelocity2();
}

double Fluid::soundSpeed() const
{
    return thermoData[A];
}

std::array<double, 3>& Fluid::getThermoData()
{
    return thermoData;
}

const std::array<double, 3>& Fluid::getThermoData() const
{
    return thermoData;
}