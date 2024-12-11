#include <cmath>
#include "PrimitiveThermoVar.hpp"

double PrimitiveThermoVar::density() const
{
    return data[RHO];
}

double PrimitiveThermoVar::internalEnergy() const
{
    return data[INT_E];
}

double PrimitiveThermoVar::soundSpeed() const
{
    return data[A];
}