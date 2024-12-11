#include <cmath>
#include "ThermoVar.hpp"

double ThermoVar::temperature() const
{
    return data[T];
}

double ThermoVar::pressure() const
{
    return data[P];
}

double ThermoVar::soundSpeed() const
{
    return data[A];
}