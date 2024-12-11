#include <cmath>
#include "TwoFluid.hpp"

/*double TwoFluid::alpha() const
{
    return data[ALPHA];
}

double TwoFluid::pressure() const
{
    return thermoData[P];
}

double TwoFluid::temperatureG() const
{
    return data[T_G];
}

double TwoFluid::temperatureL() const
{
    return data[T_L];
}

Vars<3> TwoFluid::velocityG() const
{
    return Vars<3>({data[U_G], data[V_G], data[W_G]});
}

double TwoFluid::absVelocityG() const
{
    return std::sqrt(this->absVelocity2G());
}

double TwoFluid::absVelocity2G() const
{
    return data[U_G]*data[U_G] + data[V_G]*data[V_G] + data[W_G]*data[W_G];
}

double TwoFluid::normalVelocityG(const Vars<3>& normalVector) const
{
    return dot(this->velocityG(), normalVector);
}

double TwoFluid::velocityUG() const
{
    return data[U_G];
}

double TwoFluid::velocityVG() const
{
    return data[V_G];
}

double TwoFluid::velocityWG() const
{
    return data[W_G];
}



Vars<3> TwoFluid::velocityL() const
{
    return Vars<3>({data[U_L], data[V_L], data[W_L]});
}

double TwoFluid::absVelocityL() const
{
    return std::sqrt(this->absVelocity2L());
}

double TwoFluid::absVelocity2L() const
{
    return data[U_L]*data[U_L] + data[V_L]*data[V_L] + data[W_L]*data[W_L];
}

double TwoFluid::normalVelocityL(const Vars<3>& normalVector) const
{
    return dot(this->velocityL(), normalVector);
}

double TwoFluid::velocityUL() const
{
    return data[U_L];
}

double TwoFluid::velocityVL() const
{
    return data[V_L];
}

double TwoFluid::velocityWL() const
{
    return data[W_L];
}*/


double TwoFluid::densityG() const
{
    return thermoData[RHO_G];
}

double TwoFluid::densityL() const
{
    return thermoData[RHO_L];
}

double TwoFluid::internalEnergyG() const
{
    return thermoData[INT_E_L];
}

double TwoFluid::internalEnergyL() const
{
    return thermoData[INT_E_L];
}

double TwoFluid::soundSpeedG() const
{
    return thermoData[A_G];
}

double TwoFluid::soundSpeedL() const
{
    return thermoData[A_L];
}

double TwoFluid::interfacialPressure() const
{
    constexpr double sigma = 2.0;
    constexpr double epsilonP = 0.01;

    double pInt = sigma*((data[ALPHA]*thermoData[RHO_G]*(1.0 - data[ALPHA])*thermoData[RHO_L])
        /(data[ALPHA]*thermoData[RHO_L] + (1.0 - data[ALPHA])*thermoData[RHO_G]))*norm2(this->velocityG() - this->velocityL()); //NUTNO OVERTIT

    return std::min(pInt, epsilonP*data[P]);
}

void TwoFluid::blend()
{
    //TODO implement
}

std::array<double, 6>& TwoFluid::getThermoData()
{
    return thermoData;
}

const std::array<double, 6>& TwoFluid::getThermoData() const
{
    return thermoData;
}