#include <cmath>
#include "TwoFluid.hpp"

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
    return thermoData[INT_E_G];
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

double TwoFluid::totalEnergyG() const
{
    return thermoData[INT_E_G] + 0.5*norm2sqr(this->velocityG());
}

double TwoFluid::totalEnergyL() const
{
    return thermoData[INT_E_L] + 0.5*norm2sqr(this->velocityL());
}

double TwoFluid::interfacialPressure() const
{
    constexpr double sigma = 2.0;
    constexpr double epsilonP = 0.01;

    const double delta_pInt = sigma*((data[ALPHA]*thermoData[RHO_G]*(1.0 - data[ALPHA])*thermoData[RHO_L])
        /(data[ALPHA]*thermoData[RHO_L] + (1.0 - data[ALPHA])*thermoData[RHO_G]))
        *norm2sqr(this->velocityL() - this->velocityG()); //NUTNO OVERTIT

    return data[P] - std::min(delta_pInt, epsilonP*data[P]); //pInt = p - deltaP* (ve clanku "stratified model" od Lia je to pInt = deltaP*)
}

void TwoFluid::blend()
{
    if (alphaL() <= epsilonMin)
    {
        data[ALPHA] = 1.0 - epsilonMin;
        data[U_L] = data[U_G];
        data[V_L] = data[V_G];
        data[W_L] = data[W_G];
        data[T_L] = data[T_G];
    }
    else if (alphaL() < epsilonMax)
    {
        const double xi = (alphaL() - epsilonMin)/(epsilonMax - epsilonMin);
        const double gFunc = -std::pow(xi, 2.0)*(2.0*xi - 3.0);

        data[U_L] = gFunc*data[U_L] + (1.0 - gFunc)*data[U_G];
        data[V_L] = gFunc*data[V_L] + (1.0 - gFunc)*data[V_G];
        data[W_L] = gFunc*data[W_L] + (1.0 - gFunc)*data[W_G];
        data[T_L] = gFunc*data[T_L] + (1.0 - gFunc)*data[T_G];
    }
    
    if (alphaG() <= epsilonMin)
    {
        data[ALPHA] = epsilonMin;
        data[U_G] = data[U_L];
        data[V_G] = data[V_L];
        data[W_G] = data[W_L];
        data[T_G] = data[T_L];
    }
    else if (alphaG() < epsilonMax)
    {
        const double xi = (alphaG() - epsilonMin)/(epsilonMax - epsilonMin);
        const double gFunc = -std::pow(xi, 2.0)*(2.0*xi - 3.0);

        data[U_G] = gFunc*data[U_G] + (1.0 - gFunc)*data[U_L];
        data[V_G] = gFunc*data[V_G] + (1.0 - gFunc)*data[V_L];
        data[W_G] = gFunc*data[W_G] + (1.0 - gFunc)*data[W_L];
        data[T_G] = gFunc*data[T_G] + (1.0 - gFunc)*data[T_L];
    }
}

std::array<double, 6>& TwoFluid::getThermoData()
{
    return thermoData;
}

const std::array<double, 6>& TwoFluid::getThermoData() const
{
    return thermoData;
}