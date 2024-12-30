#include <cmath>
#include "TwoFluidCompressible.hpp"

TwoFluidCompressible::TwoFluidCompressible(const TwoFluid& u) : Vars<10>()
{
    const double pInt = u.interfacialPressure();

    data[ALPHA_RHO_G]        = u.alphaG()*u.densityG();
    data[ALPHA_RHO_U_G]      = u.alphaG()*u.densityG()*u.velocityUG();
    data[ALPHA_RHO_V_G]      = u.alphaG()*u.densityG()*u.velocityVG();
    data[ALPHA_RHO_W_G]      = u.alphaG()*u.densityG()*u.velocityWG();
    data[ALPHA_RHO_E_PINT_G] = u.alphaG()*(u.densityG()*u.totalEnergyG() + pInt);

    data[ALPHA_RHO_L]        = u.alphaL()*u.densityL();
    data[ALPHA_RHO_U_L]      = u.alphaL()*u.densityL()*u.velocityUL();
    data[ALPHA_RHO_V_L]      = u.alphaL()*u.densityL()*u.velocityVL();
    data[ALPHA_RHO_W_L]      = u.alphaL()*u.densityL()*u.velocityWL();
    data[ALPHA_RHO_E_PINT_L] = u.alphaL()*(u.densityL()*u.totalEnergyL() + pInt);
}

void TwoFluidCompressible::updateInterfacialPressure(const double pInt, const TwoFluid& u)
{
    data[ALPHA_RHO_E_PINT_G] = u.alphaG()*(u.densityG()*u.totalEnergyG() + pInt);
    data[ALPHA_RHO_E_PINT_L] = u.alphaL()*(u.densityL()*u.totalEnergyL() + pInt);
}

double TwoFluidCompressible::alphaDensityG() const
{
    return data[ALPHA_RHO_G];
}

double TwoFluidCompressible::densityG(double alpha)
{
    return data[ALPHA_RHO_G]/alpha;
}

Vars<3> TwoFluidCompressible::velocityG() const
{
    return Vars<3>({data[ALPHA_RHO_U_G] / data[ALPHA_RHO_G], data[ALPHA_RHO_V_G] / data[ALPHA_RHO_G], data[ALPHA_RHO_W_G] / data[ALPHA_RHO_G]});
}

double TwoFluidCompressible::absVelocityG() const
{
    return sqrt(this->absVelocity2G());
}

double TwoFluidCompressible::absVelocity2G() const
{
    return (data[ALPHA_RHO_U_G]*data[ALPHA_RHO_U_G] + data[ALPHA_RHO_V_G]*data[ALPHA_RHO_V_G] + data[ALPHA_RHO_W_G]*data[ALPHA_RHO_W_G]) / (data[ALPHA_RHO_G]*data[ALPHA_RHO_G]);
}

double TwoFluidCompressible::normalVelocityG(const Vars<3>& normalVector) const
{
    return dot(this->velocityG(), normalVector);
}

double TwoFluidCompressible::velocityUG() const
{
    return data[ALPHA_RHO_U_G] / data[ALPHA_RHO_G];
}

double TwoFluidCompressible::velocityVG() const
{
    return data[ALPHA_RHO_V_G] / data[ALPHA_RHO_G];
}

double TwoFluidCompressible::velocityWG() const
{
    return data[ALPHA_RHO_W_G] / data[ALPHA_RHO_G];
}

///////////////

double TwoFluidCompressible::alphaDensityL() const
{
    return data[ALPHA_RHO_L];
}

double TwoFluidCompressible::densityL(double alpha)
{
    return data[ALPHA_RHO_L]/alpha;
}

Vars<3> TwoFluidCompressible::velocityL() const
{
    return Vars<3>({data[ALPHA_RHO_U_L] / data[ALPHA_RHO_L], data[ALPHA_RHO_V_L] / data[ALPHA_RHO_L], data[ALPHA_RHO_W_L] / data[ALPHA_RHO_L]});
}

double TwoFluidCompressible::absVelocityL() const
{
    return sqrt(this->absVelocity2L());
}

double TwoFluidCompressible::absVelocity2L() const
{
    return (data[ALPHA_RHO_U_L]*data[ALPHA_RHO_U_L] + data[ALPHA_RHO_V_L]*data[ALPHA_RHO_V_L] + data[ALPHA_RHO_W_L]*data[ALPHA_RHO_W_L]) / (data[ALPHA_RHO_L]*data[ALPHA_RHO_L]);
}

double TwoFluidCompressible::normalVelocityL(const Vars<3>& normalVector) const
{
    return dot(this->velocityL(), normalVector);
}

double TwoFluidCompressible::velocityUL() const
{
    return data[ALPHA_RHO_U_L] / data[ALPHA_RHO_L];
}

double TwoFluidCompressible::velocityVL() const
{
    return data[ALPHA_RHO_V_L] / data[ALPHA_RHO_L];
}

double TwoFluidCompressible::velocityWL() const
{
    return data[ALPHA_RHO_W_L] / data[ALPHA_RHO_L];
}

/*Vars<5> TwoFluidCompressible::flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const
{
    Vars<3> velocity = this->velocity();

    double normalVelocity = dot(velocity, normalVector);

    return TwoFluidCompressible({data[RHO] * normalVelocity,
                         data[RHO] * velocity[0] * normalVelocity + thermoData[1] * normalVector[0],
                         data[RHO] * velocity[1] * normalVelocity + thermoData[1] * normalVector[1],
                         data[RHO] * velocity[2] * normalVelocity + thermoData[1] * normalVector[2],
                         (data[RHO_E] + thermoData[1]) * normalVelocity});
}*/