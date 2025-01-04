#include <cmath>
#include "TwoFluidCompressible.hpp"

TwoFluidCompressible::TwoFluidCompressible(const TwoFluid& u) : Vars<10>()
{
    const double pInt = u.interfacialPressure();
    const double alphaRhoG = u.alphaG()*u.densityG();
    const double alphaRhoL = u.alphaL()*u.densityL();

    data[ALPHA_RHO_G]        = alphaRhoG;
    data[ALPHA_RHO_U_G]      = alphaRhoG*u.velocityUG();
    data[ALPHA_RHO_V_G]      = alphaRhoG*u.velocityVG();
    data[ALPHA_RHO_W_G]      = alphaRhoG*u.velocityWG();
    data[ALPHA_RHO_E_PINT_G] = u.alphaG()*(u.densityG()*u.totalEnergyG() + pInt);

    data[ALPHA_RHO_L]        = alphaRhoL;
    data[ALPHA_RHO_U_L]      = alphaRhoL*u.velocityUL();
    data[ALPHA_RHO_V_L]      = alphaRhoL*u.velocityVL();
    data[ALPHA_RHO_W_L]      = alphaRhoL*u.velocityWL();
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