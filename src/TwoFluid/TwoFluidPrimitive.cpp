#include <cmath>
#include "TwoFluidPrimitive.hpp"

double TwoFluidPrimitive::alpha() const
{
    return data[ALPHA];
}

double TwoFluidPrimitive::pressure() const
{
    return data[P];
}

double TwoFluidPrimitive::temperatureG() const
{
    return data[T_G];
}

double TwoFluidPrimitive::temperatureL() const
{
    return data[T_L];
}

double TwoFluidPrimitive::alphaG() const
{
    return data[ALPHA];
}

double TwoFluidPrimitive::alphaL() const
{
    return 1.0 - data[ALPHA];
}

Vars<3> TwoFluidPrimitive::velocityG() const
{
    return Vars<3>({data[U_G], data[V_G], data[W_G]});
}

double TwoFluidPrimitive::absVelocityG() const
{
    return std::sqrt(this->absVelocity2G());
}

double TwoFluidPrimitive::absVelocity2G() const
{
    return data[U_G]*data[U_G] + data[V_G]*data[V_G] + data[W_G]*data[W_G];
}

double TwoFluidPrimitive::normalVelocityG(const Vars<3>& normalVector) const
{
    return dot(this->velocityG(), normalVector);
}

double TwoFluidPrimitive::velocityUG() const
{
    return data[U_G];
}

double TwoFluidPrimitive::velocityVG() const
{
    return data[V_G];
}

double TwoFluidPrimitive::velocityWG() const
{
    return data[W_G];
}



Vars<3> TwoFluidPrimitive::velocityL() const
{
    return Vars<3>({data[U_L], data[V_L], data[W_L]});
}

double TwoFluidPrimitive::absVelocityL() const
{
    return std::sqrt(this->absVelocity2L());
}

double TwoFluidPrimitive::absVelocity2L() const
{
    return data[U_L]*data[U_L] + data[V_L]*data[V_L] + data[W_L]*data[W_L];
}

double TwoFluidPrimitive::normalVelocityL(const Vars<3>& normalVector) const
{
    return dot(this->velocityL(), normalVector);
}

double TwoFluidPrimitive::velocityUL() const
{
    return data[U_L];
}

double TwoFluidPrimitive::velocityVL() const
{
    return data[V_L];
}

double TwoFluidPrimitive::velocityWL() const
{
    return data[W_L];
}