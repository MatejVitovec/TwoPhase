#include <cmath>
#include "Primitive.hpp"

double Primitive::pressure() const
{
    return data[P];
}

Vars<3> Primitive::velocity() const
{
    return Vars<3>({data[U], data[V], data[W]});
}

double Primitive::absVelocity() const
{
    return sqrt(this->absVelocity2());
}

double Primitive::absVelocity2() const
{
    return data[U]*data[U] + data[V]*data[V] + data[W]*data[W];
}

double Primitive::normalVelocity(const Vars<3>& normalVector) const
{
    return dot(this->velocity(), normalVector);
}

double Primitive::velocityU() const
{
    return data[U];
}

double Primitive::velocityV() const
{
    return data[V];
}

double Primitive::velocityW() const
{
    return data[W];
}

double Primitive::temperature() const
{
    return data[T];
}

Compressible Primitive::toCompressible(const PrimitiveThermoVar& thermoVar) const
{
    const double rho = thermoVar.density();
    return Compressible({rho, rho*data[U], rho*data[V], rho*data[W], rho*(thermoVar.internalEnergy() + 0.5*norm2sqr(this->velocity()))});
}