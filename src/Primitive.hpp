#ifndef PRIMITIVE_HPP
#define PRIMITIVE_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"
#include "Compressible.hpp"
#include "ThermoVar.hpp"
#include "PrimitiveThermoVar.hpp"

class Primitive : public Vars<5>
{
    public:
        using Vars<5>::operator+=;
        using Vars<5>::operator-=;

        enum {P, U, V, W, T};

        Primitive() : Vars<5>() {}
        Primitive(const Vars<5>& varsIn) : Vars<5>(varsIn) {}
        Primitive(const std::array<double, 5>& in) : Vars<5>(in) {}

        Primitive(const Compressible& w, const ThermoVar& thermo) : Vars<5>({thermo.pressure(), w.velocityU(), w.velocityV(), w.velocityW(), thermo.temperature()}) {}

        virtual ~Primitive() {}

        double pressure() const;
        Vars<3> velocity() const;
        double absVelocity() const;
        double absVelocity2() const;
        double normalVelocity(const Vars<3>& normalVector) const;
        double velocityU() const;
        double velocityV() const;
        double velocityW() const;
        double temperature() const;

        Compressible toCompressible(const PrimitiveThermoVar& thermoVar) const;
};

#endif // PRIMITIVE_HPP