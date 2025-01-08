#ifndef FLUID_HPP
#define FLUID_HPP

#include <vector>
#include <memory>

#include "Primitive.hpp"

class Fluid : public Primitive
{
    public:
        using Vars<5>::operator+=;
        using Vars<5>::operator-=;

        enum {RHO, INT_E, A};

        Fluid() : Primitive(), thermoData() {}
        Fluid(const Vars<5>& varsIn) : Primitive(varsIn), thermoData() {}
        Fluid(const std::array<double, 5>& in) : Primitive(in), thermoData() {}

        virtual ~Fluid() {}

        double density() const;
        double internalEnergy() const;
        double totalEnergy() const;
        double soundSpeed() const;

        std::array<double, 3>& getThermoData();
        const std::array<double, 3>& getThermoData() const;
    
    protected:
        std::array<double, 3> thermoData;
};

#endif // FLUID_HPP