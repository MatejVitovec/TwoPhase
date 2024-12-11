#ifndef PRIMITIVETHERMOVAR_HPP
#define PRIMITIVETHERMOVAR_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"
#include "Compressible.hpp"
#include "ThermoVar.hpp"

class PrimitiveThermoVar : public Vars<3>
{
    public:
        enum {RHO, INT_E, A};

        PrimitiveThermoVar() : Vars<3>() {}
        PrimitiveThermoVar(const Vars<3>& in) : Vars<3>(in) {}
        PrimitiveThermoVar(const std::array<double, 3>& in) : Vars<3>(in) {}

        PrimitiveThermoVar(const Compressible& w, const ThermoVar& thermo) : Vars<3>({w.density(), w.internalEnergy(), thermo.soundSpeed()}) {}

        virtual ~PrimitiveThermoVar() {}

        double density() const;
        double internalEnergy() const;
        double soundSpeed() const;
};


#endif // PRIMITIVETHERMOVAR_HPP