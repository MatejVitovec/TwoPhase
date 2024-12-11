#ifndef THERMOVAR_HPP
#define THERMOVAR_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"

class ThermoVar : public Vars<3>
{
    public:
        enum {T, P, A};

        ThermoVar() : Vars<3>() {}
        ThermoVar(const Vars<3>& in) : Vars<3>(in) {}
        ThermoVar(const std::array<double, 3>& in) : Vars<3>(in) {}

        virtual ~ThermoVar() {}

        double temperature() const;
        double pressure() const;        
        double soundSpeed() const;

};


#endif // THERMOVAR_HPP