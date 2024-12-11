#ifndef IDEALGASTHERMO_HPP
#define IDEALGASTHERMO_HPP

#include "Thermo.hpp"

class IdealGasThermo : public Thermo
{
    public:

        IdealGasThermo() : Thermo(), gamma(1.3326), R(461.51805), energyShift(1995917.326) {}
        IdealGasThermo(double gamma_, double R_) : Thermo(), gamma(gamma_), R(R_), energyShift(1995917.326) {}

        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        Vars<3> updateThermo(const Primitive& data, const PrimitiveThermoVar& thermoData) const;
        void updateThermo(ComponentThermoVar& thermoData) const;

        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

    private:
        double gamma;
        double R;
        double energyShift;

        double p(double T, double rho) const;
        double rho(double p, double T) const;
        double T(double rho, double p) const;
        double a(double rho, double p) const;
        double e(double p, double T) const;
        
        double pFromRho_e(double rho, double e) const;

};

#endif // IDEALGASTHERMO_HPP