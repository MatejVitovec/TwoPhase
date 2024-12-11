#ifndef STIFFENEDGASTHERMO_HPP
#define STIFFENEDGASTHERMO_HPP

#include "Thermo.hpp"

class StiffenedGasThermo : public Thermo
{
    public:

        StiffenedGasThermo() : Thermo(), gamma(1.4), cp(1000.0), pInf(0.0) {}
        StiffenedGasThermo(double gamma_, double cp_, double pInf_) : Thermo(), gamma(gamma_), cp(cp_), pInf(pInf_) {}

        ////OLD
        Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const;
        Vars<3> updateThermo(const Primitive& data, const PrimitiveThermoVar& thermoData) const;
        void updateThermo(ComponentThermoVar& thermoData) const;

        Compressible primitiveToConservative(const Vars<5>& primitive) const;
        Compressible stagnationState(double TTot, double pTot) const;

        Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const;
        std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const;

        ///OLD

    private:
        double gamma;
        double cp;
        double pInf;

        double p(double T, double rho) const;
        double rho(double p, double T) const;
        double T(double rho, double p) const;
        double a(double rho, double p) const;
        double e(double p, double T) const;
        
        double pFromRho_e(double rho, double e) const;

};

#endif // STIFFENEDGASTHERMO_HPP