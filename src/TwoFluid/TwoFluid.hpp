#ifndef TWOFLUID_HPP
#define TWOFLUID_HPP

#include <vector>
#include <memory>

#include "../Vars.hpp"
#include "TwoFluidPrimitive.hpp"

class TwoFluid : public TwoFluidPrimitive
{
    public:
        //using Vars<10>::operator+=;
        //using Vars<10>::operator-=;

        //enum {ALPHA, P, T_G, T_L, U_G, V_G, W_G, U_L, V_L, W_L};
        enum {RHO_G, RHO_L, INT_E_G, INT_E_L, A_G, A_L};

        TwoFluid() : TwoFluidPrimitive(), thermoData() {}
        TwoFluid(const Vars<10>& varsIn) : TwoFluidPrimitive(varsIn), thermoData() {}
        TwoFluid(const TwoFluidPrimitive& in) : TwoFluidPrimitive(in), thermoData() {}
        TwoFluid(const std::array<double, 10>& in) : TwoFluidPrimitive(in), thermoData() {}

        virtual ~TwoFluid() {}

        /*double alpha() const;
        double pressure() const;        
        double temperatureG() const;
        double temperatureL() const;

        Vars<3> velocityG() const;
        double absVelocityG() const;
        double absVelocity2G() const;
        double normalVelocityG(const Vars<3>& normalVector) const;
        double velocityUG() const;
        double velocityVG() const;
        double velocityWG() const;

        Vars<3> velocityL() const;
        double absVelocityL() const;
        double absVelocity2L() const;
        double normalVelocityL(const Vars<3>& normalVector) const;
        double velocityUL() const;
        double velocityVL() const;
        double velocityWL() const;*/

        double densityG() const;
        double densityL() const;
        double internalEnergyG() const;
        double internalEnergyL() const;
        double soundSpeedG() const;
        double soundSpeedL() const;

        double interfacialPressure() const;
        void blend();

        std::array<double, 6>& getThermoData();
        const std::array<double, 6>& getThermoData() const;

    private:
        std::array<double, 6> thermoData;
};

#endif // TWOFLUID_HPP