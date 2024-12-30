#ifndef TWOFLUID_HPP
#define TWOFLUID_HPP

#include <vector>
#include <memory>

#include "../Vars.hpp"
#include "TwoFluidPrimitive.hpp"

class TwoFluid : public TwoFluidPrimitive
{
    public:
        using Vars<10>::operator+=;
        using Vars<10>::operator-=;

        enum {RHO_G, RHO_L, INT_E_G, INT_E_L, A_G, A_L};

        TwoFluid() : TwoFluidPrimitive(), thermoData() {}
        TwoFluid(const Vars<10>& varsIn) : TwoFluidPrimitive(varsIn), thermoData() {}
        TwoFluid(const TwoFluidPrimitive& in) : TwoFluidPrimitive(in), thermoData() {}
        TwoFluid(const std::array<double, 10>& in) : TwoFluidPrimitive(in), thermoData() {}

        virtual ~TwoFluid() {}

        double densityG() const;
        double densityL() const;
        double internalEnergyG() const;
        double internalEnergyL() const;
        double soundSpeedG() const;
        double soundSpeedL() const;

        double totalEnergyG() const;
        double totalEnergyL() const;

        double interfacialPressure() const;
        void blend();

        std::array<double, 6>& getThermoData();
        const std::array<double, 6>& getThermoData() const;

    private:

        static constexpr double epsilonMin = 0.1e-8;
        static constexpr double epsilonMax = 0.1e-4;

        std::array<double, 6> thermoData;
};

#endif // TWOFLUID_HPP