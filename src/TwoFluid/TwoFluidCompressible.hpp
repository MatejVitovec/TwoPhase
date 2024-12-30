#ifndef TWOFLUIDCOMPRESSIBLE_HPP
#define TWOFLUIDCOMPRESSIBLE_HPP

#include <vector>
#include <memory>

#include "../Vars.hpp"
#include "TwoFluid.hpp"

class TwoFluidCompressible : public Vars<10>
{
    public:
        using Vars<10>::operator+=;
        using Vars<10>::operator-=;

        enum {ALPHA_RHO_G, ALPHA_RHO_U_G, ALPHA_RHO_V_G, ALPHA_RHO_W_G, ALPHA_RHO_E_PINT_G, ALPHA_RHO_L, ALPHA_RHO_U_L, ALPHA_RHO_V_L, ALPHA_RHO_W_L, ALPHA_RHO_E_PINT_L};

        TwoFluidCompressible() : Vars<10>() {}
        TwoFluidCompressible(const Vars<10>& varsIn) : Vars<10>(varsIn) {}
        TwoFluidCompressible(const std::array<double, 10>& in) : Vars<10>(in) {}
        TwoFluidCompressible(const TwoFluid& u);
        //TwoFluidCompressible(const std::array<double, 5>& inRef ) : 

        virtual ~TwoFluidCompressible() {}

        void updateInterfacialPressure(const double pInt, const TwoFluid& u);

        double alphaDensityG() const;
        double densityG(double alpha);
        Vars<3> velocityG() const;
        double absVelocityG() const;
        double absVelocity2G() const;
        double normalVelocityG(const Vars<3>& normalVector) const;
        double velocityUG() const;
        double velocityVG() const;
        double velocityWG() const;

        double alphaDensityL() const;
        double densityL(double alpha);
        Vars<3> velocityL() const;
        double absVelocityL() const;
        double absVelocity2L() const;
        double normalVelocityL(const Vars<3>& normalVector) const;
        double velocityUL() const;
        double velocityVL() const;
        double velocityWL() const;

        //Vars<10> flux(const Vars<3>& thermoData, const Vars<3>& normalVector) const;
};

#endif // TWOFLUIDCOMPRESSIBLE_HPP