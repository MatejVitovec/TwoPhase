#ifndef TWOFLUIDPRIMITIVE_HPP
#define TWOFLUIDPRIMITIVE_HPP

#include <vector>
#include <memory>

#include "../Vars.hpp"

class TwoFluidPrimitive : public Vars<10>
{
    public:
        using Vars<10>::operator+=;
        using Vars<10>::operator-=;

        enum {ALPHA, P, T_G, T_L, U_G, V_G, W_G, U_L, V_L, W_L};

        TwoFluidPrimitive() : Vars<10>() {}
        TwoFluidPrimitive(const Vars<10>& varsIn) : Vars<10>(varsIn) {}
        TwoFluidPrimitive(const std::array<double, 10>& in) : Vars<10>(in) {}
        //TwoFluidPrimitive(const std::array<double, 5>& inRef ) : 

        virtual ~TwoFluidPrimitive() {}

        double alpha() const;
        double pressure() const;        
        double temperatureG() const;
        double temperatureL() const;

        double alphaG() const;
        double alphaL() const;

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
        double velocityWL() const;
};

#endif // TWOFLUIDPRIMITIVE_HPP