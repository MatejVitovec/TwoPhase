#include <cmath>

#include "TwoFluidHllc.hpp"

Vars<10> TwoFluidHllc::calculateFlux(const TwoFluid& ul, const TwoFluid& ur, const Vars<3>& normalVector) const
{
    double uG;
    double uL;
    double rhoG;
    double rhoL;
    double pG;
    double pL;
    double aG;
    double aL;

    double sl;
    double sr;
    
    if (ul.alphaG() < ur.alphaG())
    {
        uG = dot(ur.velocityG(), normalVector);
        uL = dot(ul.velocityL(), normalVector);
        rhoG = ur.densityG();
        rhoL = ul.densityL();
        pG = ur.pressure();
        pL = ul.pressure();
        aG = ur.soundSpeedG();
        aL = ul.soundSpeedL();
    }
    else
    {
        uG = dot(ul.velocityG(), normalVector);
        uL = dot(ur.velocityL(), normalVector);
        rhoG = ul.densityG();
        rhoL = ur.densityL();
        pG = ul.pressure();
        pL = ur.pressure();
        aG = ul.soundSpeedG();
        aL = ur.soundSpeedL();
    }
    
    double bG = std::min(uG + aG, 0.0);
    double bL = std::min(uL - aL, 0.0);


    ///TODOOOOO
}