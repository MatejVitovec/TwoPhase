#include <cmath>

#include "TwoFluidSlau2.hpp"


void TwoFluidSlau2::calculateFlux(const double alphaL, const double alphaR,
                                     const double rhoL, const double rhoR,
                                     const Vars<3> uL, const Vars<3> uR,
                                     const double pL, const double pR,
                                     const double eL, const double eR,
                                     const double aL, const double aR,
                                     const Vars<3>& normalVector,
                                     Vars<5>& fluxL, Vars<5>& fluxR)
{
    const double EL = eL + 0.5*norm2sqr(uL);
    const double ER = eR + 0.5*norm2sqr(uR);

    const double HLeft  = eL + 0.5*norm2sqr(uL)  + pL/rhoL;
    const double HRight = eR + 0.5*norm2sqr(uR) + pR/rhoR;

    const double qLeft  = dot(uL, normalVector);
    const double qRight = dot(uR, normalVector);

    const double aTilde   = 0.5*(aL   + aR);
    const double rhoTilde = 0.5*(rhoL + rhoR);

    const double sqrtUDash = std::sqrt(0.5*(qLeft*qLeft + qRight*qRight));

    const double MaRelLeft  = qLeft /aTilde;
    const double MaRelRight = qRight/aTilde;

    const double Chi = std::pow(1.0 - std::min(1.0, (1.0/aTilde)*sqrtUDash), 2.0);

    const double g = -std::max(std::min(MaRelLeft, 0.0), -1.0)*std::min(std::max(MaRelRight, 0.0), 1.0);

    const double magVnBar = (rhoL*std::abs(qLeft) + rhoR*std::abs(qRight))
        /(rhoL + rhoR);
        
    const double magVnBarPlus  = (1.0 - g)*magVnBar + g*std::abs(qLeft);
    const double magVnBarMinus = (1.0 - g)*magVnBar + g*std::abs(qRight);

    const double PPlusLeft   = ((std::abs(MaRelLeft)  >= 1.0) ?
        0.5*(1.0 + sgn(MaRelLeft))
        : (0.25*std::pow(MaRelLeft + 1.0, 2.0)*(2.0 - MaRelLeft)));
    const double PMinusRight = ((std::abs(MaRelRight) >= 1.0) ?
        0.5*(1.0 - sgn(MaRelRight))
        : (0.25*std::pow(MaRelRight - 1.0, 2.0)*(2.0 + MaRelRight)));

    const double pTilde = 0.5*(pL + pR)
        + 0.5*(PPlusLeft - PMinusRight)*(pL - pR)
        + sqrtUDash*(PPlusLeft + PMinusRight - 1.0)*rhoTilde*aTilde;


    const double mTilde = 0.5*(rhoL*(qLeft + magVnBarPlus)
        + rhoR*(qRight - magVnBarMinus)
        - (Chi/aTilde)*(pR - pL));

    const double magMTilde = std::abs(mTilde);

    fluxL = Vars<5>({0.5*(mTilde + magMTilde)*alphaL       + 0.5*(mTilde - magMTilde)*alphaR,
                     0.5*(mTilde + magMTilde)*alphaL*uL[0] + 0.5*(mTilde - magMTilde)*alphaR*uR[0],
                     0.5*(mTilde + magMTilde)*alphaL*uL[1] + 0.5*(mTilde - magMTilde)*alphaR*uR[1],
                     0.5*(mTilde + magMTilde)*alphaL*uL[2] + 0.5*(mTilde - magMTilde)*alphaR*uR[2],
                     0.5*(mTilde + magMTilde)*alphaL*HLeft + 0.5*(mTilde - magMTilde)*alphaR*HRight});

    fluxR = fluxL;

    fluxL[1] += alphaL*pTilde*normalVector[0];
    fluxL[2] += alphaL*pTilde*normalVector[1];
    fluxL[3] += alphaL*pTilde*normalVector[2];

    fluxR[1] += alphaR*pTilde*normalVector[0];
    fluxR[2] += alphaR*pTilde*normalVector[1];
    fluxR[3] += alphaR*pTilde*normalVector[2];
}




void TwoFluidSlau2::calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList, Field<Vars<10>>& fluxesl, Field<Vars<10>>& fluxesr)
{
    for (int i = 0; i < ul.size(); i++)
    {
        const Vars<3>& normalVector = faceList[i].normalVector;
        const TwoFluid& uL = ul[i];
        const TwoFluid& uR = ur[i];

        Vars<5> fluxLG;
        Vars<5> fluxRG;
        Vars<5> fluxLL;
        Vars<5> fluxRL;

        double al = uL.soundSpeedG() + uL.soundSpeedL();
        double ar = uR.soundSpeedG() + uR.soundSpeedL();

        calculateFlux(uL.alphaG(), uR.alphaG(),
                      uL.densityG(), uR.densityG(),
                      uL.velocityG(), uR.velocityG(),
                      uL.pressure(), uR.pressure(),
                      uL.internalEnergyG(), uR.internalEnergyG(),
                      al, ar,
                      normalVector,
                      fluxLG, fluxRG);

        calculateFlux(uL.alphaL(), uR.alphaL(),
                      uL.densityL(), uR.densityL(),
                      uL.velocityL(), uR.velocityL(),
                      uL.pressure(), uR.pressure(),
                      uL.internalEnergyL(), uR.internalEnergyL(),
                      al, ar,
                      normalVector,
                      fluxLL, fluxRL);

        fluxesl[i] = join(fluxLG, fluxLL);
    }
}