#include "Slau2.hpp"

Vars<5> Slau2::calculateFlux(const double rhoL, const double rhoR,
                            const Vars<3> uL, const Vars<3> uR,
                            const double pL, const double pR,
                            const double eL, const double eR,
                            const double aL, const double aR,
                            const Vars<3>& normalVector) const
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

    return Vars<5>({0.5*(mTilde + magMTilde)       + 0.5*(mTilde - magMTilde),
                    0.5*(mTilde + magMTilde)*uL[0] + 0.5*(mTilde - magMTilde)*uR[0] + pTilde*normalVector[0],
                    0.5*(mTilde + magMTilde)*uL[1] + 0.5*(mTilde - magMTilde)*uR[1] + pTilde*normalVector[1],
                    0.5*(mTilde + magMTilde)*uL[2] + 0.5*(mTilde - magMTilde)*uR[2] + pTilde*normalVector[2],
                    0.5*(mTilde + magMTilde)*HLeft + 0.5*(mTilde - magMTilde)*HRight});

}


Vars<5> Slau2::calculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    return calculateFlux(wl.density(), wr.density(), wl.velocity(), wr.velocity(), thermoL.pressure(), thermoR.pressure(), wl.internalEnergy(), wr.internalEnergy(), thermoL.soundSpeed(), thermoR.soundSpeed(), normalVector);
}

Vars<9> Slau2::calculateFlux(const CompressibleMixture& wl, const CompressibleMixture& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    return Vars<9>(); //TODO
}
