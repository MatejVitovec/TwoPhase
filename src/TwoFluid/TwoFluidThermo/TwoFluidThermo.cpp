#include <cmath>
#include <omp.h>

#include "TwoFluidThermo.hpp"

void TwoFluidThermo::updateFromConservative(VolField<TwoFluid>& u, const Field<TwoFluidCompressible>& w, const Field<double>& pInt) const
{
    #pragma omp parallel for
    for (int i = 0; i < w.size(); i++)
    {       
        const double guessAlpha = u[i].alpha();
        const double guessPressure = u[i].pressure();

        const Vars<3> uG = w[i].velocityG();
        const Vars<3> uL = w[i].velocityL();

        const double alphaRhoG = w[i][TwoFluidCompressible::ALPHA_RHO_G];
        const double alphaRhoL = w[i][TwoFluidCompressible::ALPHA_RHO_L];

        const double epsilonG = w[i][TwoFluidCompressible::ALPHA_RHO_E_PINT_G];
        const double epsilonL = w[i][TwoFluidCompressible::ALPHA_RHO_E_PINT_L];

        const double thetaG = 0.5*alphaRhoG*norm2sqr(uG);
        const double thetaL = 0.5*alphaRhoL*norm2sqr(uL);

        const double AG = (gammaG - 1.0)*pInt[i];
        const double AL = (gammaL - 1.0)*pInt[i] + gammaL*pInfL;
        const double BG = (gammaG - 1.0)*(thetaG - epsilonG);
        const double BL = (gammaL - 1.0)*(thetaL - epsilonL);

        std::pair<double, double> alphaP = newton.solve([=](double alpha, double pressure) { return alpha*(pressure + AG) + BG; },
                                                        [=](double alpha, double pressure) { return pressure + AG; },
                                                        [=](double alpha, double pressure) { return alpha; },
                                                        [=](double alpha, double pressure) { return (1.0 - alpha)*(pressure + AL) + BL; },
                                                        [=](double alpha, double pressure) { return -(pressure + AL); },
                                                        [=](double alpha, double pressure) { return 1.0 - alpha; },
                                                        guessAlpha,
                                                        guessPressure);

        const double alphaG = alphaP.first;
        const double alphaL = 1.0 - alphaG;
        const double pressure = alphaP.second;

        const double rhoG = alphaRhoG/alphaG;
        const double rhoL = alphaRhoL/alphaL;

        u[i][TwoFluid::ALPHA] = alphaG;
        u[i][TwoFluid::P] = pressure;

        //////NATVRDO ZAPSANE FUNKCE - PREDELAT
        u[i][TwoFluid::T_G] = (gammaG/cpG)*(epsilonG/alphaRhoG - 0.5*norm2sqr(uG) - pInt[i]/rhoG);
        u[i][TwoFluid::T_L] = (gammaL/cpL)*(epsilonL/alphaRhoL - 0.5*norm2sqr(uL) - pInt[i]/rhoL - pInfL/rhoL);

        u[i][TwoFluid::U_G] = uG[0];
        u[i][TwoFluid::V_G] = uG[1];
        u[i][TwoFluid::W_G] = uG[2];
        u[i][TwoFluid::U_L] = uL[0];
        u[i][TwoFluid::V_L] = uL[1];
        u[i][TwoFluid::W_L] = uL[2];
    }
}

void TwoFluidThermo::update(VolField<TwoFluid>& u) const
{
    updateInternal(u);
    updateBoundary(u);
}

void TwoFluidThermo::update(Field<TwoFluid>& u) const
{
    updateInternal(u);
}

void TwoFluidThermo::updateInternal(Field<TwoFluid>& u) const
{
    #pragma omp parallel for
    for (int i = 0; i < u.size(); i++)
    {
        const double pressure = u[i][TwoFluid::P];
        const double temperatureG = u[i][TwoFluid::T_G];
        const double temperatureL = u[i][TwoFluid::T_L];

        auto& thermoDataRef = u[i].getThermoData();

        thermoDataRef[TwoFluid::RHO_G] = (pressure*gammaG)/((gammaG - 1.0)*cpG*temperatureG);
        thermoDataRef[TwoFluid::RHO_L] = ((pressure + pInfL)*gammaL)/((gammaL - 1.0)*cpL*temperatureL);
        thermoDataRef[TwoFluid::INT_E_G] = (cpG/gammaG)*temperatureG;
        thermoDataRef[TwoFluid::INT_E_L] = (cpL/gammaL)*temperatureL + pInfL/thermoDataRef[TwoFluid::RHO_L];
        thermoDataRef[TwoFluid::A_G] = std::sqrt((gammaG*pressure)/thermoDataRef[TwoFluid::RHO_G]);
        thermoDataRef[TwoFluid::A_L] = std::sqrt((gammaL*(pressure + pInfL))/thermoDataRef[TwoFluid::RHO_L]);
    }
}


void TwoFluidThermo::updateBoundary(VolField<TwoFluid>& u) const
{
    for (int id = 0; id < u.boundarySize(); id++)
    {
        for (int i = 0; i < u.boundary(id).size(); i++)
        {
            const double pressure = u.boundary(id)[i][TwoFluid::P];
            const double temperatureG = u.boundary(id)[i][TwoFluid::T_G];
            const double temperatureL = u.boundary(id)[i][TwoFluid::T_L];

            auto& thermoDataRef = u.boundary(id)[i].getThermoData();

            thermoDataRef[TwoFluid::RHO_G] = (pressure*gammaG)/((gammaG - 1.0)*cpG*temperatureG);
            thermoDataRef[TwoFluid::RHO_L] = ((pressure + pInfL)*gammaL)/((gammaL - 1.0)*cpL*temperatureL);
            thermoDataRef[TwoFluid::INT_E_G] = (cpG/gammaG)*temperatureG;
            thermoDataRef[TwoFluid::INT_E_L] = (cpL/gammaL)*temperatureL + pInfL/thermoDataRef[TwoFluid::RHO_L];
            thermoDataRef[TwoFluid::A_G] = std::sqrt((gammaG*pressure)/thermoDataRef[TwoFluid::RHO_G]);
            thermoDataRef[TwoFluid::A_L] = std::sqrt((gammaL*(pressure + pInfL))/thermoDataRef[TwoFluid::RHO_L]);
        }
    }    
}


TwoFluid TwoFluidThermo::isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, const Vars<3>& velocityDirection, const TwoFluid& stateIn) const
{
    //TODO
    TwoFluid out = TwoFluid();

    return out;
}