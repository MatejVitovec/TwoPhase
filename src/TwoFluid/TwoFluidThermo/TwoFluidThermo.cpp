#include <cmath>
#include <omp.h>

#include "TwoFluidThermo.hpp"

void TwoFluidThermo::updateFromConservative(VolField<TwoFluid>& u, const Field<TwoFluidCompressible>& w, const Field<double>& pInt) const
{
    for (int i = 0; i < w.size(); i++)
    {
        const double guessAlpha = u[i].alpha();
        const double guessP = u[i].pressure();

        const Vars<3> uG = w[i].velocityG();
        const Vars<3> uL = w[i].velocityL();

        const double alphaRhoG = w[i][TwoFluidCompressible::ALPHA_RHO_G];
        const double alphaRhoL = w[i][TwoFluidCompressible::ALPHA_RHO_L];

        const double epsilonG = w[i][TwoFluidCompressible::ALPHA_RHO_E_PINT_G];
        const double epsilonL = w[i][TwoFluidCompressible::ALPHA_RHO_E_PINT_L];

        const double thetaG = 1.0;
        const double thetaL = 1.0;

        const double AG = (gammaG - 1.0)*pInt[i];
        const double AL = gammaL*pInfL + (gammaL - 1.0)*pInt[i];
        const double BG = (gammaG - 1.0)*(thetaG - epsilonG);
        const double BL = (gammaL - 1.0)*(thetaL - epsilonL);

        std::pair<double, double> alphaP = newton.solve([=](double val1, double val2) { return val1*(val2 + AG) + BG; },
                                                        [=](double val1, double val2) { return val2 + AG; },
                                                        [=](double val1, double val2) { return val1; },
                                                        [=](double val1, double val2) { return val1*(val2 + AL) + BL; },
                                                        [=](double val1, double val2) { return -val2 - AL; },
                                                        [=](double val1, double val2) { return 1.0 - val1; },
                                                        guessAlpha,
                                                        guessP);

        const double alphaG = alphaP.first;
        const double alphaL = 1.0 - alphaP.first;
        const double pressure = alphaP.second;

        const double rhoG = alphaRhoG/alphaG;
        const double rhoL = alphaRhoL/alphaL;

        u[i][TwoFluid::ALPHA] = alphaG;
        u[i][TwoFluid::P] = pressure;
        
        u[i][TwoFluid::RHO_G] = rhoG;
        u[i][TwoFluid::RHO_L] = rhoL;

        u[i][TwoFluid::INT_E_G] = epsilonG/alphaRhoG - 0.5*norm2sqr(uG) - pInt[i]/rhoG;
        u[i][TwoFluid::INT_E_L] = epsilonL/alphaRhoL - 0.5*norm2sqr(uL) - pInt[i]/rhoL;

        //////NATVRDO ZAPSANE FUNKCE - PREDELAT
        u[i][TwoFluid::T_G] = (gammaG*u[i][TwoFluid::INT_E_G])/cpG;
        u[i][TwoFluid::T_L] = (gammaL*u[i][TwoFluid::INT_E_L]*(pressure + pInfL))/(cpL*(pressure + gammaL*pInfL));
        u[i][TwoFluid::A_G] = std::sqrt((gammaG*pressure)/rhoG);
        u[i][TwoFluid::A_L] = std::sqrt((gammaL*(pressure + pInfL))/rhoL);
        //////

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
    update(u);
    updateBoundary(u);
}

void TwoFluidThermo::update(Field<TwoFluid>& u) const
{
    for (int i = 0; i < u.size(); i++)
    {
        const double pressure = u[i][TwoFluid::P];

        u[i][TwoFluid::RHO_G] = (pressure*gammaG)/((gammaG - 1.0)*cpG*u[i][TwoFluid::T_G]);
        u[i][TwoFluid::RHO_L] = ((pressure + pInfL)*gammaL)/((gammaL - 1.0)*cpL*u[i][TwoFluid::T_L]);
        u[i][TwoFluid::INT_E_G] = (cpG/gammaG)*u[i][TwoFluid::T_G];
        u[i][TwoFluid::INT_E_L] = (cpL/gammaL)*u[i][TwoFluid::T_L] + pInfL/u[i][TwoFluid::RHO_L];
        u[i][TwoFluid::A_G] = std::sqrt((gammaG*pressure)/u[i][TwoFluid::RHO_G]);
        u[i][TwoFluid::A_L] = std::sqrt((gammaL*(pressure + pInfL))/u[i][TwoFluid::RHO_L]);
    }
}


void TwoFluidThermo::updateBoundary(VolField<TwoFluid>& u) const
{
    for (int id = 0; id < u.boundarySize(); id++)
    {
        for (int i = 0; i < u.boundary(id).size(); i++)
        {
            const double pressure = u.boundary(id)[i][TwoFluid::P];

            u.boundary(id)[i][TwoFluid::RHO_G] = (pressure*gammaG)/((gammaG - 1.0)*cpG*u.boundary(id)[i][TwoFluid::T_G]);
            u.boundary(id)[i][TwoFluid::RHO_L] = ((pressure + pInfL)*gammaL)/((gammaL - 1.0)*cpL*u.boundary(id)[i][TwoFluid::T_L]);
            u.boundary(id)[i][TwoFluid::INT_E_G] = (cpG/gammaG)*u.boundary(id)[i][TwoFluid::T_G];
            u.boundary(id)[i][TwoFluid::INT_E_L] = (cpL/gammaL)*u.boundary(id)[i][TwoFluid::T_L] + pInfL/u.boundary(id)[i][TwoFluid::RHO_L];
            u.boundary(id)[i][TwoFluid::A_G] = std::sqrt((gammaG*pressure)/u.boundary(id)[i][TwoFluid::RHO_G]);
            u.boundary(id)[i][TwoFluid::A_L] = std::sqrt((gammaL*(pressure + pInfL))/u.boundary(id)[i][TwoFluid::RHO_L]);
        }
    }    
}


