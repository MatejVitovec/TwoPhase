#include <cmath>

#include "TwoFluidFluxSolver.hpp"

Field<Vars<10>> TwoFluidFluxSolver::calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList) const
{
    Field<Vars<10>> out(ul.size());

    #pragma omp parallel for
    for (int i = 0; i < ul.size(); i++)
    {
        const TwoFluid& stateL = ul[i];
        const TwoFluid& stateR = ur[i];

        const double thetaGG = std::min(stateL.alphaG(), stateR.alphaG());
        const double thetaGL = std::max(0.0, stateR.alphaL() - stateL.alphaL());
        const double thetaLG = std::max(0.0, stateR.alphaG() - stateL.alphaG());
        const double thetaLL = std::min(stateL.alphaL(), stateR.alphaL());

        Vars<5> ggFlux = Vars<5>();
        Vars<5> llFlux = Vars<5>();
        Vars<10> glFlux = Vars<10>();
        Vars<10> lgFlux = Vars<10>();

        if (thetaGG > 0)
        {
            ggFlux = GGFluxSolver->claculateFlux(stateL.densityG(),         stateR.densityG(),
                                                 stateL.velocityG(),        stateR.velocityG(),
                                                 stateL.pressure(),         stateR.pressure(),
                                                 stateL.internalEnergyG(),  stateR.internalEnergyG(),
                                                 stateL.soundSpeedG(),      stateR.soundSpeedG(),
                                                 faceList[i].normalVector);
        }

        if (thetaLL > 0)
        {
            llFlux = LLFluxSolver->claculateFlux(stateL.densityL(),         stateR.densityL(),
                                                 stateL.velocityL(),        stateR.velocityL(),
                                                 stateL.pressure(),         stateR.pressure(),
                                                 stateL.internalEnergyL(),  stateR.internalEnergyL(),
                                                 stateL.soundSpeedL(),      stateR.soundSpeedL(),
                                                 faceList[i].normalVector);
        }

        if (thetaGL > 0)
        {
            glFlux = Vars<10>(); //TODO
        }
        else if (thetaLG > 0)
        {
            lgFlux = Vars<10>(); //TODO
        }
        
        out[i] = (join(thetaGG*ggFlux, thetaLL*llFlux) + thetaGL*glFlux + thetaLG*lgFlux)*faceList[i].area;
    }

    return out;
}

