#include <cmath>

#include "TwoFluidFluxSolver.hpp"

Field<Vars<10>> TwoFluidFluxSolver::calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> GGFlux(ul.size());
    Field<Vars<5>> LLFlux(ul.size());
    Field<Vars<5>> GLFlux(ul.size());
    Field<Vars<5>> LGFlux(ul.size());

    Field<Vars<10>> out(ul.size());

    #pragma omp parallel for
    for (int i = 0; i < ul.size(); i++)
    {
        //out[i] = claculateFlux(wl[i], wr[i], thermoFieldL[i], thermoFieldR[i], faceList[i].normalVector)*faceList[i].area;
    }


    for (int i = 0; i < ul.size(); i++)
    {
        const TwoFluid& uL = ul[i];
        const TwoFluid& uR = ur[i];
        out[i] += join(GGFlux[i]*std::min(uL.alphaG(), uR.alphaG()), LLFlux[i]*std::min(uL.alphaL(), uR.alphaL()));
    }
    

    return out;




}

