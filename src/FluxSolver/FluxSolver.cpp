#include <cmath>

#include "FluxSolver.hpp"

Field<Vars<5>> FluxSolver::calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(wl.size());

    #pragma omp parallel for
    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], thermoFieldL[i], thermoFieldR[i], faceList[i].normalVector)*faceList[i].area;
    }
    return out;
}

Field<Vars<5>> FluxSolver::calculateFluxes(const Field<Primitive>& ul, const Field<Primitive>& ur, const Field<PrimitiveThermoVar>& thermoFieldL, const Field<PrimitiveThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const
{
    Field<Vars<5>> out(ul.size());

    #pragma omp parallel for
    for (int i = 0; i < ul.size(); i++)
    {
        out[i] = claculateFlux(thermoFieldL[i].density(),        thermoFieldR[i].density(),
                               ul[i].velocity(),                 ur[i].velocity(),
                               ul[i].pressure(),                 ur[i].pressure(),
                               thermoFieldL[i].internalEnergy(), thermoFieldR[i].internalEnergy(),
                               thermoFieldL[i].soundSpeed(),     thermoFieldR[i].soundSpeed(),
                               faceList[i].normalVector)*faceList[i].area;
    }
    
    return out;
}


Field<Vars<9>> FluxSolver::calculateFluxes(const Field<CompressibleMixture>& wl, const Field<CompressibleMixture>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const
{
    Field<Vars<9>> out(wl.size());

    #pragma omp parallel for
    for (int i = 0; i < wl.size(); i++)
    {
        out[i] = claculateFlux(wl[i], wr[i], thermoFieldL[i], thermoFieldR[i], faceList[i].normalVector)*faceList[i].area;
    }
    
    return out;
}
