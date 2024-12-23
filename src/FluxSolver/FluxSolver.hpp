#ifndef FLUXSOLVER_HPP
#define FLUXSOLVER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../Compressible.hpp"
#include "../ThermoVar.hpp"
#include "../Primitive.hpp"
#include "../PrimitiveThermoVar.hpp"

#include "../Mixture/CompressibleMixture.hpp"

class FluxSolver
{
    public:

        FluxSolver() {}

        virtual ~FluxSolver() {}

        virtual Vars<5> calculateFlux(const double rhoL, const double rhoR,
                                      const Vars<3> uL, const Vars<3> uR,
                                      const double pL, const double pR,
                                      const double eL, const double eR,
                                      const double aL, const double aR,
                                      const Vars<3>& normalVector) const = 0;

        Field<Vars<5>> calculateFluxes(const Field<Compressible>& wl, const Field<Compressible>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const;
        virtual Vars<5> calculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const = 0;
        
        Field<Vars<5>> calculateFluxes(const Field<Primitive>& wl, const Field<Primitive>& wr, const Field<PrimitiveThermoVar>& thermoFieldL, const Field<PrimitiveThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const;

        Field<Vars<9>> calculateFluxes(const Field<CompressibleMixture>& wl, const Field<CompressibleMixture>& wr, const Field<ThermoVar>& thermoFieldL, const Field<ThermoVar>& thermoFieldR, const std::vector<Face>& faceList) const;
        virtual Vars<9> calculateFlux(const CompressibleMixture& wl, const CompressibleMixture& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const = 0;

    protected:

};

#endif // FLUXSOLVER_HPP