#ifndef TwoFluidFlux_HPP
#define TwoFluidFlux_HPP

#include "../../Mesh/Mesh.hpp"
#include "../../Field.hpp"
#include "../TwoFluid.hpp"
#include "../../FluxSolver/FluxSolver.hpp"


class TwoFluidFlux
{
    public:

        TwoFluidFlux() {}

        virtual ~TwoFluidFlux() {}

        Field<Vars<10>> calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList) const;

        void calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList, Field<Vars<10>>& fluxesl, Field<Vars<10>>& fluxesr);
        

    protected:
        std::shared_ptr<FluxSolver> GGFluxSolver;
        std::shared_ptr<FluxSolver> LLFluxSolver;
        std::shared_ptr<FluxSolver> GLFluxSolver;

};

#endif // TwoFluidFlux_HPP