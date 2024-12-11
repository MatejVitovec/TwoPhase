#ifndef TWOFLUIDFLUXSOLVER_HPP
#define TWOFLUIDFLUXSOLVER_HPP

#include "../../Mesh/Mesh.hpp"
#include "../../Field.hpp"
#include "../TwoFluid.hpp"
#include "../../FluxSolver/FluxSolver.hpp"


class TwoFluidFluxSolver
{
    public:

        TwoFluidFluxSolver() {}

        virtual ~TwoFluidFluxSolver() {}

        Field<Vars<10>> calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList) const;
        

    protected:
        std::shared_ptr<FluxSolver> GGFluxSolver;
        std::shared_ptr<FluxSolver> LLFluxSolver;
        std::shared_ptr<FluxSolver> GLFluxSolver;

};

#endif // TWOFLUIDFLUXSOLVER_HPP