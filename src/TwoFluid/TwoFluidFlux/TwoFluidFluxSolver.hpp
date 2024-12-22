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

        virtual Vars<10> calculateFlux(const TwoFluid& ul, const TwoFluid& ur, const Vars<3>& normalVector) const = 0;
        

    protected:

};

#endif // TWOFLUIDFLUXSOLVER_HPP