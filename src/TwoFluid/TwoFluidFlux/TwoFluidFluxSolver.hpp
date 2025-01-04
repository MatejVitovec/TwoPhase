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

        virtual void calculateFlux(const double alphaL, const double alphaR,
            const double rhoL, const double rhoR,
            const Vars<3>& uL, const Vars<3>& uR,
            const double pL, const double pR,
            const double eL, const double eR,
            const double aL, const double aR,
            const Vars<3>& normalVector,
            Vars<5>& fluxL, Vars<5>& fluxR) = 0;
        

        virtual void calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList, Field<Vars<10>>& fluxesl, Field<Vars<10>>& fluxesr) = 0;

    protected:

};

#endif // TWOFLUIDFLUXSOLVER_HPP