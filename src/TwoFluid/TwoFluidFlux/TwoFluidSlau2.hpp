#ifndef TWOFLUIDSLAU2_HPP
#define TWOFLUIDSLAU2_HPP

#include "../../Mesh/Mesh.hpp"
#include "../../Field.hpp"
#include "../TwoFluid.hpp"
#include "TwoFluidFluxSolver.hpp"


class TwoFluidSlau2 : public TwoFluidFluxSolver
{
    public:

        TwoFluidSlau2() {}

        virtual ~TwoFluidSlau2() {}

        void calculateFlux(const double alphaL, const double alphaR,
            const double rhoL, const double rhoR,
            const Vars<3> uL, const Vars<3> uR,
            const double pL, const double pR,
            const double eL, const double eR,
            const double aL, const double aR,
            const Vars<3>& normalVector,
            Vars<5>& fluxL, Vars<5>& fluxR);

        void calculateFluxes(const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const std::vector<Face>& faceList, Field<Vars<10>>& fluxesl, Field<Vars<10>>& fluxesr);

    private:

};

#endif // TWOFLUIDSLAU2_HPP