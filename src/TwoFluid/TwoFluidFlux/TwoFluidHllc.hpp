#ifndef TWOFLUIDHLLC_HPP
#define TWOFLUIDHLLC_HPP

#include "../../Mesh/Mesh.hpp"
#include "../../Field.hpp"
#include "../TwoFluid.hpp"
#include "TwoFluidFluxSolver.hpp"


class TwoFluidHllc : public TwoFluidFluxSolver
{
    public:

        TwoFluidHllc() {}

        virtual ~TwoFluidHllc() {}

        Vars<10> calculateFlux(const TwoFluid& ul, const TwoFluid& ur, const Vars<3>& normalVector) const;

    private:

};

#endif // TWOFLUIDHLLC_HPP