#ifndef TWOFLUIDTHERMO_HPP
#define TWOFLUIDTHERMO_HPP

#include "../TwoFluid.hpp"
#include "../TwoFluidCompressible.hpp"
#include "../../Field.hpp"
#include "../../VolField.hpp"
#include "../../Mesh/Mesh.hpp"

#include "../../Thermo/StiffenedGasThermo.hpp"

#include "../../Thermo/StateEquations/NonLinearSolver/NewtonMethod.hpp"

#include <chrono>

class TwoFluidThermo
{
    public:
    
        //TwoFluidThermo() : gasThermo(StiffenedGasThermo()), liquidThermo(StiffenedGasThermo()), newton(NewtonMethod()) {}
        TwoFluidThermo() : newton(NewtonMethod()) {}

        virtual ~TwoFluidThermo() {}

        void updateFromConservative(VolField<TwoFluid>& u, const Field<TwoFluidCompressible>& w, const Field<double>& pInt) const;
        
        void update(VolField<TwoFluid>& u) const;
        void update(Field<TwoFluid>& u) const;
        void updateInternal(Field<TwoFluid>& u) const;
        void updateBoundary(VolField<TwoFluid>& u) const;

        TwoFluid isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, const Vars<3>& velocityDirection, const TwoFluid& stateIn) const;

    private:
        //StiffenedGasThermo gasThermo;
        //StiffenedGasThermo liquidThermo;

        NewtonMethod newton;

        const double gammaG = 1.4;
        const double gammaL = 2.8;

        const double cpG = 1004.5;
        const double cpL = 4186.0;

        const double pInfG = 0.0;
        const double pInfL = 8.5e8;
};

#endif // TWOFLUIDTHERMO_HPP