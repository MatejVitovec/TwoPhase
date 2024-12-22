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

        virtual void updateFromConservative(VolField<TwoFluid>& u, const Field<TwoFluidCompressible>& w, const Field<double>& pInt) const;
        
        virtual void update(VolField<TwoFluid>& u) const;
        virtual void update(Field<TwoFluid>& u) const;
        virtual void updateBoundary(VolField<TwoFluid>& u) const;

    private:
        //StiffenedGasThermo gasThermo;
        //StiffenedGasThermo liquidThermo;

        NewtonMethod newton;

        const double gammaG = 1.0;
        const double gammaL = 1.0;

        const double cpG = 1.0;
        const double cpL = 1.0;

        const double pInfG = 1.0;
        const double pInfL = 1.0;
};

#endif // TWOFLUIDTHERMO_HPP