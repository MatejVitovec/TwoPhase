#ifndef INLET_HPP
#define INLET_HPP

#include "BoundaryCondition.hpp"

class Inlet : public BoundaryCondition
{
    public:

        Inlet(Boundary meshBoundary, double temperature_, Vars<3> velocity_, int id_) : BoundaryCondition(meshBoundary, INLET, id_), temperature(temperature_), velocity(velocity_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const { return Compressible(); }
        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const { return Compressible(); }

        //TWOFLUID
        void apply(VolField<Fluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const;
        void correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const; 

        void apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const;
        void correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const; 

    private:
        double temperature;
        Vars<3> velocity;
};

#endif // INLET_HPP