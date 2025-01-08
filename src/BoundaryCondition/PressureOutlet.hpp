#ifndef PRESSUREOUTLET_HPP
#define PRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class PressureOutlet : public BoundaryCondition
{
    public:

        PressureOutlet(Boundary meshBoundary, double pressure_, int id_) : BoundaryCondition(meshBoundary, PRESSUREOUTLET, id_), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;

    //NOVY PRISTUP

        void apply(VolField<Fluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const;
        void correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const; 

        void updateState(const Fluid& stateIn, const Face& f, Fluid& u) const;

        void apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const;
        void correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const; 

        void updateState(const TwoFluid& stateIn, const Face& f, TwoFluid& u) const;
        

    private:
        double pressure;
};

#endif // PRESSUREOUTLET_HPP