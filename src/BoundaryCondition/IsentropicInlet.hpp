#ifndef ISENTROPICINLET_HPP
#define ISENTROPICINLET_HPP

#include "BoundaryCondition.hpp"

class IsentropicInlet : public BoundaryCondition
{
    public:

        IsentropicInlet(Boundary meshBoundary,
                        double totalDensity_,
                        double totalPressure_,
                        double totalTemperature_,
                        double totalEntropy_,
                        double totalEnthalpy_,
                        Vars<3> velocityDirection_,
                        int id_) : BoundaryCondition(meshBoundary, ISENTROPICINLET, id_),
                                                      totalDensity(totalDensity_),
                                                      totalPressure(totalPressure_),
                                                      totalTemperature(totalTemperature_),
                                                      totalEntropy(totalEntropy_),
                                                      totalEnthalpy(totalEnthalpy_),
                                                      velocityDirection(velocityDirection_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;

        //TOWFLUID
        void apply(VolField<Fluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const;
        void correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const; 

        void apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const;
        void correct(const VolField<TwoFluid>& u, Field<TwoFluid>& ul, Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const; 


    private:
        double totalDensity;
        double totalPressure;
        double totalTemperature;
        double totalEntropy;
        double totalEnthalpy;
        Vars<3> velocityDirection;
};

#endif // ISENTROPICINLET_HPP
