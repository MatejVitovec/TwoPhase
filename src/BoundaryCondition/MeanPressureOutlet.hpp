#ifndef MEANPRESSUREOUTLET_HPP
#define MEANPRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class MeanPressureOutlet : public BoundaryCondition
{
    public:

        MeanPressureOutlet(Boundary meshBoundary, double pressure_, int id_) : BoundaryCondition(meshBoundary, MEANPRESSUREOUTLET, id_), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        std::vector<Compressible> calc(const VolField<Compressible>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const;

        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        std::vector<CompressibleMixture> calc(const VolField<CompressibleMixture>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const;

        

    private:
        double pressure;

        double calculateCorrectionConstant(const Mesh& mesh, const Field<Compressible>& w, const Field<ThermoVar>& thermoField) const;

        double calculateCorrectionConstant(const Mesh& mesh, const Field<CompressibleMixture>& w, const Field<ThermoVar>& thermoField) const;
};



#endif // MEANPRESSUREOUTLET_HPP