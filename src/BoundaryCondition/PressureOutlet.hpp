#ifndef PRESSUREOUTLET_HPP
#define PRESSUREOUTLET_HPP

#include "BoundaryCondition.hpp"

class PressureOutlet : public BoundaryCondition
{
    public:

        PressureOutlet(Boundary meshBoundary, double pressure_, int id_) : BoundaryCondition(meshBoundary, PRESSUREOUTLET, id_), pressure(pressure_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        

    private:
        double pressure;
};

#endif // PRESSUREOUTLET_HPP