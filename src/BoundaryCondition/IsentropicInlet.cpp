#include "IsentropicInlet.hpp"

Compressible IsentropicInlet::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w, thermoVar);
}

CompressibleMixture IsentropicInlet::calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    return thermoModel->isentropicInlet(totalPressure, totalTemperature, totalDensity, totalEntropy, totalEnthalpy, velocityDirection, w, thermoVar);
}


//TWOFLUID
void IsentropicInlet::apply(VolField<TwoFluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const
{

}

void IsentropicInlet::correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{

}
