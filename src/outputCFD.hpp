#ifndef OUTPUTCFD_HPP
#define OUTPUTCFD_HPP

#include "Mesh/Mesh.hpp"
#include "Field.hpp"
#include "Compressible.hpp"
#include "Condensation/CompressibleMixture.hpp"
#include "ThermoVar.hpp"
#include "Mat.hpp"
#include <string>

namespace outputCFD
{
    void outputVTK(std::string fileName, const Mesh& mesh, const Field<Compressible>& w, const Field<ThermoVar>& thermoField);
    void outputVTK(std::string fileName, const Mesh& mesh, const Field<CompressibleMixture>& w, const Field<ThermoVar>& thermoField);

    void outputVTKPeriodicBoundary(std::string fileName, const Mesh& m, const Field<Compressible>& w, const Field<ThermoVar>& thermoField, Vars<3> shift);

    void saveData(std::string fileName, const Field<Compressible>& w);

    void saveResidual(std::string fileName, Vars<5> res);
    void saveResidual(std::string fileName, double res);
    void saveValue(std::string fileName, double val);

    void saveFieldOnBoundary(std::string fileName, std::string boundaryName, const Mesh& mesh, const Field<Compressible>& w, const Field<ThermoVar>& thermoField);

    void saveLimiters(Field<Vars<5>> phi, const Mesh& mesh);
    void saveGradients(Field<Mat<5,3>> grad, const Mesh& mesh);

    Field<Compressible> loadCompressibleFieldFromVTK(std::string fileName);
}

#endif //OUTPUTCFD_HPP