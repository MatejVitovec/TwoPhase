#ifndef THERMO_HPP
#define THERMO_HPP

#include "../Compressible.hpp"
#include "../ThermoVar.hpp"
#include "../Primitive.hpp"
#include "../PrimitiveThermoVar.hpp"
#include "../ComponentThermoVar.hpp"
#include "../TwoFluid/TwoFluid.hpp"
#include "../TwoFluid/TwoFluidCompressible.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Mesh/Mesh.hpp"

#include <chrono>

class Thermo
{
    public:
    
        Thermo() {}

        virtual ~Thermo() {}

        void updateThermo(VolField<ThermoVar>& thermoField, const VolField<Compressible>& w) const;
        void updateThermoInternal(VolField<ThermoVar>& thermoField, const VolField<Compressible>& w) const;
        void updateThermo(Field<ThermoVar>& thermoField, const Field<Compressible>& w) const;

        void updateThermo(Field<PrimitiveThermoVar>& thermoField, const Field<Primitive>& u) const;

        void updateThermo(VolField<ComponentThermoVar>& thermoField) const;
        void updateThermoInternal(VolField<ComponentThermoVar>& thermoField) const;
        void updateThermo(Field<ComponentThermoVar>& thermoField) const;

        virtual Vars<3> updateThermo(const Compressible& data, const ThermoVar& thermoData) const = 0;
        virtual Vars<3> updateThermo(const Primitive& data, const PrimitiveThermoVar& thermoData) const = 0;
        virtual void updateThermo(ComponentThermoVar& thermoData) const = 0;

        virtual Compressible primitiveToConservative(const Vars<5>& primitive) const = 0;
        virtual Compressible stagnationState(double TTot, double pTot) const = 0;

        virtual std::array<double, 3> initPressureTemperatureInlet(double pTot, double TTot) const = 0;

        virtual Compressible isentropicInlet(double pTot, double TTot, double rhoTot, double sTot, double hTot, Vars<3> velocityDirection, Compressible stateIn, ThermoVar thermoIn) const = 0;

        ////// TWOFLUID
        virtual void updateFromConservative(VolField<TwoFluid>& u, const Field<TwoFluidCompressible>& w) const;
        virtual void update(Field<TwoFluid>& u) const;

        
};

#endif // THERMO_HPP