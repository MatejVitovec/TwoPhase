#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "../Mesh/Mesh.hpp"

#include "../Mesh/Boundary.hpp"
#include "../Mesh/Face.hpp"
#include "../Compressible.hpp"
#include "../Fluid.hpp"
#include "../ThermoVar.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"


#include "../Mat.hpp"
#include "../Thermo/Thermo.hpp"

#include "../Mixture/CompressibleMixture.hpp"
#include "../TwoFluid/TwoFluid.hpp"
#include "../TwoFluid/TwoFluidThermo/TwoFluidThermo.hpp"

class BoundaryCondition
{
    public:
        enum BoundaryConditionType{ISENTROPICINLET, PRESSUREOUTLET, FREEBOUNDARY, WALL, PERIODICITY, MEANPRESSUREOUTLET, INLET};

        BoundaryCondition(BoundaryConditionType type_) : type(type_) , id(0) {}
        BoundaryCondition(Boundary meshBoundary, BoundaryConditionType type_, int id_) : boundary(meshBoundary), type(type_), id(id_) {}

        BoundaryConditionType getType() const;
        Boundary getBoundary() const;

        void updateMeshBoundary(const Mesh& mesh);

        void updateId(int id_);

        virtual Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const = 0;
        virtual CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const = 0;


        virtual std::vector<Compressible> calc(const VolField<Compressible>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const;
        virtual std::vector<CompressibleMixture> calc(const VolField<CompressibleMixture>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const;

        virtual void correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr,
                             const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi,
                             const Mesh& mesh, const Thermo * const thermoModel) const;

        virtual void correct(const Field<Primitive>& u, Field<Primitive>& ul, Field<Primitive>& ur,
                             const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi,
                             const Mesh& mesh, const Thermo * const thermoModel) const;

        virtual void correct(const Field<CompressibleMixture>& w, Field<CompressibleMixture>& wl, Field<CompressibleMixture>& wr,
                             const Field<Mat<9,3>>& grad, const Field<Vars<9>>& phi,
                             const Mesh& mesh, const Thermo * const thermoModel) const;
        
        //Jina implementace pro TWOFLUID kdyz se osvetsi predelat i pro zbytek
        virtual void apply(VolField<Fluid>& u, const Mesh& mesh, const Fluid * const thermoModel) const;
        virtual void correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const; 


        virtual void apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const;
        virtual void correct(const VolField<TwoFluid>& u, Field<TwoFluid>& ul, Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const; 


    protected:
        Boundary boundary;
        BoundaryConditionType type;
        int id;
};

#endif // BOUNDARYCONDITION_HPP