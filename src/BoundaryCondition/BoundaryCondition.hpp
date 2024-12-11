#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "../Mesh/Mesh.hpp"

#include "../Mesh/Boundary.hpp"
#include "../Mesh/Face.hpp"
#include "../Compressible.hpp"
#include "../Primitive.hpp"
#include "../ThermoVar.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../TwoFluid/TwoFluid.hpp"

#include "../Mat.hpp"
#include "../Thermo/Thermo.hpp"

#include "../Condensation/CompressibleMixture.hpp"

class BoundaryCondition
{
    public:
        enum BoundaryConditionType{ISENTROPICINLET, PRESSUREOUTLET, FREEBOUNDARY, WALL, PERIODICITY, MEANPRESSUREOUTLET};

        BoundaryCondition(BoundaryConditionType type_) : type(type_) {}
        BoundaryCondition(Boundary meshBoundary, BoundaryConditionType type_) : boundary(meshBoundary), type(type_) {}

        BoundaryConditionType getType() const;
        Boundary getBoundary() const;

        void updateMeshBoundary(const Mesh& mesh);

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
        virtual void apply(const VolField<TwoFluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const;
        virtual void correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const; 


    protected:
        Boundary boundary;
        BoundaryConditionType type;
};

#endif // BOUNDARYCONDITION_HPP