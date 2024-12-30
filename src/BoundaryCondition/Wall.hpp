#ifndef WALL_HPP
#define WALL_HPP

#include "BoundaryCondition.hpp"

class Wall : public BoundaryCondition
{
    public:

        Wall(Boundary meshBoundary, int id_) : BoundaryCondition(meshBoundary, WALL, id_) {}

        Compressible calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        void correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;
        void correct(const Field<Primitive>& u, Field<Primitive>& ul, Field<Primitive>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;


        CompressibleMixture calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const;
        void correct(const Field<CompressibleMixture>& w, Field<CompressibleMixture>& wl, Field<CompressibleMixture>& wr, const Field<Mat<9,3>>& grad, const Field<Vars<9>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const;

        //TWOFLUID
        void apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const;
        void correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const; 

        void updateState(const TwoFluid& steteIn, const Face& f, TwoFluid& u) const;

    private:

};

#endif // WALL_HPP