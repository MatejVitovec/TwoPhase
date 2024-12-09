#ifndef LEASTSQUARE_HPP
#define LEASTSQUARE_HPP

#include "GradientScheme.hpp"

class LeastSquare : public GradientScheme
{
    public:

        LeastSquare() {}

        virtual ~LeastSquare() {}

        void init(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);

        Field<Mat<5,3>> calculateGradient(const VolField<Compressible>& w, const Mesh& mesh) const;

    private:
        Field<Mat<3,3>> MInv;
        Field<Vars<3>> cellToCellDelta; //faceField

        void calculateInverseM(Field<Mat<3,3>> M);
        void calculateCellToCellDelta(const Mesh& mesh, const std::vector<std::shared_ptr<BoundaryCondition>>& boundaryConditionList);
};

#endif // LEASTSQUARE_HPP