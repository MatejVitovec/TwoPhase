#include "BoundaryCondition.hpp"


BoundaryCondition::BoundaryConditionType BoundaryCondition::getType() const
{
    return type;
}

Boundary BoundaryCondition::getBoundary() const
{
    return boundary;
}

void BoundaryCondition::updateMeshBoundary(const Mesh& mesh)
{
    const std::vector<Boundary>& boundaryList = mesh.getBoundaryList();
    for (int i = 0; i < boundaryList.size(); i++)
    {
        if(boundary.boundaryConditionName == boundaryList[i].boundaryConditionName)
        {
            boundary = boundaryList[i];
            break;
        }
    }
}

void BoundaryCondition::updateId(int id_)
{
    id = id_;
}

std::vector<Compressible> BoundaryCondition::calc(const VolField<Compressible>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<Compressible> out(boundary.facesIndex.size());

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
       out[i] = calculateState(w[ownerIndexList[boundary.facesIndex[i]]], thermoField[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], thermoModel);
    }
    
    return out;
}

void BoundaryCondition::correct(const Field<Compressible>& w, Field<Compressible>& wl, Field<Compressible>& wr, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    
}


std::vector<CompressibleMixture> BoundaryCondition::calc(const VolField<CompressibleMixture>& w, const VolField<ThermoVar>& thermoField, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<CompressibleMixture> out(boundary.facesIndex.size());

    for (int i = 0; i < boundary.facesIndex.size(); i++)
    {
       out[i] = calculateState(w[ownerIndexList[boundary.facesIndex[i]]], thermoField[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], thermoModel);
    }
    
    return out;
}

void BoundaryCondition::correct(const Field<CompressibleMixture>& w, Field<CompressibleMixture>& wl, Field<CompressibleMixture>& wr, const Field<Mat<9,3>>& grad, const Field<Vars<9>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    
}

void BoundaryCondition::correct(const Field<Primitive>& u, Field<Primitive>& ul, Field<Primitive>& ur,
                                const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi,
                                const Mesh& mesh, const Thermo * const thermoModel) const
{

}

void BoundaryCondition::apply(VolField<Fluid>& u, const Mesh& mesh, const Fluid * const thermoModel) const
{

}

void BoundaryCondition::correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{

}



void BoundaryCondition::apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{

}

void BoundaryCondition::correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{

}