#include "Inlet.hpp"


void Inlet::apply(VolField<Fluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<Fluid>& boundaryData = u.boundary(id);

    for (int i = 0; i < boundary.facesIndex.size(); i++) //Predelat - zbavid se boundary.facesIndex a nahradit indexem boundaryId a boundaryFaceList nahrat z meshe
    {
        boundaryData[i] = Fluid({u[ownerIndexList[boundary.facesIndex[i]]].pressure(), velocity[0], velocity[1], velocity[2], temperature});
    }    
}


void Inlet::correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    return; //DO NOTHING
}


//TWOFLUID
void Inlet::apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{
    constexpr double epsilon = 1.0e-7; //TODO - asi jako global

    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<TwoFluid>& boundaryData = u.boundary(id);

    for (int i = 0; i < boundary.facesIndex.size(); i++) //Predelat - zbavid se boundary.facesIndex a nahradit indexem boundaryId a boundaryFaceList nahrat z meshe
    {
        boundaryData[i] = TwoFluid({epsilon, u[ownerIndexList[boundary.facesIndex[i]]].pressure(), temperature, temperature, velocity[0], velocity[1], velocity[2], velocity[0], velocity[1], velocity[2]}); 
    }    
}

void Inlet::correct(const VolField<TwoFluid>& u, const Field<TwoFluid>& ul, const Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{
    return; //DO NOTHING
}