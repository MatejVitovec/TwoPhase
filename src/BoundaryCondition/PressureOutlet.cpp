#include "PressureOutlet.hpp"


Compressible PressureOutlet::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    if(w.normalVelocity(f.normalVector)/thermoVar.soundSpeed() >= 1.0)
    //if(w.absVelocity()/w.soundSpeed() >= 1.0)
    {
        return w;
    }

    return thermoModel->primitiveToConservative(Vars<5>({w.density(),
                                                         w.velocityU(),
                                                         w.velocityV(),
                                                         w.velocityW(),
                                                         pressure}));
}


CompressibleMixture PressureOutlet::calculateState(const CompressibleMixture& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    if(w.normalVelocity(f.normalVector)/thermoVar.soundSpeed() >= 1.0)
    //if(w.absVelocity()/w.soundSpeed() >= 1.0)
    {
        return w;
    }

    CompressibleMixture out = w;

    out.updateCompressiblePart(thermoModel->primitiveToConservative(Vars<5>({w.density(),
                                                                             w.velocityU(),
                                                                             w.velocityV(),
                                                                             w.velocityW(),
                                                                             pressure})));

    return out;
}




void PressureOutlet::apply(VolField<Fluid>& u, const Mesh& mesh, const Thermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<Fluid>& boundaryData = u.boundary(id);

    for (int i = 0; i < boundary.facesIndex.size(); i++) //Predelat - zbavid se boundary.facesIndex a nahradit indexem boundaryId a boundaryFaceList nahrat z meshe
    {
        updateState(u[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], boundaryData[i]);    
    }

    //Predat info o updatu, pokud neni update tlaku neni nutny update thermo
    //thermoModel->update(u); //TODO
}

void PressureOutlet::updateState(const Fluid& stateIn, const Face& f, Fluid& u) const
{
    u = stateIn;

    if(u.normalVelocity(f.normalVector)/u.soundSpeed() >= 1.0)
    //if(w.absVelocity()/w.soundSpeed() >= 1.0)
    {
        return;
    }

    u[Fluid::P] = pressure;
}

void PressureOutlet::correct(const VolField<Fluid>& u, const Field<Fluid>& ul, const Field<Fluid>& ur, const Field<Mat<5,3>>& grad, const Field<Vars<5>>& phi, const Mesh& mesh, const Thermo * const thermoModel) const
{
    //DO NOTHING
    return;
}



void PressureOutlet::apply(VolField<TwoFluid>& u, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{
    const std::vector<Face>& faceList = mesh.getFaceList();
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();

    std::vector<TwoFluid>& boundaryData = u.boundary(id);

    for (int i = 0; i < boundary.facesIndex.size(); i++) //Predelat - zbavid se boundary.facesIndex a nahradit indexem boundaryId a boundaryFaceList nahrat z meshe
    {
        updateState(u[ownerIndexList[boundary.facesIndex[i]]], faceList[boundary.facesIndex[i]], boundaryData[i]);    
    }

    //Predat info o updatu, pokud neni update tlaku neni nutny update thermo
    //thermoModel->update(u); //TODO
}

void PressureOutlet::updateState(const TwoFluid& stateIn, const Face& f, TwoFluid& u) const
{
    u = stateIn;

    const double commonNormalMach = 0.5*(u.normalVelocityG(f.normalVector)/u.soundSpeedG() + u.normalVelocityL(f.normalVector)/u.soundSpeedL());

    if(commonNormalMach >= 1.0)
    {
        return;
    }

    u[Fluid::P] = pressure;
}



void PressureOutlet::correct(const VolField<TwoFluid>& u, Field<TwoFluid>& ul, Field<TwoFluid>& ur, const Field<Mat<10,3>>& grad, const Field<Vars<10>>& phi, const Mesh& mesh, const TwoFluidThermo * const thermoModel) const
{
    //DO NOTHING
    return;
}

