#include <set>
#include <iostream>
#include <algorithm>
#include "TwoFluidFVMScheme.hpp"

#include "BoundaryCondition/Periodicity.hpp"

#include "outputCFD.hpp"

void TwoFluidFVMScheme::setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_)
{
    gradientScheme = std::move(gradScheme_);
}
void TwoFluidFVMScheme::setReconstructionLimiter(std::unique_ptr<Limiter> limiter_)
{
    limiter = std::move(limiter_);
}

void TwoFluidFVMScheme::setCfl(double cfl_)
{
    cfl = cfl_;
}

void TwoFluidFVMScheme::setMaxIter(int maxIter_)
{
    maxIter = maxIter_;
}

void TwoFluidFVMScheme::setSaveEveryIter(int saveEveryIter_)
{
    saveEveryIter = saveEveryIter_;
}

void TwoFluidFVMScheme::setTargetError(double targetError_)
{
    targetError = targetError_;
}

void TwoFluidFVMScheme::setLocalTimeStep(bool localTimeStep_)
{
    localTimeStep = localTimeStep_;
}

void TwoFluidFVMScheme::setReconstructionSettings(bool reconstruction_)
{
    reconstruction = reconstruction_;
}

double TwoFluidFVMScheme::getCfl() const
{
    return cfl;
}

int TwoFluidFVMScheme::getMaxIter() const
{
    return maxIter;
}

double TwoFluidFVMScheme::getTargetError() const
{
    return targetError;
}

bool TwoFluidFVMScheme::getTimeStepsettings() const
{
    return localTimeStep;
}

bool TwoFluidFVMScheme::getReconstructionSettings() const
{
    return reconstruction;
}

const Mesh& TwoFluidFVMScheme::getMesh() const
{
    return mesh;
}

const Thermo* TwoFluidFVMScheme::getThermoRef()
{
    return thermo.get();
}


void TwoFluidFVMScheme::setInitialConditions(TwoFluidPrimitive initialCondition)
{
    for (int i = 0; i < mesh.getCellsSize(); i++)
    {
        u[i] = initialCondition;
    }

    //thermo->updateThermoInternal(thermoField, w);  //TODO
}

void TwoFluidFVMScheme::init()
{
    applyFreeBoundaryCondition(); // presenout do CaseSetter, pote mouno zredukovat tuto init() funkci

    gradientScheme->init(mesh, boundaryConditionList);

    //w.getBoundaryData() = std::vector<std::vector<Compressible>>(boundaryConditionList.size());
    u.getBoundaryData() = std::vector<std::vector<TwoFluid>>(boundaryConditionList.size());

    ul = Field<TwoFluid>(mesh.getFacesSize());
    ur = Field<TwoFluid>(mesh.getFacesSize());

    grad = Field<Mat<10,3>>(mesh.getCellsSize());
    phi = Field<Vars<10>>(mesh.getCellsSize());

    for (int boundaryId = 0; boundaryId < boundaryConditionList.size(); boundaryId++)
    {
        int size = boundaryConditionList[boundaryId]->getBoundary().facesIndex.size();
        u.boundary(boundaryId) = std::vector<TwoFluid>(size);
    }    

    timeSteps = Field<double>(mesh.getCellsSize());

    fixGradient = 1000000;
    fixedGradStep = 10000000;
}

void TwoFluidFVMScheme::applyFreeBoundaryCondition()
{
    for(auto& boundaryCondition : boundaryConditionList)
    {
        if(boundaryCondition->getType() == BoundaryCondition::FREEBOUNDARY)
        {
            mesh.deleteBoundary(boundaryCondition->getBoundary().boundaryConditionName);
        }
    }

    boundaryConditionList.erase(std::remove_if(boundaryConditionList.begin(), boundaryConditionList.end(),
        [](const std::shared_ptr<BoundaryCondition> & boundary) { return (boundary->getType() == BoundaryCondition::FREEBOUNDARY); }),
        boundaryConditionList.end());

    for(auto& boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->updateMeshBoundary(mesh);
    }

    for(auto& boundaryCondition : boundaryConditionList)
    {
        if(boundaryCondition->getType() == BoundaryCondition::PERIODICITY)
        {
            static_cast<Periodicity*>(boundaryCondition.get())->init(mesh);
        }
    }
}

void TwoFluidFVMScheme::setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions)
{
    boundaryConditionList = std::move(boundaryConditions);
}

void TwoFluidFVMScheme::interpolateToFaces()
{
    const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
    const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

    if(reconstruction)
    {
        for (int boundaryId = 0; boundaryId < u.getBoundaryData().size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;
            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                ur[boundaryFacesIndexList[i]] = u.boundary(boundaryId)[i];
            }
        }

        if (iter < fixGradient || iter % fixedGradStep == 0)
        {
            grad = gradientScheme->calculateGradient(u, mesh);
            phi = limiter->calculateLimiter(u, grad, mesh);
        }

        const std::vector<Cell>& cells = mesh.getCellList();
        const std::vector<Face>& faces = mesh.getFaceList();
        const std::vector<int>& ownerIndexList = mesh.getOwnerIndexList();
        const std::vector<int>& neighborIndexList = mesh.getNeighborIndexList();

        for (int i = 0; i < faces.size(); i++)
        {
            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                Vars<10> ulDiff = dot(grad[ownerIndexList[i]], faces[i].midpoint - cells[ownerIndexList[i]].center);
                Vars<10> urDiff = dot(grad[neighborIndexList[i]], faces[i].midpoint - cells[neighborIndexList[i]].center);

                ul[i] = u[ownerIndexList[i]] + phi[ownerIndexList[i]]*ulDiff;
                ur[i] = u[neighborIndexList[i]] + phi[neighborIndexList[i]]*urDiff;
            }
            else
            {
                ul[i] = u[ownerIndexList[i]];
            }
        }

        for (auto & boundaryCondition : boundaryConditionList)
        {
            boundaryCondition->correct(u, ul, ur, grad, phi, mesh, thermo.get());
        }
        
        thermo->update(ul);
        thermo->update(ur);
    }
    else
    {
        for (int i = 0; i < mesh.getFacesSize(); i++)
        {        
            ul[i] = u[ownerIndexList[i]];

            int neighbour = neighborIndexList[i];
            if(neighbour >= 0)
            {
                ur[i] = u[neighbour];
            }
        }

        for (int boundaryId = 0; boundaryId < u.getBoundaryData().size(); boundaryId++)
        {
            const std::vector<int>& boundaryFacesIndexList = boundaryConditionList[boundaryId]->getBoundary().facesIndex;

            for (int i = 0; i < boundaryFacesIndexList.size(); i++)
            {
                ur[boundaryFacesIndexList[i]] = u.boundary(boundaryId)[i];
            }
        }
    }
}


void TwoFluidFVMScheme::applyBoundaryConditions() 
{
    for(auto& boundaryCondition : boundaryConditionList)
    {
        boundaryCondition->apply(u, mesh, thermo.get());
    }
}

void TwoFluidFVMScheme::boundField()
{
    /*constexpr double minDensity = 0.1;
    constexpr double maxMach = 2.2;
    constexpr double minInternalEnergy = 2100000.0;
    constexpr double maxInternalEnergy = 2800000.0;
    

    const std::vector<Cell>& cells = mesh.getCellList();
    for (int i = 0; i < cells.size(); i++)
    {
        if (w[i].density() < minDensity)
        {
            //std::cout << "rho: " << w[i].density() << " U: " << w[i].absVelocity() << " energy: " << w[i].internalEnergy() << std::endl;
            w[i][Compressible::RHO_U] = w[i].velocityU()*minDensity;
            w[i][Compressible::RHO_V] = w[i].velocityV()*minDensity;
            w[i][Compressible::RHO_W] = w[i].velocityW()*minDensity;
            w[i][Compressible::RHO_E] = w[i].totalEnergy()*minDensity;
            w[i][Compressible::RHO] = minDensity;
        }

        //if (w[i].machNumber() > maxMach)
        //{
        //    double velocityCoeff = (maxMach*w[i].soundSpeed())/(w[i].absVelocity());
        //    w[i][Compressible::RHO_U] = w[i][Compressible::RHO_U]*velocityCoeff;
        //    w[i][Compressible::RHO_V] = w[i][Compressible::RHO_V]*velocityCoeff;
        //    w[i][Compressible::RHO_W] = w[i][Compressible::RHO_W]*velocityCoeff;
        //    //std::cout << "M i: " << std::endl;
        //}

        if (w[i].internalEnergy() < minInternalEnergy)
        {
            w[i][Compressible::RHO_E] = minInternalEnergy*w[i].density();
            //std::cout << "e min i: " << i << std::endl;
        }

        //if (w[i].internalEnergy() > maxInternalEnergy)
        //{
        //    w[i][Compressible::RHO_E] = maxInternalEnergy*w[i].density();
        //    //std::cout << "e min i: " << i << std::endl;
        //}
    }*/
}

void TwoFluidFVMScheme::updateTimeStep()
{
    const std::vector<Cell>& cells = mesh.getCellList();

    for (int i = 0; i < w.size(); i++)
    {
        Vars<3> projectedArea = cells[i].projectedArea;
        timeSteps[i] = cfl*(cells[i].volume/sum(projectedArea*(abs(u[i].velocityG()) + Vars<3>(u[i].soundSpeedG())))); // pridat vypocet cfl i z liquid zatim pouze z gas phase
        if (timeSteps[i] < 0.0)
        {
            std::cout << "cell index: " << i << std::endl;
        }
    }

    if (!localTimeStep)
    {
        double timeStep = min(timeSteps);
        timeSteps = Field<double>(timeSteps.size(), timeStep);

        time += timeStep;
    }
}

void TwoFluidFVMScheme::calculateFluxes()
{
    fluxes = fluxSolver->calculateFluxes(ul, ur, mesh.getFaceList());
}

Field<Vars<10>> TwoFluidFVMScheme::calculateResidual()
{
    const std::vector<Cell>& cells = mesh.getCellList();
    const std::vector<Face>& faces = mesh.getFaceList();
    const std::vector<int>& owners = mesh.getOwnerIndexList();
    const std::vector<int>& neighbors = mesh.getNeighborIndexList();

    Field<Vars<10>> res(w.size());

    for (int i = 0; i < faces.size(); i++)
    {
        int owner = owners[i];
        int neighbor = neighbors[i];
        
        res[owner] -= fluxes[i];
        if (neighbor >= 0)
        {
            res[neighbor] += fluxes[i];
        }
    }

    for (int i = 0; i < cells.size(); i++)
    {
        res[i] = res[i]/cells[i].volume;
    }
    
    return res;
}

Field<TwoFluid> TwoFluidFVMScheme::getResults() const
{
    return u;
}