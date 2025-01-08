#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>

#include "outputCFD.hpp"
#include "CaseSetter.hpp"

#include "TwoFluid/TwoFluidFVMScheme.hpp"
#include "TwoFluid/Explicit.hpp"
#include "TwoFluid/TwoFluidFlux/TwoFluidSlau2.hpp"

#include "BoundaryCondition/Inlet.hpp"
#include "BoundaryCondition/PressureOutlet.hpp"
#include "BoundaryCondition/Wall.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"

Field<TwoFluidPrimitive> riemannInitialCond(const TwoFluidPrimitive& leftState, const TwoFluidPrimitive& rightState, double x, const Mesh& mesh)
{
    int N = mesh.getCellsSize();
    const std::vector<Cell>& cellList = mesh.getCellList();

    Field<TwoFluidPrimitive> out = Field<TwoFluidPrimitive>(N);

    for (int i = 0; i < N; i++)
    {
        if (cellList[i].center[0] <= x)
        {
            out[i] = leftState;
        }
        else
        {
            out[i] = rightState;
        }
    }
    
    return out;
}

void dropletAlphaInitialCond(double epsilon, double radius, const double delta, const Mesh& mesh, Field<TwoFluidPrimitive>& u)
{
    const std::vector<Cell>& cellList = mesh.getCellList();

    for (int i = 0; i < u.size(); i++)
    {
        double cellRadius = norm2(cellList[i].center);
        if (radius - 2.0*delta <= cellRadius && cellRadius <= radius + 2.0*delta)
        {
            double xi = (cellRadius - (radius - 2.0*delta))/(4.0*delta);
            double G = -std::pow(xi, 2.0)*(2.0*xi - 3.0);
            u[i][TwoFluid::ALPHA] = G*epsilon + (1.0 - G)*(1.0 - epsilon);
        }
        else if(cellRadius < radius)
        {
            u[i][TwoFluid::ALPHA] = 1.0 - epsilon;
        }
    }
}



int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    Mesh tempMesh = Mesh();
    tempMesh.loadGmsh2("../meshes/dropletMeshSmall.msh");

    std::unique_ptr<TwoFluidFVMScheme> solver = std::make_unique<Explicit>(std::move(tempMesh), std::make_unique<TwoFluidSlau2>(), std::make_unique<TwoFluidThermo>());

    solver->setSavePath("../results/Droplet");

    solver->setReconstructionGradient(std::make_unique<LeastSquaresPS>());
    solver->setReconstructionLimiter(std::make_unique<Venkatakrishnan>());
    solver->setReconstructionSettings(true);

    solver->setCfl(0.1);
    solver->setMaxIter(24000);
    solver->setSaveEveryIter(100);

    std::vector<std::shared_ptr<BoundaryCondition>> bcList(solver->getMesh().getBoundarySize());
    bcList[0] = std::make_shared<Inlet>(solver->getMesh().getBoundaryList()[0], 381.85, Vars<3>({225.86, 0.0, 0.0}), 0);
    bcList[1] = std::make_shared<PressureOutlet>(solver->getMesh().getBoundaryList()[1], 1.0e5, 1);
    bcList[2] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[2], 2);
    bcList[3] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[3], 3);
    //bcList[4] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[4], 4);
    //bcList[5] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[5], 5);

    
    solver->setBoundaryConditions(bcList);

    const double epsilon = 1.0e-7;

                                                                                   //alpha,       p,     tg,     tl,     ug,  vg,  wg,     ul,  vl,  wl
    Field<TwoFluidPrimitive> initialCond = riemannInitialCond(TwoFluidPrimitive({epsilon, 2.35438e5, 381.85, 381.85, 225.86, 0.0, 0.0, 225.86, 0.0, 0.0}),
                                                              TwoFluidPrimitive({epsilon,     1.0e5, 293.15, 293.15,    0.0, 0.0, 0.0,    0.0, 0.0, 0.0}),
                                                              -0.004, solver->getMesh());

    dropletAlphaInitialCond(epsilon, 0.0032, 0.000025, solver->getMesh(), initialCond);

    //Field<TwoFluidPrimitive> initialCond = riemannInitialCond(TwoFluidPrimitive({epsilon,       1e7, 308.15, 308.15, 100.0, 0.0, 0.0, 100.0, 0.0, 0.0}),
    //                                                          TwoFluidPrimitive({1.0 - epsilon, 1e6, 308.15, 308.15, 100.0, 0.0, 0.0, 100.0, 0.0, 0.0}), solver->getMesh().getCellsSize());

    //Field<TwoFluidPrimitive> initialCond = riemannInitialCond(TwoFluidPrimitive({1.0 - epsilon, 1e5, 300.0, 300.0, 100.0, 0.0, 0.0, 100.0, 0.0, 0.0}),
    //                                                          TwoFluidPrimitive({epsilon,       1e5, 300.0, 300.0, 100.0, 0.0, 0.0, 100.0, 0.0, 0.0}), 200);

    solver->setInitialConditions(initialCond);

    outputCFD::outputVTK("../results/Droplet/results/results.0.vtk", solver->getMesh(), solver->getResults());

    solver->solve();

    outputCFD::outputVTK("../results/RiemannTest1/results/resultsFinal.vtk", solver->getMesh(), solver->getResults());

    return 0;
}