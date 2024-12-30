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

#include "BoundaryCondition/Wall.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"

Field<TwoFluidPrimitive> riemannInitialCond(const TwoFluidPrimitive& leftState, const TwoFluidPrimitive& rightState, const int N)
{
    Field<TwoFluidPrimitive> out = Field<TwoFluidPrimitive>(N);

    int i = 0;
    for ( ; i < N/2; i++)
    {
        out[i] = leftState;
    }
    for ( ; i < N; i++)
    {
        out[i] = rightState;
    }
    
    return out;
}


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    Mesh tempMesh = Mesh();
    tempMesh.loadGmsh2("../meshes/riemannMesh.msh");

    std::unique_ptr<TwoFluidFVMScheme> solver = std::make_unique<Explicit>(std::move(tempMesh), std::make_unique<TwoFluidSlau2>(), std::make_unique<TwoFluidThermo>());

    solver->setSavePath("../results/RiemannTest1");


    solver->setCfl(0.2);

    std::vector<std::shared_ptr<BoundaryCondition>> bcList(solver->getMesh().getBoundarySize());
    bcList[0] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[0], 0);
    bcList[1] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[1], 1);
    bcList[2] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[2], 2);
    bcList[3] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[3], 3);
    bcList[4] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[4], 4);
    bcList[5] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[5], 5);

    
    solver->setBoundaryConditions(bcList);

    const double epsilon = 1.0e-7;

                                                                                   //alpha,      p,   tg,    tl,   ug,  ul,  vg,  vl,  wg,  wl
    Field<TwoFluidPrimitive> initialCond = riemannInitialCond(TwoFluidPrimitive({1.0 - epsilon, 1e5, 300.0, 300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
                                                              TwoFluidPrimitive({epsilon,       1e5, 300.0, 300.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}), 200);

    solver->setInitialConditions(initialCond);

    //outputCFD::outputVTK("../results/RiemannTest1/results/results.0.vtk", solver->getMesh(), solver->getResults());

    solver->solve();

    outputCFD::outputVTK("../results/RiemannTest1/results/resultsFinal.vtk", solver->getMesh(), solver->getResults());

    return 0;
}