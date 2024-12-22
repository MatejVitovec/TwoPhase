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
#include "TwoFluid/TwoFluidFlux/TwoFluidFlux.hpp"

#include "BoundaryCondition/Wall.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"



int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    Mesh tempMesh = Mesh();
    tempMesh.loadGmsh2("../results/TwoFluid/test1.msh");


    std::unique_ptr<TwoFluidFVMScheme> solver = std::make_unique<Explicit>(std::move(tempMesh), std::make_unique<TwoFluidFlux>(), std::make_unique<StiffenedGasThermo>());

    solver->setCfl(0.2);

    std::vector<std::shared_ptr<BoundaryCondition>> bcList(solver->getMesh().getBoundarySize());
    bcList[0] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[0], 0);
    bcList[1] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[1], 1);
    bcList[2] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[2], 2);
    bcList[3] = std::make_shared<FreeBoundary>(solver->getMesh().getBoundaryList()[3], 3);
    bcList[4] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[4], 4);
    bcList[5] = std::make_shared<Wall>(solver->getMesh().getBoundaryList()[5], 5);


    solver->setBoundaryConditions(bcList);
    solver->setInitialConditions(TwoFluidPrimitive({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));

    solver->solve();

    return 0;
}