#include "Explicit.hpp"
#include <iostream>
#include "../outputCFD.hpp"

#include "../VolField.hpp"


void Explicit::solve()
{
    init();

    outputCFD::outputVTK(savePath + "/results/results." + std::to_string(0) + ".vtk", mesh, u);

    iter = 0;

    bool exitLoop = false;

    Vars<10> resNorm;

    auto startAll = std::chrono::high_resolution_clock::now();

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        applyBoundaryConditions();

        pInt = calculateInterfacialPressure();
        updateConservative();

        interpolateToFaces();
        fluxSolver->calculateFluxes(ul, ur, mesh.getFaceList(), fluxesl, fluxesr);

        Field<Vars<10>> res = calculateResidual();

        w += res*timeSteps;

        thermo->updateFromConservative(u, w, pInt);
        blend();        
        thermo->updateInternal(u);
        //blend();

        if(iter % 100 == 0)
        {
            outputCFD::saveValue(savePath + "/time.txt", std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startAll).count());

            resNorm = res.norm();
            outputCFD::saveResidual(savePath + "/residuals.txt", resNorm[0]);
            std::cout << "iter: " << iter << ", density res: " << resNorm[0] << ", time: " << time << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;
        }

        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, u);
        }

        if(time > 0.002)
        {
            exitLoop = true;
        }
    }

    outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, u);
    std::cout << "iter: " << iter << std::endl;

    //std::cout << "time: " << time << std::endl;
}