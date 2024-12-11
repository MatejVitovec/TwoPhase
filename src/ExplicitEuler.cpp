#include "ExplicitEuler.hpp"
#include <iostream>
#include "outputCFD.hpp"

#include "VolField.hpp"


void ExplicitEuler::solve()
{
    init();

    thermo->updateThermoInternal(thermoField, w);

    iter = 0;

    bool exitLoop = false;

    Vars<5> resNorm;

    auto startAll = std::chrono::high_resolution_clock::now();

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        calcBoundaryConditionFields();

        //interpolateToFaces();
        interpolateToFacesPrimitive();

        calculateFluxes();

        Field<Vars<5>> res = calculateResidual();

        w += res*timeSteps;

        thermo->updateThermo(thermoField, w);


        if(iter % 100 == 0)
        {
            outputCFD::saveValue(savePath + "/time.txt", std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startAll).count());

            resNorm = res.norm();
            outputCFD::saveResidual(savePath + "/residuals.txt", resNorm);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;
        }

        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, thermoField);
        }
    }

    outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, thermoField);
    std::cout << "iter: " << iter << std::endl;

    //std::cout << "time: " << time << std::endl;
}