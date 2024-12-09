#include "ExplicitEuler.hpp"
#include <iostream>
#include "outputCFD.hpp"

#include "VolField.hpp"


void ExplicitEuler::solve()
{
    init();

    double lastTmdUpdateTime = 0.0;
    long thermoUpdateFieldTime = 0.0;

    //w = thermo->updateField(w);
    //w.update(thermo->updateField(w));

    thermo->updateThermoInternal(thermoField, w);

    iter = 0;

    bool exitLoop = false;

    //VolField<Compressible> wn = w;
    Vars<5> resNorm;

    auto startAll = std::chrono::high_resolution_clock::now();

    while (iter < maxIter && !exitLoop)
    {
        iter++;

        //boundField();

        updateTimeStep();

        calcBoundaryConditionFields();

        interpolateToFaces();

        calculateFluxes();

        Field<Vars<5>> res = calculateResidual();

        w += res*timeSteps;
        //wn = w;
        //wn += res*timeSteps;

        //auto start = std::chrono::high_resolution_clock::now();

        //wn = thermo->updateField(wn);
        //wn.update(thermo->updateField(wn));

        thermo->updateThermo(thermoField, w);

        //auto stop = std::chrono::high_resolution_clock::now();
        //thermoUpdateFieldTime += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

        if(iter % 100 == 0)
        {
            outputCFD::saveValue(savePath + "/time.txt", std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startAll).count());

            resNorm = res.norm();
            outputCFD::saveResidual(savePath + "/residuals.txt", resNorm);
            std::cout << "iter: " << iter << " density res: " << resNorm[0] << std::endl;
            if(resNorm[0] < targetError) exitLoop = true;

            //outputCFD::saveValue(savePath + "/tmdUpdateTime.txt", thermoUpdateFieldTime/1000.0 - lastTmdUpdateTime);
            //lastTmdUpdateTime = thermoUpdateFieldTime/1000.0;
        }

        if(iter % saveEveryIter == 0)
        {
            outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, thermoField);
        }
    }

    outputCFD::outputVTK(savePath + "/results/results." + std::to_string(iter) + ".vtk", mesh, w, thermoField);
    std::cout << "iter: " << iter << std::endl;

    //std::cout << "time: " << time << std::endl;

    //std::cout << "timeForCalculation: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startAll).count() << std::endl;
}