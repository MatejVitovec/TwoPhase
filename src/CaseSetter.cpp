#include <fstream>
#include <iostream>
#include <algorithm> 

#include "CaseSetter.hpp"

//#include "BoundaryCondition/PressureTemperatureInlet.hpp"
//#include "BoundaryCondition/PressureDensityInlet.hpp"
#include "BoundaryCondition/IsentropicInlet.hpp"
#include "BoundaryCondition/PressureOutlet.hpp"
#include "BoundaryCondition/MeanPressureOutlet.hpp"
#include "BoundaryCondition/FreeBoundary.hpp"
#include "BoundaryCondition/Wall.hpp"
#include "BoundaryCondition/Periodicity.hpp"


void CaseSetter::loadSettingFile(std::string fileName)
{
    std::ifstream stream;
    std::string line;

    stream.open(fileName, std::ios_base::in);

    if (stream.is_open())
    {
        while (std::getline(stream, line))
        {
            if(line != "\r")
            {
                data_.push_back(line);
            }            
        }
    }
}


std::unique_ptr<FVMScheme> CaseSetter::createSolver()
{
    std::string solverName = toLower(findParameterByKey("Solver: ", data_));

    Mesh tempMesh = createMesh();
    std::unique_ptr<FluxSolver> tempFluxSolver = createFluxSolver();
    std::unique_ptr<Thermo> tempThermo = createThermoModel();

    if(solverName == "expliciteuler")
    {
        return std::make_unique<ExplicitEuler>(std::move(tempMesh), std::move(tempFluxSolver), std::move(tempThermo));
    }
    else if(solverName == "heunscheme")
    {
        return std::make_unique<HeunScheme>(std::move(tempMesh), std::move(tempFluxSolver), std::move(tempThermo));
    }

    errorMessage("Error: neznamy solver");
    
    std::cout << "nastaven defaultni solver ExplicitEuler" << std::endl;
    return std::make_unique<ExplicitEuler>(std::move(tempMesh), std::move(tempFluxSolver), std::move(tempThermo));
}


std::unique_ptr<FVMScheme> CaseSetter::createSolver(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_)
{
    std::string solverName = toLower(findParameterByKey("Solver: ", data_));

    if(solverName == "expliciteuler")
    {
        return std::make_unique<ExplicitEuler>(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_));
    }
    else if(solverName == "heunscheme")
    {
        return std::make_unique<HeunScheme>(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_));
    }

    errorMessage("Error: neznamy solver");

    std::cout << "nastaven defaultni solver ExplicitEuler" << std::endl;
    return std::make_unique<ExplicitEuler>(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_));
}


std::unique_ptr<FVMScheme> CaseSetter::createAndSetSolver()
{
    std::unique_ptr<FVMScheme> tmpSolver = createSolver(createMesh(), createFluxSolver(), createThermoModel());
    
    if (isParameterDefined("Reconstruction", data_))
    {
        tmpSolver->setReconstructionSettings(true);
        tmpSolver->setReconstructionGradient(createReconstructionGradient());
        tmpSolver->setReconstructionLimiter(createReconstructionLimiter());
    }

    tmpSolver->setCfl(getCfl());
    tmpSolver->setMaxIter(getMaxIter());
    tmpSolver->setTargetError(getTargetError());
    tmpSolver->setLocalTimeStep(getLocalTimeStepSetting());
    tmpSolver->setSaveEveryIter(getSaveIterInterval());

    tmpSolver->setBoundaryConditions(createBoundaryCondition(tmpSolver->getMesh(), tmpSolver->getThermoRef()));
    tmpSolver->setInitialConditions(getInitialCondition(tmpSolver->getThermoRef()));

    return std::move(tmpSolver);
}


Mesh CaseSetter::createMesh()
{
    std::vector<std::string> data = findParametersByKey("Mesh", data_);

    std::string type = toLower(findParameterByKey("type: ", data));

    std::string meshPath = findParameterByKey("fileName: ", data);

    Mesh tempMesh = Mesh();

    if(type == "gmsh2")
    {
        tempMesh.loadGmsh2(meshPath);

        return tempMesh;
    }
    
    std::cout << "Error: neznamy format meshe" << std::endl;
    return tempMesh;
}


std::unique_ptr<FluxSolver> CaseSetter::createFluxSolver()
{
    std::string fluxSolverName = toLower(findParameterByKey("RiemannSolver: ", data_));

    if(fluxSolverName == "hllc")
    {
        return std::make_unique<Hllc>();
    }
    else if(fluxSolverName == "hll")
    {
        return std::make_unique<Hll>();
    }
    else if(fluxSolverName == "ausm+")
    {
        return std::make_unique<Hllc>(); //TODO
    }

    errorMessage("neznamy Rieman solver");

    std::cout << "Nastaven defaultni Rieman solver HLLC" << std::endl;
    return std::make_unique<Hllc>();
}


std::unique_ptr<Thermo> CaseSetter::createThermoModel()
{
    std::vector<std::string> data = findParametersByKey("ThermoModel", data_);

    std::string name = toLower(findParameterByKey("name: ", data));
    
    if(name == "idealgas")
    {
        std::string gamma = findParameterByKey("gamma: ", data);
        std::string specGasConst = findParameterByKey("specificgasconstant: ", data);

        if(gamma == "" || specGasConst == "")
        {
            return std::make_unique<IdealGasThermo>();
        }
        else
        {
            return std::make_unique<IdealGasThermo>(std::stod(gamma), std::stod(specGasConst));
        }
    }
    else if(name == "iapws95")
    {
        std::string interpolationName = findParameterByKey("interpolation: ", data);

        if(interpolationName == "bilinear")
        {
            std::cout << "Nacitani biLin interpolace" << std::endl;
            return std::make_unique<Iapws95InterpolationThermo<BiLinearInterpolation>>();
        }
        else if(interpolationName == "biquadratic")
        {
            std::cout << "Nacitani biQuad interpolace" << std::endl;
            return std::make_unique<Iapws95InterpolationThermo<BiQuadraticInterpolation>>();
        }

        return std::make_unique<Iapws95Thermo>();
    }
    else if(name == "specialgasequation")
    {
        std::cout << "Nacitani specialGasEquation modelu" << std::endl;
        return std::make_unique<SpecialGasThermo>();
    }
    /*else if(name == "iapws95interpolation")
    {
        return std::make_unique<Iapws95InterpolationThermo<BiQuadraticInterpolation>>();
    }*/

    errorMessage("neznamy termo solver");

    std::cout << "Nastaven defaultni thermo model idealGas(gamma = 1.4, R = 287.05)" << std::endl;
    return std::make_unique<IdealGasThermo>();
}

std::unique_ptr<GradientScheme> CaseSetter::createReconstructionGradient()
{
    std::vector<std::string> data = findParametersByKey("Reconstruction", data_);

    std::string name = toLower(findParameterByKey("gradientSolver: ", data));

    if (name == "leastsquares")
    {
        return std::make_unique<LeastSquare>();
    }
    else if (name == "leastsquaresps")
    {
        return std::make_unique<LeastSquaresPS>();
    }
    
    
    std::cout << "Gradient solver not defined" << std::endl;
    return std::make_unique<GradientScheme>();
}

std::unique_ptr<Limiter> CaseSetter::createReconstructionLimiter()
{
    std::vector<std::string> data = findParametersByKey("Reconstruction", data_);

    std::string name = toLower(findParameterByKey("limiter: ", data));

    if (name == "barthjespersen")
    {
        return std::make_unique<BarthJespersen>();
    }
    else if(name == "venkatakrishnan")
    {
        return std::make_unique<Venkatakrishnan>();
    }
    else if(name == "cubic")
    {
        return std::make_unique<CubicLimiter>(2);
    }


    std::cout << "Gradient limiter not defined" << std::endl;
    return std::make_unique<Limiter>();
}


bool CaseSetter::getLocalTimeStepSetting()
{
    return stringToBool(findParameterByKey("LocalTimeStep:", data_));
}

double CaseSetter::getTargetError()
{
    return std::stod(findParameterByKey("TargetError: ", data_));
}

double CaseSetter::getMaxIter()
{
    return std::stod(findParameterByKey("MaxIter:", data_));
}

double CaseSetter::getCfl()
{
    return std::stod(findParameterByKey("cfl:", data_));
}

int CaseSetter::getSaveIterInterval()
{
    return std::stoi(findParameterByKey("SaveEveryIter:", data_));
}


Compressible CaseSetter::getInitialCondition(const Thermo * const thermoModel)
{
    std::vector<std::string> data = findParametersByKey("InitialCondition", data_);

    std::string name = toLower(findParameterByKey("type: ", data));
    std::vector<double> value = vectorStringToDouble(stringArrayToVectorOfStrings(findParameterByKey("val: ", data)));

    if(value.size() != 5)
    {
        std::cout << "Error: spatne zadane initialCondition" << std::endl;
    }


    if(name == "primitive")
    {
        //std::unique_ptr<Thermo> tempThermo = createThermoModel();
        //return tempThermo->primitiveToConservative(Vars<5>({value[0], value[1], value[2], value[3], value[4]}));
        return thermoModel->primitiveToConservative(Vars<5>({value[0], value[1], value[2], value[3], value[4]}));
    }
    else if(name == "conservative")
    {
        return Compressible({value[0], value[1], value[2], value[3], value[4]});
    }
    else if(name == "stagnationstate")
    {
        //std::unique_ptr<Thermo> tempThermo = createThermoModel();
        //return tempThermo->stagnationState(value[0], value[4]);
        return thermoModel->stagnationState(value[0], value[4]);
    }
    else
    {
        std::cout << "Error: neznama pocatecni podminka" << std::endl;
    }

    return Compressible();
}


std::vector<std::shared_ptr<BoundaryCondition>> CaseSetter::createBoundaryCondition(const Mesh& mesh, const Thermo * const thermoModel)
{
    const std::vector<Boundary>& meshBoundaryList = mesh.getBoundaryList();

    std::vector<std::shared_ptr<BoundaryCondition>> out(meshBoundaryList.size());

    std::vector<std::string> boundarydata = findParametersByKey("BoundaryCondition", data_);
    std::vector<std::string> boundaryNames = findObjectNamesInGroup(boundarydata);

    for (size_t i = 0; i < boundaryNames.size(); i++)
    {
        std::vector<std::string> data = findParametersByKey(boundaryNames[i], boundarydata);
        std::string type = toLower(findParameterByKey("type: ", data));

        Boundary auxBoundary;
		bool exist = false;
        int boundaryId;
        for (boundaryId = 0; boundaryId < meshBoundaryList.size(); boundaryId++)
        {
			if(meshBoundaryList[boundaryId].boundaryConditionName == boundaryNames[i])
			{
				auxBoundary = meshBoundaryList[boundaryId];
				exist = true;
				break;
			}
		}
		/*for (const auto & meshBoundary : meshBoundaryList)
		{
			if(meshBoundary.boundaryConditionName == boundaryNames[i])
			{
				auxBoundary = meshBoundary;
				exist = true;
				break;
			}
		}*/
		if(!exist)
		{
			std::cout << "Chyba BC" << std::endl;
			continue;
		}

        if (type == "pressuretemperatureinlet")
        {
            std::string totalPressure = findParameterByKey("totalPressure: ", data);
            std::string totalTemperature = findParameterByKey("totalTemperature: ", data);
            std::string xyAngle = findParameterByKey("xyAngle: ", data);
            std::string xzAngle = findParameterByKey("xzAngle: ", data);

            if(totalPressure != "" && totalTemperature != "" && xyAngle != "" && xzAngle != "")
            {
                double pTot = std::stod(totalPressure);
                double TTot = std::stod(totalTemperature);

                std::array<double, 3> state = thermoModel->initPressureTemperatureInlet(pTot, TTot);

                out[boundaryId] = std::make_shared<IsentropicInlet>(auxBoundary,
                                                                state[0],
                                                                pTot,
                                                                TTot,
                                                                state[1],
                                                                state[2],
                                                                angleAngleToUnit(std::stod(xyAngle), std::stod(xzAngle)));
            }
            else { errorMessage("spatne zadane parametry BC"); }
        }
        else if (type == "pressureoutlet")
        {
            std::string pressure = findParameterByKey("pressure: ", data);

            if(pressure != "")
            {
                out[boundaryId] = std::make_shared<PressureOutlet>(auxBoundary, std::stod(pressure));
            }
            else { errorMessage("spatne zadane parametry BC"); }
        }
        else if (type == "meanpressureoutlet")
        {
            std::string pressure = findParameterByKey("pressure: ", data);

            if(pressure != "")
            {
                out[boundaryId] = std::make_shared<MeanPressureOutlet>(auxBoundary, std::stod(pressure));
            }
            else { errorMessage("spatne zadane parametry BC"); }
        }
        else if (type == "wall")
        {
            out[boundaryId] = std::make_shared<Wall>(auxBoundary);
        }
        else if (type == "periodicity")
        {
            std::vector<double> shiftArray = vectorStringToDouble(stringArrayToVectorOfStrings(findParameterByKey("shift: ", data)));

            if(shiftArray.size() == 3)
            {
                out[boundaryId] = std::make_shared<Periodicity>(auxBoundary, Vars<3>({shiftArray[0], shiftArray[1], shiftArray[2]}), "000", mesh);
            }
            else { errorMessage("spatne zadane parametry BC"); }
        }
        else if (type == "symmetry")
        {
            out[boundaryId] = std::make_shared<Wall>(auxBoundary);
        }
        else if (type == "free")
        {
            out[boundaryId] = std::make_shared<FreeBoundary>(auxBoundary);
        }
        else
        {
            errorMessage("chyba v zadani okrajove podmink");
            out[boundaryId] = std::make_shared<FreeBoundary>(auxBoundary);
        }        
    }

    return out;
}

//PRIVATE

void CaseSetter::errorMessage(std::string str)
{
    std::cout << "Error: " << str << std::endl;
}

std::string CaseSetter::toLower(std::string str)
{
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    
    return str;
}

std::string CaseSetter::removeWhiteSpaces(std::string str)
{
    str.erase(std::remove_if(str.begin(), str.end(),
            [](char c) {
                return (c == ' ' || c == '\n' || c == '\r' ||
                    c == '\t' || c == '\v' || c == '\f');
            }),
        str.end());

    return str;
}

bool CaseSetter::stringToBool(std::string str)
{
    transform(str.begin(), str.end(), str.begin(), ::tolower);

    if(str == "true")
        return true;
    else if(str == "false")
        return false;
    else
        std::cout << "Eror konverze stringu na bool" << std::endl;

    return false;
}

std::vector<std::string> CaseSetter::stringArrayToVectorOfStrings(std::string str)
{
    std::vector<std::string> out;
    out.clear();

    std::size_t start = 0;
    std::size_t end = 0;

    std::size_t found = str.find("[");
    if (found != std::string::npos)
        start = found;
    else
        std::cout << "Error: spatny zapis pole hodnot" << std::endl;

    found = str.find("]");
    if (found != std::string::npos)
        end = found;
    else
        std::cout << "Error: spatny zapis pole hodnot" << std::endl;


    str = str.substr(start + 1, end - start - 1);

    bool stop = false;
    while(!stop)
    {
        found = str.find(",");
        if (found != std::string::npos)
        {
            out.push_back(str.substr(0, found));
            str = str.substr(found + 1);
        }
        else
        {
            stop = true;
            out.push_back(str);
        }
    }

    return out;
}

bool CaseSetter::isParameterDefined(std::string key, std::vector<std::string> data)
{
    int level = 0;    

    for (size_t i = 0; i < data.size(); i++)
    {
        if (data[i].find('{') != std::string::npos)
        {
            level++;
        }
        if (data[i].find('}') != std::string::npos)
        {
            level--;
        }

        if(level == 0)
        {
            std::size_t found = data[i].find(key);
            if (found != std::string::npos)
            {
                return true;
            }
        }
    }

    return false;
}

std::string CaseSetter::findParameterByKey(std::string key, std::vector<std::string> data)
{
    int level = 0;    

    for (size_t i = 0; i < data.size(); i++)
    {
        if (data[i].find('{') != std::string::npos)
        {
            level++;
        }
        if (data[i].find('}') != std::string::npos)
        {
            level--;
        }

        if(level == 0)
        {
            std::size_t found = data[i].find(key);
            if (found != std::string::npos)
            {
                return removeWhiteSpaces(data[i].substr(found + key.size()));
            }
        }
    }

    return "";
}

std::vector<double> CaseSetter::vectorStringToDouble(std::vector<std::string> in)
{
    std::vector<double> out;
    out.clear();

    for (size_t i = 0; i < in.size(); i++)
    {
        out.push_back(std::stod(in[i]));
    }
    
    return out;
}

std::vector<std::string> CaseSetter::findParametersByKey(std::string key, std::vector<std::string> data)
{
    std::vector<std::string> out;
    out.clear();

    int level = 0;
    int startIndex = 0;

    if (data[0].find('{') != std::string::npos)
    {
        level++;
    }

    for (size_t i = 0; i < data.size(); i++)
    {
        if(level == 0)
        {
            std::size_t found = data[i].find(key);
            if (found != std::string::npos)
            {
                if((data[i].substr(key.size())).find('{') != std::string::npos)
                {
                    startIndex = i + 1;
                    break;
                }
                if(data[i+1].find('{') != std::string::npos)
                {
                    startIndex = i + 2;
                    break;
                }
            }
        }

        if (data[i].find('{') != std::string::npos)
        {
            level++;
        }
        if (data[i].find('}') != std::string::npos)
        {
            level--;
        }
    }

    level = 0;
    for (size_t i = startIndex; i < data.size(); i++)
    {
        if (data[i].find('{') != std::string::npos)
        {
            level++;
        }

        if (level == 0)
        {
            if (data[i].find('}') != std::string::npos)
            {
                break;
            }
        }
        else
        {
            if (data[i].find('}') != std::string::npos)
            {
                level--;
            }
        }

        out.push_back(data[i]);
    }
    
    return out;
}

std::vector<std::string> CaseSetter::findObjectNamesInGroup(std::vector<std::string> data)
{
    std::vector<std::string> out;
    out.clear();

    int level = 0;

    if (data[0].find('{') != std::string::npos)
    {
        level++;
    }

    for (size_t i = 0; i < data.size() - 1; i++)
    {
        std::string rwsData = removeWhiteSpaces(data[i]);
        if(level == 0)
        {
            if(rwsData != "" && data[i].find('{') == std::string::npos && removeWhiteSpaces(data[i+1])[0] == '{')
            {
                out.push_back(rwsData);
            }
            else if(rwsData != "" && rwsData[rwsData.size()-1] == '{' && rwsData.size() > 1)
            {
                out.push_back(rwsData.substr(0, rwsData.size()-1));
            }
        }

        if (data[i].find('{') != std::string::npos)
        {
            level++;
        }
        if (data[i].find('}') != std::string::npos)
        {
            level--;
        }
    }

    return out;
}