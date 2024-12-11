#ifndef FVMSCHEME_HPP
#define FVMSCHEME_HPP

#include <vector>
#include <memory>

#include "TwoFluidCompressible.hpp"
#include "TwoFluidPrimitive.hpp"
#include "TwoFluid.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Mesh/Mesh.hpp"
#include "TwoFluidFluxSolver/TwoFluidFluxSolver.hpp"
#include "../Thermo/Thermo.hpp"
#include "../GradientScheme/GradientScheme.hpp"
#include "../Limiter/Limiter.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"

#include "../GradientScheme/LeastSquare.hpp"
#include "../GradientScheme/LeastSquaresPS.hpp"
#include "../Limiter/BarthJespersen.hpp"
#include "../Limiter/Venkatakrishnan.hpp"


class TwoFluidFVMScheme
{
    public:

        TwoFluidFVMScheme(Mesh&& mesh_, std::unique_ptr<TwoFluidFluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : mesh(std::move(mesh_)),
            fluxSolver(std::move(fluxSolver_)),
            thermo(std::move(thermo_)),
            w(Field<TwoFluidCompressible>(mesh.getCellsSize())),
            u(VolField<TwoFluid>(mesh)),
            ul(Field<TwoFluid>()),
            ur(Field<TwoFluid>()),
            cfl(0.8), maxIter(10000000),
            targetError(0000005),
            localTimeStep(false),
            reconstruction(false),
            time(0.0),
            gradientScheme(std::make_unique<GradientScheme>()),
            limiter(std::make_unique<Limiter>()) {}


        virtual ~TwoFluidFVMScheme() {}

        void setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_);
        void setReconstructionLimiter(std::unique_ptr<Limiter> limiter_);

        void setCfl(double cfl_);
        void setMaxIter(int maxIter_);
        void setSaveEveryIter(int saveEveryIter_);
        void setTargetError(double targetError_);
        void setLocalTimeStep(bool localTimeStep_);
        void setReconstructionSettings(bool reconstruction_);
        void setSavePath(std::string path) {savePath = path;}

        void setThermoModel(std::unique_ptr<Thermo> thermo_) {thermo = std::move(thermo_);}
        
        double getCfl() const;
        int getMaxIter() const;
        double getTargetError() const;
        bool getTimeStepsettings() const;
        bool getReconstructionSettings() const;
        std::string getSavePath() const { return savePath; }
        const Mesh& getMesh() const;
        const Thermo* getThermoRef();

        void setInitialConditions(TwoFluidPrimitive initialCondition);

        void setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions);
        void applyFreeBoundaryCondition();

        void init();
        virtual void solve() = 0;

        Field<TwoFluid> getResults() const;
        
        
    protected:
        std::unique_ptr<TwoFluidFluxSolver> fluxSolver;
        std::unique_ptr<Thermo> thermo;
        std::unique_ptr<GradientScheme> gradientScheme;
        std::unique_ptr<Limiter> limiter;

        Mesh mesh;
        std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditionList;

        Field<TwoFluidCompressible> w; //cell size
        VolField<TwoFluid> u;

        Field<TwoFluid> ul; //faces size
        Field<TwoFluid> ur;

        Field<Vars<10>> fluxes;

        Field<Mat<10,3>> grad;
        Field<Vars<10>> phi;

        double cfl;
        int maxIter;
        int saveEveryIter;
        double targetError;
        bool localTimeStep;
        bool reconstruction;
        std::string savePath;

        int fixGradient;
        int fixedGradStep;

        Field<double> timeSteps;

        double time;
        int iter;

        void updateTimeStep();
        void interpolateToFaces();
        void applyBoundaryConditions(); //drive se jmenovalo calcBoundaryConditionFields()
        void calculateFluxes();
        Field<Vars<10>> calculateResidual();
        void boundField();

    private:

};



#endif // FVMSCHEME_HPP