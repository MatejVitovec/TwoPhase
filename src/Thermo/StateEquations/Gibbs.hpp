#ifndef GIBBS
#define GIBBS

#include <array>
#include <vector>

#include "NonLinearSolver/NewtonMethod.hpp"

class Gibbs
{
    public:

        Gibbs();

        virtual ~Gibbs() {}

        double rho(double p, double T) const;
        double e(double p, double T) const;
        double s(double p, double T) const;
        double h(double p, double T) const;
        double a2(double p, double T) const;
        double a(double p, double T) const;

        double pFromTRho(double T, double rho, double guessP) const;
        double pFromTE(double T, double e, double guessP) const;
        double pFromTS(double T, double s, double guessP) const;
        double pFromTH(double T, double h, double guessP) const;

        double tFromPRho(double p, double rho, double guessT) const;
        double tFromPE(double p, double e, double guessT) const;
        double tFromPS(double p, double s, double guessT) const;
        double tFromPH(double p, double h, double guessT) const;

        std::pair<double, double> PTFromRhoE(double rho, double e, double guessP, double guessT) const;

        inline double R() const { return specGasConst; }

    protected:
        static constexpr double critT = 647.096;
        static constexpr double critRho = 322.0;
        static constexpr double reducedP_inv = 1000000.0;
        static constexpr double reducedT = 540.0;
        static constexpr double specGasConst = 461.51805;

        NewtonMethod nonLinearSolver;

        inline double tauF(double T) const { return reducedT/T; }
        inline double piF(double p) const { return p*reducedP_inv; }

        virtual double gamma0(double pi, double tau) const = 0;
        virtual double gamma0p(double pi, double tau) const = 0;
        virtual double gamma0p(double pi, double tau) const = 0;
        virtual double gamma0t(double pi, double tau) const = 0;
        virtual double gamma0tt(double pi, double tau) const = 0;
        virtual double gamma0pt(double pi, double tau) const = 0;

        virtual double gammar(double pi, double tau) const = 0;
        virtual double gammarp(double pi, double tau) const = 0;
        virtual double gammarpp(double pi, double tau) const = 0;
        virtual double gammart(double pi, double tau) const = 0;
        virtual double gammartt(double pi, double tau) const = 0;
        virtual double gammarpt(double pi, double tau) const = 0;

        double vFunc(double pi, double tau) const;
        double eFunc(double pi, double tau) const;
        double sFunc(double pi, double tau) const;
        double hFunc(double pi, double tau) const;

        double vDPiFunc(double pi, double tau) const;
        double eDPiFunc(double pi, double tau) const;
        double sDPiFunc(double pi, double tau) const;
        double hDPiFunc(double pi, double tau) const;


        double vDTauFunc(double pi, double tau) const;
        double eDTauFunc(double pi, double tau) const;
        double sDTauFunc(double pi, double tau) const;
        double hDTauFunc(double pi, double tau) const;

};

#endif // GIBBS