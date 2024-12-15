#include <cmath>
#include <iostream>
#include <fstream>

#include "Gibbs.hpp"

Gibbs::Gibbs()
{
    nonLinearSolver = NewtonMethod();
}

double Gibbs::rho(double p, double T) const
{
    return 1.0/vFunc(piF(p), tauF(T));
}

double Gibbs::e(double p, double T) const
{
    return eFunc(piF(p), tauF(T));
}

double Gibbs::s(double p, double T) const
{
    return sFunc(piF(p), tauF(T));
}

double Gibbs::h(double p, double T) const
{
    return hFunc(piF(p), tauF(T));
}

double Gibbs::a2(double p, double T) const
{
    double pi = piF(p);
    double tau = tauF(T);

    double gammarpAux = gammarp(pi, tau);

    return specGasConst*T*(1.0 + 2*pi*gammarpAux + std::pow(pi*gammarpAux, 2.0))/(1.0 - std::pow(pi, 2.0)*gammarpp(pi, tau) + std::pow(1.0 + pi*gammarpAux - tau*pi*gammarpt(pi, tau), 2.0)/(std::pow(tau, 2.0)*(gamma0tt(pi, tau) + gammartt(pi, tau))));
}

double Gibbs::a(double p, double T) const
{
    return std::sqrt(a2(p, T));
}

double Gibbs::pFromTRho(double T, double rho, double guessP) const
{
    double tau = tauF(T);
    double v = 1.0/rho;

    double pi = nonLinearSolver.solve([=](double val) { return vFunc(val, tau) - v; },
                                      [=](double val) { return vDPiFunc(val, tau);},
                                      piF(guessP));
    
    return pi*reducedP_inv;
}

double Gibbs::pFromTE(double T, double e, double guessP) const
{
    double tau = tauF(T);

    double pi = nonLinearSolver.solve([=](double val) { return eFunc(val, tau) - e; },
                                      [=](double val) { return eDPiFunc(val, tau);},
                                      piF(guessP));
    
    return pi*reducedP_inv;
}

double Gibbs::pFromTS(double T, double s, double guessP) const
{
    double tau = tauF(T);

    double pi = nonLinearSolver.solve([=](double val) { return sFunc(val, tau) - s; },
                                      [=](double val) { return sDPiFunc(val, tau);},
                                      piF(guessP));
    
    return pi*reducedP_inv;
}

double Gibbs::pFromTH(double T, double h, double guessP) const
{
    double tau = tauF(T);

    double pi = nonLinearSolver.solve([=](double val) { return hFunc(val, tau) - h; },
                                      [=](double val) { return hDPiFunc(val, tau);},
                                      piF(guessP));
    
    return pi*reducedP_inv;
}


double Gibbs::tFromPRho(double p, double rho, double guessT) const
{
    double pi = piF(p);
    double v = 1.0/rho;

    double tau = nonLinearSolver.solve([=](double val) { return vFunc(pi, val) - v; },
                                      [=](double val) { return vDTauFunc(pi, val);},
                                      tauF(guessT));
    
    return reducedT/tau;
}

double Gibbs::tFromPE(double p, double e, double guessT) const
{
    double pi = piF(p);

    double tau = nonLinearSolver.solve([=](double val) { return eFunc(pi, val) - e; },
                                      [=](double val) { return eDTauFunc(pi, val);},
                                      tauF(guessT));
    
    return reducedT/tau;
}

double Gibbs::tFromPS(double p, double s, double guessT) const
{
    double pi = piF(p);

    double tau = nonLinearSolver.solve([=](double val) { return sFunc(pi, val) - s; },
                                      [=](double val) { return sDTauFunc(pi, val);},
                                      tauF(guessT));
    
    return reducedT/tau;
}

double Gibbs::tFromPH(double p, double h, double guessT) const
{
    double pi = piF(p);

    double tau = nonLinearSolver.solve([=](double val) { return hFunc(pi, val) - h; },
                                      [=](double val) { return hDTauFunc(pi, val);},
                                      tauF(guessT));
    
    return reducedT/tau;
}


std::pair<double, double> Gibbs::PTFromRhoE(double rho, double e, double guessP, double guessT) const
{
    double v = 1.0/rho;
    std::pair<double, double> result = nonLinearSolver.solve([=](double val1, double val2) { return vFunc(val1, val2) - v; },
                                                             [=](double val1, double val2) { return vDPiFunc(val1, val2); },
                                                             [=](double val1, double val2) { return vDTauFunc(val1, val2); },
                                                             [=](double val1, double val2) { return eFunc(val1, val2) - e; },
                                                             [=](double val1, double val2) { return eDPiFunc(val1, val2); },
                                                             [=](double val1, double val2) { return eDTauFunc(val1, val2); },
                                                             piF(guessP),
                                                             tauF(guessT));

    return std::make_pair(piF(result.first), (result.second));
}


double Gibbs::vFunc(double pi, double tau) const
{
    return ((specGasConst*reducedT*reducedP_inv)/tau)*(gamma0p(pi, tau) + gammarp(pi, tau));
}

double Gibbs::eFunc(double pi, double tau) const
{
    return specGasConst*(reducedT/tau)*(tau*(gamma0t(pi, tau) + gammart(pi, tau)) - pi*(gamma0p(pi, tau) + gammarp(pi, tau)));
}

double Gibbs::sFunc(double pi, double tau) const
{
    return specGasConst*(tau*(gamma0t(pi, tau) + gammart(pi, tau)) - (gamma0(pi, tau) + gammar(pi, tau)));
}

double Gibbs::hFunc(double pi, double tau) const
{
    specGasConst*reducedT*(gamma0t(pi, tau) + gammart(pi, tau));
}

double Gibbs::vDPiFunc(double pi, double tau) const
{
    
}

double Gibbs::eDPiFunc(double pi, double tau) const
{
    
}

double Gibbs::sDPiFunc(double pi, double tau) const
{
    
}

double Gibbs::hDPiFunc(double pi, double tau) const
{
    
}

double Gibbs::vDTauFunc(double pi, double tau) const
{
    
}

double Gibbs::eDTauFunc(double pi, double tau) const
{
    
}

double Gibbs::sDTauFunc(double pi, double tau) const
{
    
}

double Gibbs::hDTauFunc(double pi, double tau) const
{
    
}