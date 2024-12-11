#ifndef HLLC_HPP
#define HLLC_HPP

#include "FluxSolver.hpp"

class Hllc : public FluxSolver
{
    public:

        Hllc() {}

        virtual ~Hllc() {}

        Vars<5> claculateFlux(const double rhoL, const double rhoR,
                              const Vars<3> uL, const Vars<3> uR,
                              const double pL, const double pR,
                              const double eL, const double eR,
                              const double aL, const double aR,
                              const Vars<3>& normalVector) const;

        Vars<5> claculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const;
        Vars<9> claculateFlux(const CompressibleMixture& wl, const CompressibleMixture& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const;

    private:
        Vars<3> waveSpeedsEstimate2(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const;

};

#endif // HLLC_HPP