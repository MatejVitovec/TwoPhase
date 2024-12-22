#ifndef EXPLICIT_HPP
#define EXPLICIT_HPP

#include "TwoFluidFVMScheme.hpp"

class Explicit : public TwoFluidFVMScheme
{
    public:

        //Explicit() : FVMScheme() {}
        Explicit(Mesh&& mesh_, std::unique_ptr<TwoFluidFlux> fluxSolver_, std::unique_ptr<TwoFluidThermo> thermo_) : TwoFluidFVMScheme(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_)) {}

        virtual ~Explicit() {}

        void solve();        
        

    private:
    

};



#endif // EXPLICIT_HPP