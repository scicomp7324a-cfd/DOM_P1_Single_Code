# pragma once

# include <iostream>
# include <vector>
# include "FVMatrix.hpp"
# include "SparseLinearVector.hpp"

using namespace std ;


/* *******************************************************************
    UnderRelaxation.hpp -  This class contains functions that apply 
    velocity and pressure under-relaxation operations used in the 
    SIMPLE iteration procedure.
******************************************************************* */

class UnderRelaxation{

    private:
        double alphaU_ ;
        double alphaP_ ;

    public:
    
    // Constructs the under-relaxation handler with separate relaxation factors for momentum and pressure.
    UnderRelaxation(double alphaU, double alphaP): alphaU_(alphaU),
                                                   alphaP_(alphaP){}

    // Under-relaxes the momentum equation by scaling the diagonal and adding the deferred correction to the source.
    void underRelaxMomentum(FVMatrix& U){
        int numCells = U.sys().size() ;
        for(int i=0; i<numCells; ++i){
            U.addToSource(i, U.sys().X()(i)*((1-alphaU_)/alphaU_)*U.sys().A()(i,i)) ;
            U.sys().A()(i,i) = U.sys().A()(i,i)/alphaU_ ;
        }
    }

    // Updates the pressure field using pressure under-relaxation with the previous pressure solution.
    void underRelaxPressure(FVMatrix& P, SparseLinearVector& pOld){
        int numCells = P.sys().size() ;
        for(int i=0; i<numCells; ++i){
            P.sys().X()(i) = pOld(i) + alphaP_*P.sys().X()(i) ;
        }
    }
} ;