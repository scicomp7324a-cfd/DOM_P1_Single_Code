# pragma once

# include<iostream>
# include<vector>
# include "SparseLinearSystem.hpp"

using namespace std ;

/* *******************************************************************
    GaussSeidelSolver.hpp - This class contains the smoother 
    implementation.
******************************************************************* */

class GaussSeidelSmoother{

    private:
    SparseLinearSystem& sys_ ;
    
    public:
    // Constructs the gauss seidel smoother object using the sparse linear system 
    explicit GaussSeidelSmoother(SparseLinearSystem& sys) : sys_(sys) {}

    // Function containing smoother implementation
    void smooth(int nSweeps){
        
       const vector<int> rowArr_ = sys_.A().addressing().rowArray() ;
       const vector<int> colIndArr_ = sys_.A().addressing().columnIndexArray() ; 
       int n = 0, numEq = sys_.A().dim();
       double sum = 0.0;

       while(n!=nSweeps){

        for(int i=0; i<numEq; ++i){
            for(int j = rowArr_[i]; j < rowArr_[i+1]; ++j){
                sum += sys_.A()(i,colIndArr_[j])*sys_.X()(colIndArr_[j]) ;
            }
            sys_.X()(i) = (sys_.b()(i) - sum)/(sys_.A()(i,i)) ;
            sum = 0.0 ;
        }
        
        ++n ;
       }

    }

} ;