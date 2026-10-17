# pragma once 

# include <iostream> 
# include "SparseLinearSystem.hpp"
# include "GaussSeidelSolver.hpp"
# include "SparseAddress.hpp"
# include "SparseLinearVector.hpp"


using namespace std ;

/* *******************************************************************
    FVMatrix.hpp - The FVMatrix is the linear system manager which 
    controls the updates and initializations to the discretization
    matrix, source vector and solution vector.   
******************************************************************* */

class FVMatrix{

    private:
    SparseLinearSystem sys_ ;
    unique_ptr<GaussSeidelSolver> solver_ ;

    public:

    // Constructs the Fvmatrix by constructing the linear system. The FV Matrix is directly constructed without the need to independently construct a linear system
    FVMatrix(const SparseAddress& addr): sys_(addr)
    {
        solver_ = make_unique<GaussSeidelSolver>(sys_) ;
    }

    // non-constant acces to sparse linear system
    SparseLinearSystem& sys(){
        return sys_ ;
    }

    const SparseLinearSystem& sys() const{
         return sys_ ;
    }

    // Function to clear the discretization matrix and source vector. Required to construct new matrix coefficients at the start of each loop
    void clear(){

        const vector<int> rowArr = sys_.A().addressing().rowArray() ;
        const vector<int> colIndArr = sys_.A().addressing().columnIndexArray() ;

        for(int i = 0; i < sys_.A().dim(); ++i){
            sys_.A()(i,i) = 0 ;
            for(int j = rowArr[i]; j < rowArr[i+1]; ++j){
                sys_.A()(i, colIndArr[j]) = 0 ;
            }
        }

        for(int i = 0; i < sys_.A().dim(); ++i){
            sys_.b()(i) = 0 ; 
        }

    }

    // Function to initialize the solution vector 
    void initializeSolution(SparseLinearVector& uInit){
        for(int i = 0; i < sys_.A().dim(); ++i){
            sys_.X()(i) = uInit(i) ; 
        }
    }

    // Funtion to add to the diagonal coefficient of the sparse matrix
    void addDiagonal(int cellIndex, double val){
        sys_.A()(cellIndex, cellIndex) += val ;
    }

    // Function to add to the non-diagonal coefficient of the sparse matrix
    void addCoefficient(int ownerIndex, int neighbourIndex, double val){
        sys_.A()(ownerIndex, neighbourIndex) += val ;
    }

    // Function to add to the source vector element of the matrix
    void addToSource(int cellIndex, double val){
        sys_.b()(cellIndex) += val ;
    }

    // Function to set the diagonal coefficient value of the matrix
    void setDiagonal(int cellIndex, double val){
        sys_.A()(cellIndex, cellIndex) = val ;
    }

    // Funciton to set the non-diagonal coefficient of the matrix 
    void setCoefficient(int ownerIndex, int neighbourIndex, double val){
        sys_.A()(ownerIndex, neighbourIndex) = val ;
    }

    // Function to set the source vector element of the matrix 
    void setSource(int cellIndex, double val){
        sys_.b()(cellIndex) = val ;
    }

    // Function call to gauss-seidel solver to solve the sparse linear system
    void solve(double relTol, bool verbose, ostream& outInner, ostream& outFinal, int loop_count){
        solver_->solve(relTol, verbose, outInner, outFinal, loop_count) ;
    }

    FVMatrix(const FVMatrix&) = delete;
    FVMatrix& operator=(const FVMatrix&) = delete;

    FVMatrix(FVMatrix&&) = delete;
    FVMatrix& operator=(FVMatrix&&) = delete;

} ;

