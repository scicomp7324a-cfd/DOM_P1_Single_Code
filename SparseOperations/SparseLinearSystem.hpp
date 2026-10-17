# pragma once

# include<iostream>
# include<vector>
# include "SparseMatrix.hpp"
# include "SparseLinearVector.hpp"

/* *******************************************************************
    SparseLinearSystem.hpp - This class contains all the algebraic 
    operations for sparse linear system and residual calculation  
    functions.
******************************************************************* */

class SparseLinearSystem{

    private:

        SparseMatrix A_ ; 
        int n_ ;

        SparseLinearVector X_ ;
        SparseLinearVector b_ ;
        SparseLinearVector r_ ;

    public:

        // Constructor builds the sparse linear system using sparse addressing
        SparseLinearSystem(const SparseAddress& addr): A_(addr),
                                                 n_(A_.dim()), 
                                                 X_(n_),
                                                 b_(n_),
                                                 r_(n_) {}

        // Returns the dimension of the sparse matrix
        int size() const{
            return n_ ;
        }

        // Non-constant reference to the sparse matrix
        SparseMatrix& A(){ return A_ ;}

        // Constant reference to the sparse matrix
        const SparseMatrix& A() const{ return A_ ;}

        // Non-constant access to the solution vector
        SparseLinearVector& X(){ return X_ ;}

        // Constant access to the solution vector
        const SparseLinearVector& X() const{ return X_ ;}

        // Non-constant access to the source vector
        SparseLinearVector& b(){ return b_ ;}

        // Constant access to the source vector
        const SparseLinearVector& b() const{ return b_ ;}

        // Non-constant access to the residual
        SparseLinearVector& r(){ return r_ ;}

        // Constant access to the residual
        const SparseLinearVector& r() const{ return r_ ;}

        // Operator overload for addition of two sparse linear system
        SparseLinearSystem operator+(SparseLinearSystem& L){

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A() + L.A() ;
            for(int i = 0; i < n_; ++i){
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b() + L.b() ;
            Lr.X() = (*this).X() + L.X() ;

            return Lr ;
        }

        // Operator overload for substraction of two sparse linear system
        SparseLinearSystem operator-(SparseLinearSystem& L){

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A() - L.A() ;
            for(int i = 0; i < n_; ++i){
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b() - L.b() ;
            Lr.X() = (*this).X() - L.X() ;

            return Lr ;
            
        }

        // Operator overload for multiplication of sparse linear system with scalar
        SparseLinearSystem operator*(const double s){

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A()*s ;
            for(int i = 0; i < n_; ++i){
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b()*s ;
            Lr.X() = (*this).X()*s ;

            return Lr ;
        }

        // Operator overload for multiplication of sparse matrix with linear vector in the linear system
        SparseLinearVector operator*(SparseLinearVector& x){
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A()*x ;
        }

        // Function to calculate residual
        void calculateResidual(){
            int sizeOfArray = A_.dim() ;
            
            SparseLinearVector b_star_ = (*this).A()*((*this).X_) ;
            for(int i=0; i < sizeOfArray; ++i){
                r_(i) = b_(i) - b_star_(i) ;
            }
        }

        // Function to calculate residual using optimised algorithm
        void calculateResidualOptimised(){
            for(int i=0; i<n_; ++i){
                (*this).r()(i) = (*this).b()(i) ;
                for(int j=(*this).A().addressing().rowArray()[i]; j<(*this).A().addressing().rowArray()[i+1]; ++j){
                    (*this).r()(i) -= (*this).A()(i,(*this).A().addressing().columnIndexArray()[j])*(*this).X()((*this).A().addressing().columnIndexArray()[j]) ;
                }
            }
        }

        // Function to calculate the average of all the elements in x
        double average(SparseLinearVector x){
            int size = x.size() ;
            double sum = 0 ;
            for(int i=0; i<size; ++i){
                sum += x(i) ;
            }
            return sum/size ;
        }

        // Function to return one norm 
        double oneNorm(){
            return (r_.oneNorm()) ;
        }

        // Function to return two norm
        double twoNorm(){
            return (r_.twoNorm()) ;
        }

        // Function to return infinity norm 
        double infinityNorm(){
            return (r_.infinityNorm()) ;
        }

} ;

// Operator to write to mtx files
inline ostream& operator<<(ostream& out, SparseLinearSystem& sys) {
    sys.A().writeMatrix(out);
    return out;
}

// Operator to read the mtx files
inline istream& operator>>(istream& in, SparseLinearSystem& sys) {
    sys.A().readMatrix(in);
    return in;
}
