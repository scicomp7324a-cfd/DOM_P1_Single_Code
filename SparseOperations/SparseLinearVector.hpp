# pragma once

# include<iostream>
# include<vector>
# include<cmath>
# include "SparseAddress.hpp"

/* *******************************************************************
    SparseLinearVector.hpp - This class contains all the algebraic 
    operations for sparse linear vector and norm calculation  
    functions.
******************************************************************* */

class SparseLinearVector{

    private:
        int n_ ;
        vector<double> x_ ;

    public:

        // Constructs the sparse linear vector using the given dimension
        SparseLinearVector(int n): n_(n){
            x_.resize(n_) ;
        }

        // Non-constant operator overload to access an element in the vector 
        double& operator()(int j){
            return x_[j] ;
        }

        // Constant operator overload to access an element in the vector
        const double& operator()(int j) const{
            return x_[j] ;
        }

        // Operator overload to return the whole linear vector
        vector<double>& operator()(){
            return x_ ;
        }


        // Operator overload to add two vectors
        inline SparseLinearVector operator+(SparseLinearVector& y_){

            if (n_ != y_.size()) throw std::runtime_error("Vector size mismatch");

            SparseLinearVector z_(n_) ;
            for(int i=0; i < n_; ++i){
                z_(i) = (*this)(i) + y_(i) ;
            }

            return z_ ;
        }

        // Operator overload to substract two vectors
        inline SparseLinearVector operator-(SparseLinearVector& y_){

            if (n_ != y_.size()) throw std::runtime_error("Vector size mismatch");

            SparseLinearVector z_(n_) ;
            for(int i=0; i < n_; ++i){
                z_(i) = (*this)(i) - y_(i) ;
            }

            return z_ ;
        }

        // Operator overload to multiply linear vector to scalar
        inline SparseLinearVector operator*(double s){

            SparseLinearVector z_(n_) ;
            for(int i=0; i < n_; ++i){
                z_(i) = (*this)(i)*s ;
            }
            return z_ ;
        }

        // Operator overload for element wise multiplication of sparse linear vector 

        inline SparseLinearVector operator*(SparseLinearVector& slv){

            SparseLinearVector z_(n_) ;
            for(int i=0; i<n_; ++i){
                z_(i) = (*this)(i)*slv(i) ;
            }
            return z_ ;
        }

        // Function to calculate one norm of the linear vector
        double oneNorm(){
            double sum = 0;
            for(int i = 0; i < n_; ++i){
                sum += abs((*this)(i)) ;
            }
            return sum ;
        }

        // Function to calculate two norm of the linear vector
        double twoNorm(){
            double sum = 0;
            for(int i=0; i< n_ ; ++i){
                sum += pow((*this)(i),2) ;
            }     
            sum = sqrt(sum) ;
            return sum ;
        }

        // Function to calculate infinity norm of the linear vector
        double infinityNorm(){
            double max = 0 ;
            for(int i=0; i< n_; ++i){
                if(abs((*this)(i)) > max){
                    max = abs((*this)(i)) ;
                }
            }
            return max ;
        }

        // Function to return size of linear vector
        int size() const{
            return n_ ;
        }

        // Function to return the maximum absolute value of an element in the linear vector 
        double maxAbsVal(){
            double max = 1e-8 ;
            for(int i=0; i<(*this).size(); ++i){
                if(max < abs((*this)(i))){
                    max = abs((*this)(i)) ;
                }
            }
            return max ;
        }

        // Function to return the minimum absolute value of an element in the linear vector 
        double minAbsVal(){
            double min = 1e8 ;
            for(int i=0; i<(*this).size(); ++i){
                if(min > abs((*this)(i))){
                    min = abs((*this)(i)) ;
                }
            }
            return min ;
        }

        // Function to return the maximum value of an element in the linear vector 
        double maxVal(){
            double max = 1e-8 ;
            for(int i=0; i<(*this).size(); ++i){
                if(max < (*this)(i)){
                    max = (*this)(i) ;
                }
            }
            return max ;
        }

        // Function to return the minimum value of an element in the linear vector 
        double minVal(){
            double min = 1e8 ;
            for(int i=0; i<(*this).size(); ++i){
                if(min > (*this)(i)){
                    min = (*this)(i) ;
                }
            }
            return min ;
        }

        void clear(){
            for(int i=0; i<(*this).size(); ++i){
                (*this)(i) = 0.0 ;
            }
        }

} ;