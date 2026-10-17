# pragma once

# include<iostream>
# include<vector>
# include<cmath>
# include "SparseAddress.hpp"
# include "SparseLinearVector.hpp"

using namespace std ;

/* *******************************************************************
    SparseMatrix.hpp - This class contains all the algebraic 
    operations for sparse matrix and calculation ofother matrix 
    related metrics functions.
******************************************************************* */

class SparseMatrix{

    private:
        
        const SparseAddress& sparseAddr ;
        int n_ ;
        vector<double> A_ ;
        vector<double> d_ ;

    public:

        // Constructs the sparse matrix using the sparse addressing
        SparseMatrix(const SparseAddress& addr): sparseAddr(addr),
                                                 n_(sparseAddr.numEquations()){
            A_.resize(sparseAddr.numOffDiagonals()) ;
            d_.resize(n_) ;
        }

        // Non-constant operator to access elements of sparse matrix

        /*
        double& operator()(int i, int j){
            if(i != j){
                return A_[sparseAddr.index(i,j)];
            }
            return d_[i] ;
        }
        */
        

        double& operator()(int i, int j)
        {
            if (i == j)
            {
                return d_.at(i);
            }

            const int coefficientIndex =
                sparseAddr.index(i, j);

            if (coefficientIndex < 0)
            {
                throw std::runtime_error(
                    "SparseMatrix: coefficient does not exist for row "
                    + std::to_string(i)
                    + ", column "
                    + std::to_string(j)
                );
            }

            return A_.at(coefficientIndex);
        }

        // Constant operator to access elements of sparse matrix
        /*

        const double& operator()(int i, int j) const{
            if(i != j){
                return A_[sparseAddr.index(i,j)];
            }
            return d_[i] ;
        }
        
        */
        
        const double& operator()(int i, int j) const
        {
            if (i == j)
            {
                return d_.at(i);
            }

            const int coefficientIndex =
                sparseAddr.index(i, j);

            if (coefficientIndex < 0)
            {
                throw std::runtime_error(
                    "SparseMatrix: coefficient does not exist for row "
                    + std::to_string(i)
                    + ", column "
                    + std::to_string(j)
                );
            }

            return A_.at(coefficientIndex);
        }

        // Function to return dimension of the sparse matrix
        const int& dim() const{
            return n_ ;
        }

        // Operator overload for adding sparse matrices
        SparseMatrix operator+(SparseMatrix& B){

            SparseMatrix R(sparseAddr) ; 
            
            for(int i=0; i<n_; ++i){
                R(i,i) = (*this)(i,i) + B(i,i) ;
                for(int j=sparseAddr.rowArray()[i]; j<sparseAddr.rowArray()[i+1]; ++j){
                    R(i,sparseAddr.columnIndexArray()[j]) = (*this)(i,sparseAddr.columnIndexArray()[j]) + B(i,sparseAddr.columnIndexArray()[j]) ;
                }
            }

            return R ;
        }

        // Operator overload for substracting sparse matrices
        SparseMatrix operator-(SparseMatrix& B){

            SparseMatrix R(sparseAddr) ; 
            
            for(int i=0; i<n_; ++i){
                R(i,i) = (*this)(i,i) - B(i,i) ;
                for(int j=sparseAddr.rowArray()[i]; j<sparseAddr.rowArray()[i+1]; ++j){
                    R(i,sparseAddr.columnIndexArray()[j]) = (*this)(i,sparseAddr.columnIndexArray()[j]) - B(i,sparseAddr.columnIndexArray()[j]) ;
                }
            }

            return R ;
        }

        // Operator overload for multiplying sparse matrix with scalar
        SparseMatrix operator*(double s){

            SparseMatrix R(sparseAddr) ; 

            for(int i=0; i<n_; ++i){
                R(i,i) = (*this)(i,i)*s ;
                for(int j=sparseAddr.rowArray()[i]; j<sparseAddr.rowArray()[i+1]; ++j){
                    R(i,sparseAddr.columnIndexArray()[j]) = (*this)(i,sparseAddr.columnIndexArray()[j])*s ;
                }
            }

            return R ;

        }

        // Function to check if the sparse matrix is symmetric
        bool symmetric() const{

            for(int i = 0; i < n_; ++i){
                for(int j = sparseAddr.rowArray()[i] ; j < sparseAddr.rowArray()[i+1]; ++j){
                    if((*this)(i,sparseAddr.columnIndexArray()[j]) != (*this)(sparseAddr.columnIndexArray()[j],i)){
                        return false ;
                    }
                }
            }
            return true ;
        }

        // Function to check strict diagonal dominance of a sparse matrix
        bool strictDiagonalDominance() const{

            double sum = 0;
            for(int i=0; i < n_; ++i){
                for(int j=sparseAddr.rowArray()[i]; j < sparseAddr.rowArray()[i+1]; ++j){
                    sum += (*this)(i, sparseAddr.columnIndexArray()[j]) ;
                }
                if(d_[i] <= sum){
                    return false ;
                }
                sum = 0;
            }
            return true ;

        }

        // Function to check diagonal dominance of a sparse matrix
        bool diagonalDominance() const{

            double sum = 0;
            for(int i=0; i < n_; ++i){
                for(int j=sparseAddr.rowArray()[i]; j < sparseAddr.rowArray()[i+1]; ++j){
                    sum += (*this)(i, sparseAddr.columnIndexArray()[j]) ;
                }
                if(d_[i] < sum){
                    return false ;
                }
                sum = 0 ;
            }
            return true ;
        }

        // Function to check if all the diagonal elements of the matrix are positive
        bool allPositiveDiagonal() const{
            for(int i=0; i < n_; ++i){
                if(d_[i] < 0){
                    return false ;
                }
            }
            return true ;
        }
        
        // Function to check if the sparse matric is symmetric positive definite
        bool spd() const{
            if((*this).symmetric() && (*this).allPositiveDiagonal() && (*this).strictDiagonalDominance()){
                return true ;
            }
            return false ;
        }

        // Operator overload to multiply sparse matrix with a linear vector
        SparseLinearVector operator*(SparseLinearVector& x){

            SparseLinearVector b(n_) ;

            for(int i=0; i < n_; ++i){
                b(i) = (*this)(i,i)*x(i) ;
                for(int j = sparseAddr.rowArray()[i]; j < sparseAddr.rowArray()[i+1]; ++j){
                    b(i) += (*this)(i, sparseAddr.columnIndexArray()[j])*x(sparseAddr.columnIndexArray()[j]) ;
                }
            }

            return b ;
        }

        // Function for vector-matrix multiplication using a nominal algorithm
        SparseLinearVector multiplyNominal(SparseLinearVector& x){

            SparseLinearVector b(n_) ;

            for(int i=0; i < n_; ++i){
                b(i) = (*this)(i,i)*x(i) ;
                for(int j = sparseAddr.rowArray()[i]; j < sparseAddr.rowArray()[i+1]; ++j){
                    b(i) += (*this)(i, sparseAddr.columnIndexArray()[j])*x(sparseAddr.columnIndexArray()[j]) ;
                }
            }

            return b ;
        }

        // Function for vector-matrix multiplication using a optimised algorithm
        SparseLinearVector multiplyOptimised(SparseLinearVector& x){

            SparseLinearVector b(n_) ;

            for(int faceIndex: sparseAddr.mesh().internalFaces()){
                array<int,2> cellON = sparseAddr.mesh().faceOwnerNeighbour(faceIndex) ;
                b(cellON[0]) += A_[sparseAddr.ownerSlot()[faceIndex]]*x(cellON[1]) ;
                b(cellON[1]) += A_[sparseAddr.neighbourSlot()[faceIndex]]*x(cellON[0]) ;
            }

            for(int i=0; i < sparseAddr.mesh().numCells(); ++i){
                b(i) += d_[i]*x(i) ;
            }

            return b ;
        } 

        // Constant reference to sparse address
        const SparseAddress& addressing() const{
            return sparseAddr ;
        }

        // Utility function to split a string at white space
        vector<string> splitWhiteSpace(string s) {
            istringstream iss(s);
            vector<string> out;
            for (string tok; iss >> tok; ) out.push_back(tok);
            return out;
        }

        // Utility function to strip white saces at the stat and end of the string
        string stripWhiteSpace(string s){
            auto first = s.find_first_not_of(" \t\r");
            if (first == string::npos) return ""; 
            auto last  = s.find_last_not_of(" \t\r");
            return s.substr(first, last - first + 1);
        }

        // Function to read matrix coefficients from mtx file
        void readMatrix(istream& in){

            string line ;
            vector<string> strVec ;

            while(getline(in,line)){
                
                if(line.find("%") != string::npos){
                    continue ;
                }else{
                    getline(in,line) ;
                    line = stripWhiteSpace(line) ;
                    strVec = splitWhiteSpace(line) ;

                    if(!strVec[0].compare(strVec[1]) || stoi(strVec[0])!=n_) std::__throw_runtime_error("inconsistent matrix dimensions") ;
                    
                    int nonZeros  = stoi(strVec[2]) ;
                    for(int i=0; i<nonZeros; ++i){

                        getline(in,line) ;
                        line = stripWhiteSpace(line) ;
                        strVec = splitWhiteSpace(line) ;

                        (*this)(stoi(strVec[0])-1, stoi(strVec[1])-1) = stod(strVec[2]) ;
                    }
                }
            }
        }

        // Function to write matrix coefficients to mtx file
        void writeMatrix(ostream& out){

            out << "%%MatrixMarket matrix coordinate real general\n" << endl ;
            out << "%Comments" << endl ;
            out << "There are " << sparseAddr.numCoefficients() << " non zero entries in the following matrix" << endl ;
            out << n_ << " " << n_ << " " << sparseAddr.numCoefficients() << endl ;

            for(int i=0; i<n_; ++i){
                out << i+1 << " " << i+1 << " " << (*this)(i,i) << endl ;
                for(int j=sparseAddr.rowArray()[i]; j<sparseAddr.rowArray()[i+1]; ++j){
                    out << i+1 << " " << sparseAddr.columnIndexArray()[j]+1 << " " << (*this)(i,sparseAddr.columnIndexArray()[j]) << endl ;
                }
            }

            out << endl ;

        }

} ;