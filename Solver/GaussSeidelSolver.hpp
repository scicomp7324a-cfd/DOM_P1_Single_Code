# pragma once

# include "PathConfig.hpp"
# include <iostream>
# include <vector>
# include <cmath>
# include <algorithm>
# include <limits>
# include <stdexcept>
# include "SparseLinearSystem.hpp"
# include "GaussSeidelSmoother.hpp"

using namespace std ;

/* *******************************************************************
    GaussSeidelSolver.hpp - This class retains the original solver
    interface used throughout the code, but now solves the sparse
    linear system using ILU(0)-preconditioned BiCGSTAB.
******************************************************************* */

class GaussSeidelSolver: public GaussSeidelSmoother{

    private:

        SparseLinearSystem& sys_ ;
        int maxIter = 10000 ;
        double resInitial = 0.0 ;
        double resFinal = 0.0 ;
        double smallNumber = 1.0e-30 ;
        double pivotTolerance = 1.0e-14 ;

        int numEq_ = 0 ;
        vector<int> rowArray_ ;
        vector<int> columnIndexArray_ ;
        vector<double> diagonal_ ;
        vector<double> offDiagonal_ ;

        vector<double> iluDiagonal_ ;
        vector<double> iluOffDiagonal_ ;

        double dotProduct(const vector<double>& x, const vector<double>& y) const{
            double sum = 0.0 ;
            for(int i=0; i<numEq_; ++i){
                sum += x[i]*y[i] ;
            }
            return sum ;
        }

        double twoNorm(const vector<double>& x) const{
            return sqrt(max(dotProduct(x,x),0.0)) ;
        }

        int coefficientPosition(int row, int column) const{

            int first = rowArray_[row] ;
            int last = rowArray_[row+1] ;

            auto iterator = lower_bound(
                columnIndexArray_.begin()+first,
                columnIndexArray_.begin()+last,
                column
            ) ;

            if(iterator == columnIndexArray_.begin()+last || *iterator != column){
                return -1 ;
            }

            return static_cast<int>(iterator-columnIndexArray_.begin()) ;
        }

        void cacheMatrix(){

            numEq_ = sys_.A().dim() ;
            rowArray_ = sys_.A().addressing().rowArray() ;
            columnIndexArray_ = sys_.A().addressing().columnIndexArray() ;

            diagonal_.assign(numEq_,0.0) ;
            offDiagonal_.assign(columnIndexArray_.size(),0.0) ;

            for(int i=0; i<numEq_; ++i){
                diagonal_[i] = sys_.A()(i,i) ;

                for(int coefficient=rowArray_[i]; coefficient<rowArray_[i+1]; ++coefficient){
                    offDiagonal_[coefficient] = sys_.A()(i,columnIndexArray_[coefficient]) ;
                }
            }
        }

        void matrixVectorProduct(const vector<double>& x, vector<double>& result) const{

            fill(result.begin(),result.end(),0.0) ;

            for(int i=0; i<numEq_; ++i){
                double sum = diagonal_[i]*x[i] ;

                for(int coefficient=rowArray_[i]; coefficient<rowArray_[i+1]; ++coefficient){
                    sum += offDiagonal_[coefficient]*x[columnIndexArray_[coefficient]] ;
                }

                result[i] = sum ;
            }
        }

        void calculateResidual(const vector<double>& x, vector<double>& residual) const{

            vector<double> Ax(numEq_,0.0) ;
            matrixVectorProduct(x,Ax) ;

            for(int i=0; i<numEq_; ++i){
                residual[i] = sys_.b()(i)-Ax[i] ;
            }
        }

        double protectedPivot(double pivot) const{

            if(abs(pivot) >= pivotTolerance){
                return pivot ;
            }

            if(pivot < 0.0){
                return -pivotTolerance ;
            }

            return pivotTolerance ;
        }

        void buildILU0(){

            iluDiagonal_ = diagonal_ ;
            iluOffDiagonal_ = offDiagonal_ ;

            for(int row=0; row<numEq_; ++row){

                for(int coefficient=rowArray_[row]; coefficient<rowArray_[row+1]; ++coefficient){

                    int lowerColumn = columnIndexArray_[coefficient] ;

                    if(lowerColumn >= row){
                        break ;
                    }

                    double pivot = protectedPivot(iluDiagonal_[lowerColumn]) ;
                    double lowerFactor = iluOffDiagonal_[coefficient]/pivot ;
                    iluOffDiagonal_[coefficient] = lowerFactor ;

                    for(int upperCoefficient=rowArray_[lowerColumn];
                        upperCoefficient<rowArray_[lowerColumn+1];
                        ++upperCoefficient){

                        int upperColumn = columnIndexArray_[upperCoefficient] ;

                        if(upperColumn <= lowerColumn){
                            continue ;
                        }

                        double correction = lowerFactor*iluOffDiagonal_[upperCoefficient] ;

                        if(upperColumn == row){
                            iluDiagonal_[row] -= correction ;
                        }
                        else{
                            int rowPosition = coefficientPosition(row,upperColumn) ;
                            if(rowPosition >= 0){
                                iluOffDiagonal_[rowPosition] -= correction ;
                            }
                        }
                    }
                }

                iluDiagonal_[row] = protectedPivot(iluDiagonal_[row]) ;
            }
        }

        void applyPreconditioner(const vector<double>& input,
                                 vector<double>& output,
                                 vector<double>& work) const{

            for(int i=0; i<numEq_; ++i){
                double sum = input[i] ;

                for(int coefficient=rowArray_[i]; coefficient<rowArray_[i+1]; ++coefficient){
                    int column = columnIndexArray_[coefficient] ;
                    if(column >= i){
                        break ;
                    }
                    sum -= iluOffDiagonal_[coefficient]*work[column] ;
                }

                work[i] = sum ;
            }

            for(int i=numEq_-1; i>=0; --i){
                double sum = work[i] ;

                for(int coefficient=rowArray_[i]; coefficient<rowArray_[i+1]; ++coefficient){
                    int column = columnIndexArray_[coefficient] ;
                    if(column > i){
                        sum -= iluOffDiagonal_[coefficient]*output[column] ;
                    }
                }

                output[i] = sum/iluDiagonal_[i] ;
            }
        }

    public:

        GaussSeidelSolver(SparseLinearSystem& sys): GaussSeidelSmoother(sys),
                                                    sys_(sys){}

        void solve(double relTol, bool verbose, ostream& outInner, ostream& outFinal, int loop_count){

            cacheMatrix() ;
            buildILU0() ;

            vector<double> x(numEq_,0.0) ;
            vector<double> residual(numEq_,0.0) ;
            vector<double> shadowResidual(numEq_,0.0) ;
            vector<double> p(numEq_,0.0) ;
            vector<double> v(numEq_,0.0) ;
            vector<double> s(numEq_,0.0) ;
            vector<double> t(numEq_,0.0) ;
            vector<double> preconditionedP(numEq_,0.0) ;
            vector<double> preconditionedS(numEq_,0.0) ;
            vector<double> preconditionerWork(numEq_,0.0) ;

            for(int i=0; i<numEq_; ++i){
                x[i] = sys_.X()(i) ;
            }

            calculateResidual(x,residual) ;
            shadowResidual = residual ;

            resInitial = twoNorm(residual) ;
            resFinal = resInitial ;

            if(resInitial <= smallNumber){
                outInner << endl ;

                if(verbose){
                    cout << "Initial Residual: " << resInitial
                         << ", Final Residual: " << resFinal << endl ;
                }

                outFinal << loop_count << "," << resInitial << ","
                         << resFinal << "," << 0.0 << endl ;
                return ;
            }

            double rhoOld = 1.0 ;
            double alpha = 1.0 ;
            double omega = 1.0 ;
            bool converged = false ;
            int iteration = 0 ;

            for(iteration=1; iteration<=maxIter; ++iteration){

                double rhoNew = dotProduct(shadowResidual,residual) ;

                if(abs(rhoNew) <= smallNumber){
                    break ;
                }

                if(iteration == 1){
                    p = residual ;
                }
                else{
                    if(abs(omega) <= smallNumber){
                        break ;
                    }

                    double beta = (rhoNew/rhoOld)*(alpha/omega) ;

                    for(int i=0; i<numEq_; ++i){
                        p[i] = residual[i] + beta*(p[i]-omega*v[i]) ;
                    }
                }

                fill(preconditionedP.begin(),preconditionedP.end(),0.0) ;
                applyPreconditioner(p,preconditionedP,preconditionerWork) ;
                matrixVectorProduct(preconditionedP,v) ;

                double shadowV = dotProduct(shadowResidual,v) ;
                if(abs(shadowV) <= smallNumber){
                    break ;
                }

                alpha = rhoNew/shadowV ;

                for(int i=0; i<numEq_; ++i){
                    s[i] = residual[i]-alpha*v[i] ;
                }

                double sNorm = twoNorm(s) ;

                if(sNorm/resInitial <= relTol){
                    for(int i=0; i<numEq_; ++i){
                        x[i] += alpha*preconditionedP[i] ;
                    }

                    residual = s ;
                    resFinal = sNorm ;
                    converged = true ;

                    outInner << iteration << ","
                             << resInitial << ","
                             << resFinal << ","
                             << resFinal/resInitial << endl ;
                    break ;
                }

                fill(preconditionedS.begin(),preconditionedS.end(),0.0) ;
                applyPreconditioner(s,preconditionedS,preconditionerWork) ;
                matrixVectorProduct(preconditionedS,t) ;

                double tDotT = dotProduct(t,t) ;
                if(abs(tDotT) <= smallNumber){
                    break ;
                }

                omega = dotProduct(t,s)/tDotT ;

                for(int i=0; i<numEq_; ++i){
                    x[i] += alpha*preconditionedP[i] + omega*preconditionedS[i] ;
                    residual[i] = s[i]-omega*t[i] ;
                }

                resFinal = twoNorm(residual) ;

                outInner << iteration << ","
                         << resInitial << ","
                         << resFinal << ","
                         << resFinal/resInitial << endl ;

                if(resFinal/resInitial <= relTol){
                    converged = true ;
                    break ;
                }

                if(abs(omega) <= smallNumber){
                    break ;
                }

                rhoOld = rhoNew ;
            }

            for(int i=0; i<numEq_; ++i){
                sys_.X()(i) = x[i] ;
            }

            // Recalculate the true residual using the original matrix object.
            sys_.calculateResidual() ;
            resFinal = sys_.twoNorm() ;

            outInner << endl ;

            if(verbose){
                cout << "Initial Residual: " << resInitial
                     << ", Final Residual: " << resFinal << endl ;

                if(!converged){
                    cout << "BiCGSTAB stopped before reaching the requested tolerance." << endl ;
                }
            }

            outFinal << loop_count << "," << resInitial << ","
                     << resFinal << "," << resFinal/resInitial << endl ;
        }

} ;
