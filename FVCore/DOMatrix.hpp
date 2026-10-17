# pragma once

# include <iostream>
# include <array>
# include <vector>
# include <map>
# include <chrono>

# include "PathConfig.hpp"
# include "FVMatrix.hpp"
# include "OrdinateField.hpp"
# include "SparseAddress.hpp"
# include "FaceFluxField.hpp"
# include "Field.hpp"
# include "SparseAddress.hpp"
# include "FaceAddressedMesh3D.hpp"
# include "BoundaryCondition.hpp"
# include "WallProperty.hpp"
# include "SparseLinearVector.hpp"
# include "PathConfig.hpp"

class DOMatrix{

    private:

        const SparseAddress& addr_ ;
        const FaceAddressedMesh3D& mesh_ ;
        const PathConfig& paths_ ;

        FVMatrix fvMatrix_ ;
        FaceFluxField flux_ ;
        OrdinateS S_ ;
        SparseLinearVector sourceDiagonalI_ ;
        SparseLinearVector sourceTermI_ ;

        Field intensityField_ ;

        int ordinateNumber_ ;

        ofstream outInnerI_ ;
        ofstream outFinalI_ ;
        ofstream outXI_ ;

    public:

        DOMatrix(const SparseAddress& addr, const FaceAddressedMesh3D& mesh, const PathConfig& paths, OrdinateS S, int ordinateNumber):
                                                         addr_(addr),
                                                         mesh_(mesh),
                                                         paths_(paths),
                                                         fvMatrix_(addr),
                                                         flux_(addr),
                                                         S_(S),
                                                         sourceDiagonalI_(mesh_.numCells()),
                                                         sourceTermI_(mesh_.numCells()),
                                                         ordinateNumber_(ordinateNumber)

                                                         
        {
            flux_.updateInitialFlux(S) ;
        }

        FVMatrix& I(){ return fvMatrix_ ;}
        const FVMatrix& I() const{ return fvMatrix_ ;} 

        OrdinateS& S(){ return S_ ;}
        const OrdinateS& S() const{ return S_ ;}

        FaceFluxField& F(){ return flux_ ;}
        const FaceFluxField& F() const{ return flux_ ;}

        Field& iabf(){ return intensityField_ ;}
        const Field& iabf() const{return intensityField_;}

        SparseLinearVector& sourceDiagonalI(){ return sourceDiagonalI_ ;}
        const SparseLinearVector& sourceDiagonalI() const{ return sourceDiagonalI_ ;}

        SparseLinearVector& sourceTermI(){ return sourceTermI_ ;}
        const SparseLinearVector& sourceTermI() const{ return sourceTermI_ ;}

        ofstream& outInnerI(){ return outInnerI_ ;}
        const ofstream& outInnerI() const{ return outInnerI_ ;}

        ofstream& outFinalI(){ return outFinalI_ ;}
        const ofstream& outFinalI() const{ return outFinalI_ ;}
        
        ofstream& outXI(){ return outXI_ ;}
        const ofstream& outXI() const{ return outXI_ ;}

        void getIntensityField(Field intensityField){
            intensityField_ = std::move(intensityField) ;
        }

        void openFile(string caseName, string subCaseName){

            string innerResPath_ = "innerResidualI" + to_string(ordinateNumber_) + ".csv" ;
            outInnerI_.open(paths_.getOutputDirectory(caseName, subCaseName) / innerResPath_) ;
            outInnerI_ << "iteration,initial_residual,final_residual,ratio" << endl ;

            string outerFinalPath_ = "finalRedidualI" + to_string(ordinateNumber_) + ".csv" ;
            outFinalI_.open(paths_.getOutputDirectory(caseName, subCaseName) / outerFinalPath_) ;
            outFinalI_ << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            string XIFinalPath_ = "XIFinal" + to_string(ordinateNumber_) + ".csv" ;
            outXI_.open(paths_.getOutputDirectory(caseName, subCaseName) / XIFinalPath_) ;
        }

        DOMatrix(const DOMatrix&) = delete;
        DOMatrix& operator=(const DOMatrix&) = delete;

        DOMatrix(DOMatrix&&) = delete;
        DOMatrix& operator=(DOMatrix&&) = delete;
     
} ;