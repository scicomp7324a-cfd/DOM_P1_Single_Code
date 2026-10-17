# pragma once

# include <iostream>
# include <vector>
# include <array>
# include <math.h>
# include <map>
# include "FaceAddressedMesh3D.hpp"
# include "SparseLinearVector.hpp"
# include "FaceFluxField.hpp"
# include "FVMatrix.hpp"
# include "Field.hpp"
# include "Patch.hpp"
# include "ScalarProperty.hpp"
# include "Patch.hpp"

using namespace std ;

/* *******************************************************************
    GaussGradient.hpp - This class contains the function to compute
    the gauss gradient of pressure at the cell centers required for the 
    calculation of corrected velocity.   
******************************************************************* */

class GaussGradient{

    private:

    FVMatrix& P_ ;
    Field& PInitF_ ;
    const FaceAddressedMesh3D& mesh_ ;
    SparseLinearVector gradX ;
    SparseLinearVector gradY ;
    SparseLinearVector gradZ ;
    map<int, double> boundaryGradients_ ;
    
    public:

    // Constructs the gauss gradient object using the mesh, fvmatrix and initial pressure field
    GaussGradient(FVMatrix& P, Field& PInitF, const FaceAddressedMesh3D& mesh):P_(P),
                                                                               PInitF_(PInitF),
                                                                               mesh_(mesh),
                                                                               gradX(mesh_.numCells()),
                                                                               gradY(mesh_.numCells()),
                                                                               gradZ(mesh_.numCells()){
        getBoundaryGradients() ;
    }

    // non-constant access to the x component of the gradient
    SparseLinearVector& X(){
        return gradX ;
    }

    // constant access to the x component of the gradient
    const SparseLinearVector& X() const{
        return gradX ;
    }

    // non-constant access to the y component of the gradient
    SparseLinearVector& Y(){
        return gradY ;
    }

    // constant access to the y component of the gradient
    const SparseLinearVector& Y() const{
        return gradY ;
    }

    SparseLinearVector& Z(){
        return gradZ ;
    }

    const SparseLinearVector& Z() const{
        return gradZ ;
    }

    // Function to compute the gauss gradient of pressure
    void computeGradient(){

        clearGradients() ;

        for(int i=0; i<mesh_.numCells(); ++i){
            
            double cellVolume = mesh_.cellVolume(i) ;
            const vector<int> faces = mesh_.cellFaces(i)  ;

            for(int faceIndex: faces){

                array<int,2> cellOwnerNeighbour = mesh_.faceOwnerNeighbour(faceIndex) ;
                array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
                double pf ;

                if(cellOwnerNeighbour[0] != i){
                    faceAreaVector[0] = -faceAreaVector[0] ;
                    faceAreaVector[1] = -faceAreaVector[1] ;
                    faceAreaVector[2] = -faceAreaVector[2] ;
                }

                if(cellOwnerNeighbour[1] != -1){
                    
                    pf = (P_.sys().X()(cellOwnerNeighbour[0]) + P_.sys().X()(cellOwnerNeighbour[1]))/2 ;
                }
                else{
                    double df = distanceFaceToCell(faceIndex) ;
                    pf = P_.sys().X()(cellOwnerNeighbour[0]) + boundaryGradients_.at(faceIndex)*df;
                }

                gradX(i) += (pf/cellVolume)*faceAreaVector[0] ;
                gradY(i) += (pf/cellVolume)*faceAreaVector[1] ;
                gradZ(i) += (pf/cellVolume)*faceAreaVector[2] ;
            }
        }
    }

    // Function to clear the gradient at the each iteration of the SIMPLE loop for the recalculation of the gradients
    void clearGradients(){
        for(int i=0; i<mesh_.numCells(); ++i){
            gradX(i) = 0.0 ;
            gradY(i) = 0.0 ;
            gradZ(i) = 0.0 ;
        }
    }

    // Function to extract the boundary gradient values of pressure as it is prescribed as a Neumann Condition
    void getBoundaryGradients(){
        for(Patch p: mesh_.boundaryPatches()){
            map<int, double> faceVal_ = PInitF_.boundary()[p.index()]->boundaryValues() ;
            for(const auto& [faceIndex, val]: faceVal_){
                boundaryGradients_.insert({faceIndex, val}) ;
            }
        }
    }

    // Function to calculate the face to cell distance inside a control volume
    double distanceFaceToCell(const int& faceIndex){
        array<double,3> faceAreaVec = mesh_.faceAreaVector(faceIndex) ;
        double sfMag = sqrt(faceAreaVec[0]*faceAreaVec[0] + faceAreaVec[1]*faceAreaVec[1] + faceAreaVec[2]*faceAreaVec[2]) ;
        faceAreaVec = {faceAreaVec[0]/sfMag, faceAreaVec[1]/sfMag, faceAreaVec[2]/sfMag};
        array<double,3> p1 = mesh_.faceCentre(faceIndex) ;
        array<double,3> p2 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[0]) ;
        double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] + (p1[2]-p2[2])*faceAreaVec[2] ;
        return abs(d);
    }

} ;