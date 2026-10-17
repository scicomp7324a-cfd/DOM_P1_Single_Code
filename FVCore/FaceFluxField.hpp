# pragma once

# include<iostream>
# include<vector>
# include<array>
# include "SparseAddress.hpp"
# include "SparseLinearVector.hpp"
# include "Patch.hpp"
# include "Field.hpp"
# include "FVMatrix.hpp"
# include "OrdinateField.hpp"

using namespace std ;

/* *******************************************************************
    FaceFluxField.hpp -  This class contains all the functions to 
    manage face fluxes in the mesh including initialisation, update
    and correction of fluxes. 
******************************************************************* */

class FaceFluxField{

    private:
    const SparseAddress& addr_ ;
    vector<double> flux_ ;
    const FaceAddressedMesh3D& mesh_ ;
    
    public:

    // Constructs the face flux field object using the sparse addressing, flux values and mesh
    FaceFluxField(const SparseAddress& addr):addr_(addr),
                                       flux_(addr_.mesh().numFaces(), 0.0),
                                       mesh_(addr_.mesh())
    {}

    // Operaor overload to access flux value at a given face index.
    double& operator()(int faceIndex){
        return flux_[faceIndex] ;
    }

    const double& operator() (int faceIndex) const{
        return flux_[faceIndex] ;
    }

    // Function to update initial flux value to the faces depending upon the initial velociy components
    void updateInitialFlux(Field& Ux, Field& Uy, Field& Uz){

        array<int,2> arr ;
        array<double,3> faceVec ;
        double Uxf, Uyf, Uzf ;

        for(int faceIndex: mesh_.internalFaces()){
            arr = mesh_.faceOwnerNeighbour(faceIndex) ;
            faceVec = mesh_.faceAreaVector(faceIndex) ;
            Uxf = (Ux.internal()(arr[0]) + Ux.internal()(arr[1]))/2 ;
            Uyf = (Uy.internal()(arr[0]) + Uy.internal()(arr[1]))/2 ; 
            Uzf = (Uz.internal()(arr[0]) + Uz.internal()(arr[1]))/2 ;
            flux_[faceIndex] = Uxf*faceVec[0] + Uyf*faceVec[1] + Uzf*faceVec[2];
        }

        updateBoundaryFlux(Ux, Uy, Uz) ;

    }

    // Function to update fixed ordinate fluxes based on ordinate directions (Required in DOM)
    void updateInitialFlux(OrdinateS S){

        array<double,3> faceVec ;
        for(int faceIndex: mesh_.internalFaces()){
            faceVec = mesh_.faceAreaVector(faceIndex) ;
            flux_[faceIndex] = (S.x)*faceVec[0] + (S.y)*faceVec[1] + (S.z)*faceVec[2];
        }

        updateBoundaryFlux(S) ;
    }

    // Function to update flux values in the SIMPLE loop after each iteration 
    void updateFlux(const SparseLinearVector& uX, const SparseLinearVector& uY, const SparseLinearVector& uZ){

        array<int,2> arr ;
        array<double,3> faceVec ;
        double Uxf, Uyf, Uzf ;

        for(int faceIndex: mesh_.internalFaces()){
            arr = mesh_.faceOwnerNeighbour(faceIndex) ;
            faceVec = mesh_.faceAreaVector(faceIndex) ;
            Uxf = (uX(arr[0]) + uX(arr[1]))/2 ;
            Uyf = (uY(arr[0]) + uY(arr[1]))/2 ; 
            Uzf = (uZ(arr[0]) + uZ(arr[1]))/2 ;
            flux_[faceIndex] = Uxf*faceVec[0] + Uyf*faceVec[1] + Uzf*faceVec[2];
        }
    }

    // Function to update flux values at the boundary faces
    void updateBoundaryFlux(Field& Ux, Field& Uy, Field& Uz){

        array<double,3> faceVec ;

        for(Patch p: mesh_.boundaryPatches()){

            map<int, double> faceValX = Ux.boundary()[p.index()]->boundaryValues() ;
            map<int, double> faceValY = Uy.boundary()[p.index()]->boundaryValues() ;
            map<int, double> faceValZ = Uz.boundary()[p.index()]->boundaryValues() ;

            for(int faceIndex: p.faces()){
                faceVec = mesh_.faceAreaVector(faceIndex) ;
                flux_[faceIndex] = faceValX.at(faceIndex)*faceVec[0] + faceValY.at(faceIndex)*faceVec[1] + faceValZ.at(faceIndex)*faceVec[2];
            }
        }
    }

    // Function to update fixed flux values at the boundary faces (Required in DOM)
    void updateBoundaryFlux(OrdinateS S){

        array<double,3> faceVec ;

        for(Patch p: mesh_.boundaryPatches()){
            for(int faceIndex: p.faces()){
                faceVec = mesh_.faceAreaVector(faceIndex) ;
                flux_[faceIndex] = (S.x)*faceVec[0] + (S.y)*faceVec[1] + (S.z)*faceVec[2];
            }
        }
    }

    // Function to correct flux values after the pressure correction step in the SIMPLE loop
    void correctFlux(FVMatrix& Ux_, FVMatrix& Uy_, FVMatrix& Uz_, FVMatrix& P_){

        for(int faceIndex: mesh_.internalFaces()){

            array<int,2> cellOwnerNeighbour = mesh_.faceOwnerNeighbour(faceIndex) ;
            array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
            double df = distanceCellToCell(faceIndex) ;
            double sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1] + faceAreaVector[2]*faceAreaVector[2]) ;

            double invAxf = (1/Ux_.sys().A()(cellOwnerNeighbour[0],cellOwnerNeighbour[0]) + 1/Ux_.sys().A()(cellOwnerNeighbour[1], cellOwnerNeighbour[1]))*0.5 ;
            double invAyf = (1/Uy_.sys().A()(cellOwnerNeighbour[0],cellOwnerNeighbour[0]) + 1/Uy_.sys().A()(cellOwnerNeighbour[1], cellOwnerNeighbour[1]))*0.5 ;
            double invAzf = (1/Uz_.sys().A()(cellOwnerNeighbour[0],cellOwnerNeighbour[0]) + 1/Uz_.sys().A()(cellOwnerNeighbour[1], cellOwnerNeighbour[1]))*0.5 ;

            double sfx2 = faceAreaVector[0]*faceAreaVector[0] ;
            double sfy2 = faceAreaVector[1]*faceAreaVector[1] ;
            double sfz2 = faceAreaVector[2]*faceAreaVector[2] ; 
            
            double Df = (sfx2*invAxf + sfy2*invAyf + sfz2*invAzf)/(df*sf) ;
            double dPn = P_.sys().X()(cellOwnerNeighbour[1]) - P_.sys().X()(cellOwnerNeighbour[0]) ;
            flux_[faceIndex] = flux_[faceIndex] - Df*dPn ;
        }
    }

    // Utility function to calculate distance between the centre of two cells
    double distanceCellToCell(const int& faceIndex){

        array<double,3> faceAreaVec = mesh_.faceAreaVector(faceIndex) ;
        double sfMag = sqrt(faceAreaVec[0]*faceAreaVec[0] + faceAreaVec[1]*faceAreaVec[1] + faceAreaVec[2]*faceAreaVec[2]) ;
        faceAreaVec = {faceAreaVec[0]/sfMag, faceAreaVec[1]/sfMag, faceAreaVec[2]/sfMag};
        array<double,3> p1 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[0]) ;
        array<double,3> p2 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[1]) ;
        double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] + (p1[2]-p2[2])*faceAreaVec[2] ;
        return abs(d);
        
    }

} ;