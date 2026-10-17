# pragma once

# include <iostream>
# include <vector>
# include <map>
# include "BoundaryCondition.hpp"
# include "FVMatrix.hpp"
# include "FaceAddressedMesh3D.hpp"
# include "ScalarProperty.hpp"
# include "FaceFluxField.hpp"

using namespace std ;

/* *******************************************************************
    PeriodicBC.hpp -  This class is built over the boundary 
    condition interface class and redefines the functions to apply 
    periodic value boundary condition in case of diffusion and 
    convection.   
******************************************************************* */


class PeriodicBC: public BoundaryCondition{

    private:
        map<int,double> faceVal_ ;
        FVMatrix& fvMatrix_ ;
        const FaceAddressedMesh3D& mesh_ ;
        FaceFluxField& flux_ ;
        

    public:

    // Constructs the periodic boundary condition object using the mesh, fvmatrix, associated scalar property for diffusion, and boundary values
    PeriodicBC(FVMatrix& fvMatrix, const FaceAddressedMesh3D& mesh, FaceFluxField& flux, map<int,double> faceVal):faceVal_(faceVal),
                             fvMatrix_(fvMatrix),
                             mesh_(mesh),
                             flux_(flux){}

    // Redefinition of the interface applyDiffusion() function to cater to the fixed value case
    void applyDiffusionBC() override{}

    // Redefinition of the interface applyConvection() function to cater to the fixed value case
    void applyConvectionBC() override{

        if (faceVal_.empty()){return;}

        int sampleFaceIndex = faceVal_.begin()->first;
        if(flux_(sampleFaceIndex) > 0){
                for(const auto& [faceIndex, val]: faceVal_){

                int oppositeFaceIndex = mesh_.oppositeFaceIndex()[faceIndex] ;
                int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                int oppositeCellOwner = mesh_.faceOwnerNeighbour(oppositeFaceIndex)[0] ;

                fvMatrix_.addDiagonal(cellOwner, max(flux_(faceIndex),0.0)) ;
                fvMatrix_.addCoefficient(oppositeCellOwner, cellOwner, min(flux_(oppositeFaceIndex),0.0)) ;
            }
        }
    }

    // Redefinition of the interface updateBoundaryValues() function
    void updateBoundaryValues() override{}

    // Redefinition of type name 
    string typeName() const override{
        return "periodic" ;
    }

    // non-constant access function to the boundary face values
    map<int, double>& boundaryValues() override{
        return faceVal_ ;
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

    // Function to calculate the cell to cell distance inside the control volume
    double distanceCellToCell(const int& faceIndex){
        array<double,3> faceAreaVec = mesh_.faceAreaVector(faceIndex) ;
        double sfMag = sqrt(faceAreaVec[0]*faceAreaVec[0] + faceAreaVec[1]*faceAreaVec[1] + faceAreaVec[2]*faceAreaVec[2]) ;
        faceAreaVec = {faceAreaVec[0]/sfMag, faceAreaVec[1]/sfMag, faceAreaVec[2]/sfMag};
        array<double,3> p1 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[0]) ;
        array<double,3> p2 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[1]) ;
        double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] + (p1[2]-p2[2])*faceAreaVec[2] ;
        return abs(d);
    }

};