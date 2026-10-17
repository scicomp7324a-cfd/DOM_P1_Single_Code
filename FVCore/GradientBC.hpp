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
    GradientBC.hpp -  This class is built over the boundary 
    condition interface class and redefines the functions to apply 
    gradient value boundary condition in case of diffusion and 
    convection.   
******************************************************************* */

class GradientBC: public BoundaryCondition{

    private:
        map<int,double> faceVal_ ;
        FVMatrix& fvMatrix_ ;
        const FaceAddressedMesh3D& mesh_ ;
        FaceFluxField& flux_ ;
        const ScalarProperty& gamma_ ;

    public:

    // Constructs the gradient boundary condition object using the mesh, flux field, associated scalar field for diffusion and boundary face values
    GradientBC(FVMatrix& fvMatrix, const FaceAddressedMesh3D& mesh, FaceFluxField& flux, const ScalarProperty& gamma, map<int,double> faceVal):faceVal_(faceVal),
                             fvMatrix_(fvMatrix),
                             mesh_(mesh),
                             flux_(flux),
                             gamma_(gamma){}


    // Redefinition of the interface applyDiffusion() function to cater to the gradient case
    void applyDiffusionBC() override{
        for(const auto& [faceIndex, val]: faceVal_){
            int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
            array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
            double Sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1] + faceAreaVector[2]*faceAreaVector[2]) ;
            double gammaF = gamma_(cellOwner) ;
            double Df = Sf*gammaF*val ;
            fvMatrix_.addToSource(cellOwner, Df) ;
        }
    }

    // Redefinition of the interface applyConvection() function to cater to the gradient case
    void applyConvectionBC() override{
        for(const auto& [faceIndex, val]: faceVal_){
            array<int,2> cellOwner = mesh_.faceOwnerNeighbour(faceIndex) ;
            double d = distanceFaceToCell(faceIndex) ;
            if(flux_(faceIndex) >= 0){
                fvMatrix_.addDiagonal(cellOwner[0], flux_(faceIndex)) ;
            }
            else{
                fvMatrix_.addDiagonal(cellOwner[0], flux_(faceIndex)) ;
                fvMatrix_.addToSource(cellOwner[0], -flux_(faceIndex)*val*d) ;
            }
        }
    }

    // Redefinition of the interface updateBoundaryValues() function
    void updateBoundaryValues() override{}

    // Redefinition of type name 
    string typeName() const override{
        return "fixedgradient" ;
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
        double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] + (p1[2]-p2[2])*faceAreaVec[2];
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