# pragma once

# include <iostream>
# include <vector>
# include <map>
# include <memory>

# include "BoundaryCondition.hpp"
# include "FVMatrix.hpp"
# include "FaceAddressedMesh3D.hpp"
# include "ScalarProperty.hpp"
# include "FaceFluxField.hpp"
# include "DOMatrix.hpp"

using namespace std ;

/* *******************************************************************
    DiscreteOrdinateBC.hpp -  This class is built over the boundary 
    condition interface class and redefines the functions to apply 
    DiscreteOrdinate value boundary condition in case of diffusion and 
    convection.   
******************************************************************* */

class DiscreteOrdinateBC: public BoundaryCondition{

    private:
        map<int,double> faceVal_ ;
        FVMatrix& fvMatrix_ ;
        const FaceAddressedMesh3D& mesh_ ;
        FaceFluxField& flux_ ;
        const double emissivity_ ;
        vector<unique_ptr<DOMatrix>>& intensityField_;
        map<int,double>& wallTemperature_ ;

        double sigma_ = 5.670374419e-8  ;
        double pi_ = 3.1415926535897932 ;

 
    public:

    // Constructs the periodic boundary condition object using the mesh, fvmatrix, associated scalar property for diffusion, and boundary values
    DiscreteOrdinateBC(FVMatrix& fvMatrix, const FaceAddressedMesh3D& mesh, FaceFluxField& flux, const double emissivity, vector<unique_ptr<DOMatrix>>& intensityField, map<int,double> faceVal, map<int,double>& wallTemperature):faceVal_(faceVal),
                             fvMatrix_(fvMatrix),
                             mesh_(mesh),
                             flux_(flux),
                             emissivity_(emissivity),
                             intensityField_(intensityField),
                             wallTemperature_(wallTemperature)
                             {}

    // Redefinition of the interface applyDiffusion() function to cater to the fixed value case
    void applyDiffusionBC() override{}

    // Redefinition of the interface applyConvection() function to cater to the fixed value case
    void applyConvectionBC() override {
        for(const auto& [faceIndex, val]: faceVal_){
            array<int,2> cellOwner = mesh_.faceOwnerNeighbour(faceIndex) ;
            fvMatrix_.addDiagonal(cellOwner[0], max(flux_(faceIndex), 0.0)) ;
            fvMatrix_.addToSource(cellOwner[0], -min(flux_(faceIndex), 0.0)*val) ;
        }
    }

    // Redefinition of the interface updateBoundaryValues() function
    void updateBoundaryValues() override {
        for(auto& [faceIndex, val]: faceVal_){
            int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
            faceVal_[faceIndex] = emissivity_*(sigma_*pow(wallTemperature_[faceIndex],4)/pi_) ;
            for(const unique_ptr<DOMatrix>& intensityPointer : intensityField_){
                DOMatrix& IF = *intensityPointer ;
                faceVal_[faceIndex] += ((1-emissivity_)/pi_)* max(IF.F()(faceIndex)/mesh_.faceArea(faceIndex), 0.0) * (IF.S().w) * (IF.I().sys().X()(cellOwner)) ;
            }
        }
    }

    // Redefinition of type name 
    string typeName() const override {
        return "discreteordinate" ;
    }

    // non-constant access function to the boundary face values
    map<int, double>& boundaryValues() override {
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