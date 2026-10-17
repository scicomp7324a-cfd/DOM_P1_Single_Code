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
    FixedValueBC.hpp -  This class is built over the boundary 
    condition interface class and redefines the functions to apply 
    fixed value boundary condition in case of diffusion and 
    convection.   
******************************************************************* */

class MarshakRadiationBC: public BoundaryCondition{

    private:
        map<int,double> faceVal_ ;
        map<int,double> wallTemperature_ ;
        FVMatrix& fvMatrixG_ ;
        const FaceAddressedMesh3D& mesh_ ;
        FaceFluxField& flux_ ;
        const ScalarProperty& gamma_ ;
        const double emissivity_ ;
        double sigma = 5.670374419e-8  ;

    public:

    // Constructs the fixed value boundary condition object using the mesh, fvmatrix, associated scalar property for diffusion, and boundary values
    MarshakRadiationBC(FVMatrix& fvMatrixG, const FaceAddressedMesh3D& mesh, FaceFluxField& flux, const ScalarProperty& gamma, double emissivity, map<int,double> faceVal, map<int,double> wallTemperature):faceVal_(faceVal),
                             wallTemperature_(std::move(wallTemperature)),
                             fvMatrixG_(fvMatrixG),
                             mesh_(mesh),
                             flux_(flux),
                             gamma_(gamma),
                             emissivity_(emissivity){}

    // Redefinition of the interface applyDiffusion() function to cater to the fixed value case
    /*
    void applyDiffusionBC() override{
        for(const auto& [faceIndex, val]: faceVal_){
            //cout << "Face Index: " << faceIndex << endl ;
            //cout << "Value: " << val << endl ;
            int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
            array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
            double Sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1] + faceAreaVector[2]*faceAreaVector[2]) ;
            double gammaF = gamma_(cellOwner) ;
            double df = distanceFaceToCell(faceIndex) ;
            double Df = gammaF*Sf/df ;
            fvMatrixG_.addDiagonal(cellOwner, Df) ;
            fvMatrixG_.addToSource(cellOwner, Df*val) ;
        }
    }
    */

    void applyDiffusionBC() override{

        for(const auto& [faceIndex,val]: faceVal_){
            int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
            array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
            double Sf = sqrt(
                faceAreaVector[0]*faceAreaVector[0]
                + faceAreaVector[1]*faceAreaVector[1]
                + faceAreaVector[2]*faceAreaVector[2]
            ) ;
            double gammaF = gamma_(cellOwner) ;
            double df = distanceFaceToCell(faceIndex) ;
            double Bw = gammaF/df ;
            double Ew = emissivity_/(2*(2-emissivity_)) ;
            double coefficient = Sf*Bw*Ew/(Bw+Ew) ;

            fvMatrixG_.addDiagonal(cellOwner,coefficient) ;
            fvMatrixG_.addToSource(cellOwner,coefficient*4*sigma*powl(wallTemperature_[faceIndex],4)) ;
        }
    }

    // Redefinition of the interface applyConvection() function
    void applyConvectionBC() override{}

    // Redefinition of the interface updateBoundaryValues() function 
    void updateBoundaryValues() override{
        for(auto& [faceIndex, val]: faceVal_){
            int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
            double Ew = emissivity_/(2*(2-emissivity_)) ;
            double df = distanceFaceToCell(faceIndex) ;
            double Bw = gamma_(cellOwner)/df ;
            faceVal_[faceIndex] = (Bw*fvMatrixG_.sys().X()(cellOwner) + 4*Ew*sigma*powl(wallTemperature_[faceIndex],4))/(Bw + Ew) ;
        }
    }

    // non-constant access function to the boundary face values
    map<int, double>& boundaryValues() override{
        return faceVal_ ;
    }

    // Redefinition of type name 
    string typeName() const override{
        return "marshak" ;
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