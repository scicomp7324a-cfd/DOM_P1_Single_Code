# pragma once

# include <iostream>
# include <vector>
# include <array>
# include <math.h>
# include <map>
# include "FaceAddressedMesh3D.hpp"
# include "SparseLinearVector.hpp"
# include "ScalarProperty.hpp"
# include "FaceFluxField.hpp"
# include "FVMatrix.hpp"
# include "Patch.hpp"
# include "Field.hpp"

using namespace std ;

/* *******************************************************************
    Diffusion.hpp - This class contains the diffusion assembly 
    function. 
******************************************************************* */

class Diffusion{

    private:

        FVMatrix& fvMatrix_ ;
        const FaceAddressedMesh3D& mesh_ ;
        Field& fieldInitBC_ ;
        const ScalarProperty& gamma_ ;


    public:

        // Constructs the diffusion term assembler using the target matrix, mesh, boundary condition field data, and diffusivity field.
        Diffusion(FVMatrix& fvMatrix, const FaceAddressedMesh3D& mesh, Field& fieldInitBC, const ScalarProperty& gamma):fvMatrix_(fvMatrix),
                                              mesh_(mesh),
                                              fieldInitBC_(fieldInitBC),
                                              gamma_(gamma){}

        // This function assembles the diagonal, off-diagonal and boundary contributions arising out of diffusion
        void assemble(){

            for(int faceIndex: mesh_.internalFaces()){

                array<int,2> cellOwnerNeighbour = mesh_.faceOwnerNeighbour(faceIndex) ;
                array<double,3> faceAreaVector = mesh_.faceAreaVector(faceIndex) ;
                double Sf = sqrt(faceAreaVector[0]*faceAreaVector[0] + faceAreaVector[1]*faceAreaVector[1] + faceAreaVector[2]*faceAreaVector[2]) ;
                double gammaF = (gamma_(cellOwnerNeighbour[0]) + gamma_(cellOwnerNeighbour[1]))/2 ;
                double df = distanceCellToCell(faceIndex) ;
                double Df = gammaF*Sf/df ;

                fvMatrix_.addDiagonal(cellOwnerNeighbour[0], Df) ;
                fvMatrix_.addCoefficient(cellOwnerNeighbour[0], cellOwnerNeighbour[1], -Df) ;
                fvMatrix_.addDiagonal(cellOwnerNeighbour[1], Df) ;
                fvMatrix_.addCoefficient(cellOwnerNeighbour[1], cellOwnerNeighbour[0], -Df) ;

            }

            for(Patch p: mesh_.boundaryPatches()){
                fieldInitBC_.boundary()[p.index()]->updateBoundaryValues() ;
                fieldInitBC_.boundary()[p.index()]->applyDiffusionBC() ;
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