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

using namespace std ;

/* *******************************************************************
    Convection.hpp - This class contains the convection assembly 
    function using the upwaind differencing scheme. 
******************************************************************* */

class Convection{

    private:
    FVMatrix& fvMatrix_ ;
    FaceFluxField& flux_ ;
    Field& fieldInitBC_ ;
    const FaceAddressedMesh3D& mesh_ ;
    
    public:
    
    // Constructs the convection term assembler using the target matrix, face flux field, and boundary condition field data.
    Convection(FVMatrix& fvMatrix, FaceFluxField& flux, Field& fieldInitBC):fvMatrix_(fvMatrix),
                                                        flux_(flux),
                                                        fieldInitBC_(fieldInitBC),
                                                        mesh_(fvMatrix.sys().A().addressing().mesh()){}

    
    // The function carries out upwind assembly using faces. 
    void assemble(){

        for(int faceIndex: mesh_.internalFaces()){
            array<int,2> cellOwnerNeighbour = mesh_.faceOwnerNeighbour(faceIndex) ;
            fvMatrix_.addDiagonal(cellOwnerNeighbour[0], max(flux_(faceIndex),0.0)) ;
            fvMatrix_.addCoefficient(cellOwnerNeighbour[0], cellOwnerNeighbour[1], min(flux_(faceIndex),0.0)) ;
            fvMatrix_.addDiagonal(cellOwnerNeighbour[1], max(-flux_(faceIndex),0.0)) ;
            fvMatrix_.addCoefficient(cellOwnerNeighbour[1], cellOwnerNeighbour[0], min(-flux_(faceIndex),0.0)) ;
        }

        for(Patch p: mesh_.boundaryPatches()){
            fieldInitBC_.boundary()[p.index()]->updateBoundaryValues() ;
            fieldInitBC_.boundary()[p.index()]->applyConvectionBC() ;
        }

    }


} ;