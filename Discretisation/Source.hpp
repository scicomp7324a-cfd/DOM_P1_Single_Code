#pragma once

#include <iostream>
#include <map>
#include <array>
#include "FaceAddressedMesh3D.hpp"
#include "FVMatrix.hpp"
#include "SparseLinearVector.hpp"

using namespace std;

/* *******************************************************************
    Source.hpp - This class contains the source assembly 
    function. 
******************************************************************* */

class Source
{

private:

    FVMatrix &fvMatrix_;
    const FaceAddressedMesh3D &mesh_;
    const SparseLinearVector &diagonalTerm_ ;
    const SparseLinearVector &sourceTerm_ ;
    

public:
    
    // Constructs the source term assembler using the target matrix, mesh, and per-cell source coefficients.
    Source(FVMatrix &fvMatrix, const FaceAddressedMesh3D &mesh, const SparseLinearVector &diagonalTerm, const SparseLinearVector &sourceTerm) : fvMatrix_(fvMatrix),
                                                            mesh_(mesh),
                                                            diagonalTerm_(diagonalTerm),
                                                            sourceTerm_(sourceTerm){}

    // This function assembles. Source is defined as Su - Sp*PhiP
    void assemble()
    {
        for (int i = 0; i < mesh_.numCells(); ++i)
        {
            double volume = mesh_.cellVolume(i) ;
            fvMatrix_.addDiagonal(i, diagonalTerm_(i)*volume) ;
            fvMatrix_.addToSource(i, sourceTerm_(i)*volume) ;
        }
    }

};