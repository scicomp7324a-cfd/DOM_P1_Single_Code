# pragma once

# include<iostream>
# include "FaceAddressedMesh3D.hpp"

using namespace std ;

/* *******************************************************************
    ScalarProperty.hpp -  This class stores a cell-centred scalar 
    property field with one scalar value associated with each mesh cell.
******************************************************************* */

class ScalarProperty{

    private:
    const FaceAddressedMesh3D& mesh_ ;
    vector<double> val_ ;

    public:

    // Constructs a scalar property field and initializes all cell values to the same constant.
    ScalarProperty(const FaceAddressedMesh3D& mesh, double val):mesh_(mesh),
                                              val_(mesh_.numCells(), val){}
    
    // Returns writable access to the scalar value stored at the given cell index.
    double& operator()(int cellIndex){
        return val_[cellIndex] ;
    }

    // Returns read-only access to the scalar value stored at the given cell index.
    const double& operator()(int cellIndex) const{
        return val_[cellIndex] ;
    }
    
} ;

