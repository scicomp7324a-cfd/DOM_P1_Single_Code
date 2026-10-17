# pragma once

# include<iostream>
# include "FaceAddressedMesh3D.hpp"

using namespace std ;

/* *******************************************************************
    ScalarProperty.hpp -  This class stores a cell-centred scalar 
    property field with one scalar value associated with each mesh cell.
******************************************************************* */

class WallProperty{

    private:
    const FaceAddressedMesh3D& mesh_ ;
    const vector<Patch>& boundaryPatches_ ;
    vector<double> val_ ;
    
    public:

    // Constructs a scalar property field and initializes all cell values to the same constant.
    WallProperty(const FaceAddressedMesh3D& mesh, vector<double> val):mesh_(mesh),
                                                                      boundaryPatches_(mesh_.boundaryPatches()),
                                                                      val_(val)
    {
        if (val_.size() != boundaryPatches_.size()){
            throw runtime_error("WallProperty size does not match number of boundary patches");
        }                        
    }
    
    WallProperty(const FaceAddressedMesh3D& mesh):mesh_(mesh),
                                                  boundaryPatches_(mesh_.boundaryPatches()){}
    
    // Returns writable access to the scalar value stored at the given cell index.
    double& operator()(int patchIndex){
        return val_[patchIndex] ;
    }

    // Returns read-only access to the scalar value stored at the given cell index.
    const double& operator()(int patchIndex) const{
        return val_.at(patchIndex) ;
    }
    
} ;

