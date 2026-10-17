# pragma once 

# include <iostream>
# include <vector>

using namespace std ;

/* *******************************************************************
    Patch.hpp - This class stores the identity and face connectivity 
    of a boundary patch in the mesh.
******************************************************************* */


class Patch{

    private:

        string patchName_ ;
        int patchIndex_ ;
        const vector<int> patchFaces_ ;
        const int nFaces_ ;

    public:

        // Constructs a patch from its name, index, list of face indices, and total number of faces.
        Patch(string patchName, int patchIndex, vector<int> patchFaces, int nFaces):patchName_(patchName),
                                                                                     patchIndex_(patchIndex),
                                                                                     patchFaces_(patchFaces),
                                                                                     nFaces_(nFaces){}

        // Returns the name of the patch.                                                                             
        string name(){return patchName_ ;}
        const string name() const{return patchName_ ;}

        // Returns the integer index of the patch.
        int index(){return patchIndex_; }
        
        // Returns the list of face indices that belong to this patch.
        const vector<int>& faces() const{ return patchFaces_ ;}

        // Returns the total number of faces in this patch.
        int numFaces() const{ return nFaces_ ; } 

} ;