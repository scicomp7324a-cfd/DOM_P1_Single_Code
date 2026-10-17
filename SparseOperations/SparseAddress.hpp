# pragma once

# include<iostream>
# include<vector>
# include<map>
# include "FaceAddressedMesh3D.hpp"

using namespace std ;

/* *******************************************************************
    SparseAddress.hpp -  This class provides the functions to build 
    and access the modified compressed row storage format (CSR) row 
    and column arrays. It also builds the index mappings required to 
    access the position of a coefficient in the sparse matrix. 
******************************************************************* */

class SparseAddress{

    private:

    const FaceAddressedMesh3D mesh_ ;
    int numEquations_ ;
    int numOffDiagonals_ ;
    int numCoefficients_ ;
    vector<int> rowArray_ ;
    vector<int> colIndexArray_ ;
    vector<int> diagonalCoeff_ ;
    map<array<int,2>, int> indexMap_ ;

    // Optimisation related array
    vector<int> ownerSlot_ ;
    vector<int> neighbourSlot_ ;

    public:

    // This constructor builds up the row, column index and slot arrays 
    SparseAddress(vector<string>& fileNames): mesh_(fileNames),
                                             numEquations_(mesh_.numCells())
    {
        buildRowAndColIndexArray() ;
        numOffDiagonals_ = colIndexArray_.size() ;
        numCoefficients_ = colIndexArray_.size() + numEquations_;

        ownerSlot_.resize(mesh_.numFaces()) ;
        neighbourSlot_.resize(mesh_.numFaces()) ;
        buildSlotArrays() ;
      
    }

    // Returns the number of equations i.e. number of rows and columns in the equation.
    const int& numEquations() const{
        return numEquations_ ;
    }

    // Returns the number of off diagonal coefficients in the matrix.
    const int& numOffDiagonals() const{
        return numOffDiagonals_ ;
    }

    // Returns the number of total coefficients in the matrix.
    const int& numCoefficients() const{
        return numCoefficients_ ;
    }

    // Returns the row array which records the start and end of each row in the column array 
    inline const vector<int>& rowArray() const{
        return rowArray_ ;
    }

    // Returns the column array which records the column index for each coefficients
    inline const vector<int>& columnIndexArray() const{
        return colIndexArray_ ;
    }

    // Returns reference to owner slot array for optimised multiplication
    inline const vector<int>& ownerSlot() const{
        return ownerSlot_ ;
    }

    // Returns reference to neighbour slot array for optimised multiplication
    inline const vector<int>& neighbourSlot() const{
        return neighbourSlot_ ;
    }

    // Return the position of a coefficient in the sparse matrix
    int index(int i, int j) const{
        try{
            return indexMap_.at({i,j}) ;
        }
        catch(exception& e){
            return -1 ;
        }
    }

    // Returns a reference to the mesh 
    const FaceAddressedMesh3D& mesh() const{
        return mesh_ ;
    }

    // Builds the column array, row array and index mapping
    /*
    void buildRowAndColIndexArray(){
        vector<int> neighbours ;
        rowArray_.push_back(0) ;
        int count = 0;
        for(int i = 0; i < numEquations_; ++i){
            neighbours = mesh_.cellNeighbours(i) ;
            for(int var: neighbours){
                if(var != -1){
                    colIndexArray_.push_back(var) ;
                    int size = colIndexArray_.size() ;
                    indexMap_.insert({{i,var}, size-1}) ;
                    ++count ;
                } 
            }
            rowArray_.push_back(count) ;
        }
    }
    */

    void buildRowAndColIndexArray(){
        rowArray_.clear();
        colIndexArray_.clear();
        indexMap_.clear();

        rowArray_.push_back(0);

        int coefficientIndex = 0;

        for (int cellIndex = 0;
            cellIndex < numEquations_;
            ++cellIndex)
        {
            set<int> rowColumns;

            const vector<int> ordinaryNeighbours =
                mesh_.cellNeighbours(cellIndex);

            for (int neighbourCell : ordinaryNeighbours)
            {
                if (
                    neighbourCell >= 0
                    && neighbourCell != cellIndex
                )
                {
                    rowColumns.insert(neighbourCell);
                }
            }

            for (const Patch& patch : mesh_.boundaryPatches())
            {
                for (int faceIndex : patch.faces())
                {
                    const int ownerCell =
                        mesh_.faceOwnerNeighbour(faceIndex)[0];

                    if (ownerCell != cellIndex)
                    {
                        continue;
                    }

                    const int oppositeFaceIndex =
                        mesh_.oppositeFaceIndex()[faceIndex];

                    if (oppositeFaceIndex < 0)
                    {
                        continue;
                    }

                    const int oppositeOwnerCell =
                        mesh_.faceOwnerNeighbour(
                            oppositeFaceIndex
                        )[0];

                    if (
                        oppositeOwnerCell >= 0
                        && oppositeOwnerCell != cellIndex
                    )
                    {
                        rowColumns.insert(oppositeOwnerCell);
                    }
                }
            }

            for (int columnCell : rowColumns)
            {
                colIndexArray_.push_back(columnCell);

                indexMap_[{cellIndex, columnCell}] =
                    coefficientIndex;

                ++coefficientIndex;
            }

            rowArray_.push_back(coefficientIndex);
        }
    }

    // Function to build owner and neighbour slot arrays for optimised multiplication
    void buildSlotArrays(){
        for(int faceIndex: mesh_.internalFaces()){
            array<int,2> cellON = mesh_.faceOwnerNeighbour(faceIndex) ;
            ownerSlot_[faceIndex] = index(cellON[0], cellON[1]) ;
            neighbourSlot_[faceIndex] = index(cellON[1], cellON[0]) ;
        }
    }

    // Returns the matrix band  
    int matrixBand(){
        int colIndex ;
        int matrixBand_ ;
        for(size_t i=0; i<rowArray_.size(); ++i){
            for(int j = rowArray_[i]; j < rowArray_[i+1]; ++j){
                colIndex = colIndexArray_[j] ;
                if(abs(colIndex- static_cast<int>(i)) > matrixBand_){
                    matrixBand_ = abs(colIndex-static_cast<int>(i)) ;
                }
            }
        }
        return matrixBand_ ;
    }
}; 