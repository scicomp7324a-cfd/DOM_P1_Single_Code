# pragma once

# include<iostream>
# include <string>
# include <vector>
# include <array>
# include <set>
# include <map>
# include <ostream>
# include <algorithm>
# include <cmath>

# include "MeshFileReader3D.hpp"
# include "Patch.hpp"


using namespace std ;

class FaceAddressedMesh3D{

    private:

        const MeshFileReader3D mfr ;
        const vector<array<double,3>> points_ ;
        const vector<vector<int>> faces_ ;
        const vector<vector<int>> cells_ ;
        const vector<Patch> boundaryPatches_ ;

        const int numPoints_ ;
        const int numFaces_ ;
        const int numCells_ ;
        const int numBoundaryPatches_ ;
 
        vector<array<int,2>> faceOwnerNeighbour_ ;
        vector<int> internalFaces_ ;
        vector<int> boundaryFaces_ ;
        vector<int> oppositeFaceIndex_ ;

        vector<int> numCellXYZ_ ;
        vector<double> dimensionsXYZ_ ;

    public:

        FaceAddressedMesh3D(vector<string> fileNames): mfr(fileNames),
                                                       points_(mfr.points()),
                                                       faces_(mfr.faces()),
                                                       cells_(mfr.cells()),
                                                       boundaryPatches_(mfr.boundaryPatches()),
                                                       numPoints_(mfr.numPoints()),
                                                       numFaces_(mfr.numFaces()),
                                                       numCells_(mfr.numCells()),
                                                       numBoundaryPatches_(mfr.numBoundaryPatches()),
                                                       numCellXYZ_(mfr.numCellsXYZ()),
                                                       dimensionsXYZ_(mfr.dimensionsXYZ())
        {
            createFaceOwnerNeighbourArray() ;
            createInternalAndBoundaryFaceList() ;
            createOppositePatchFaceIndexList() ;
            //testFunctions() ;
        }

        // Utility functions 

        void createFaceOwnerNeighbourArray(){

            int l = -1;
            array<int,2> arr ;
            for(size_t i=0; i < faces_.size(); ++i){

                for(size_t j=0; j < cells_.size(); ++j){
                    if(find(cells_[j].begin(), cells_[j].end(), i)!=cells_[j].end()){
                        arr[++l] = j ;
                    }
                }

                if(l==0){
                    arr[++l] = -1 ;
                }

                faceOwnerNeighbour_.push_back(arr) ;
                l = -1 ;
            }
        }

        double calculateVolumeOfTetrahedra(array<double,3> point1, array<double,3> point2, array<double,3> point3, array<double,3> point4) const{

            double firstElem = (point2[0] - point1[0])*((point3[1] - point1[1])*(point4[2] - point1[2]) - (point4[1]-point1[1])*(point3[2] - point1[2])) ;
            double secondElem = (point2[1] - point1[1])*((point3[0] - point1[0])*(point4[2] - point1[2]) - (point4[0] - point1[0])*(point3[2] - point1[2])) ;
            double thirdElem = (point2[2] - point1[2])*((point3[0] - point1[0])*(point4[1] - point1[1]) - (point4[0] - point1[0])*(point3[1] - point1[1])) ;
            double volume = (firstElem - secondElem + thirdElem)/6 ;

            return abs(volume) ;
        }

        double calculateAreaOfTriangle(array<double,3> point1, array<double,3> point2, array<double,3> point3) const{

            array<double,3> u = {point2[0]-point1[0], point2[1]-point1[1], point2[2]-point1[2]} ;
            array<double,3> v = {point3[0]-point1[0], point3[1]-point1[1], point3[2]-point1[2]} ;

            array<double,3> n = { u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] } ;
            double area = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])/2; 
            return abs(area) ;

        }
           
        const array<double,3> returnNormalVector(array<double,3> p1, array<double,3> p2, array<double,3> p3) const{

            array<double,3> u = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]} ;
            array<double,3> v = {p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]} ;

            array<double,3> n = { u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0] } ;
            double mod = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

            return {n[0]/mod, n[1]/mod, n[2]/mod} ;
        }

        void createInternalAndBoundaryFaceList(){

            array<int,2> arr;
            for(size_t i=0; i<faces_.size(); ++i){
                arr = faceOwnerNeighbour_[i] ;
                if(arr[1] != -1){
                    internalFaces_.push_back(i) ;
                }
                else{
                    boundaryFaces_.push_back(i) ;
                }
            }

        }

        void createOppositePatchFaceIndexList(){
            const double tolerance = 1.0e-12;

            /*
                The vector is indexed using global face indices.
                Therefore, it must contain one entry for every mesh face.
                Internal and unpaired faces remain -1.
            */
            oppositeFaceIndex_.assign(numFaces_, -1);

            map<string, const vector<int>*> patchMap;

            for (const Patch& patch : boundaryPatches_)
            {
                patchMap[patch.name()] = &patch.faces();
            }

            /*
                Check that all six expected patches exist before attempting
                to dereference their face lists.
            */
            const vector<string> requiredPatches =
            {
                "east",
                "west",
                "north",
                "south",
                "top",
                "bottom"
            };

            for (const string& patchName : requiredPatches)
            {
                if (patchMap.find(patchName) == patchMap.end())
                {
                    throw runtime_error(
                        "FaceAddressedMesh3D: boundary patch '"
                        + patchName
                        + "' was not found while constructing opposite-face pairs."
                    );
                }
            }

            /*
                Helper function for pairing two opposite patches.

                normalDirection:
                    0 -> x direction
                    1 -> y direction
                    2 -> z direction
            */
            auto pairOppositePatches =
            [&](const string& firstPatchName,
                const string& secondPatchName,
                int normalDirection,
                double domainLength)
            {
                const vector<int>& firstPatchFaces =
                    *patchMap.at(firstPatchName);

                const vector<int>& secondPatchFaces =
                    *patchMap.at(secondPatchName);

                if (firstPatchFaces.size() != secondPatchFaces.size())
                {
                    throw runtime_error(
                        "FaceAddressedMesh3D: opposite patches '"
                        + firstPatchName
                        + "' and '"
                        + secondPatchName
                        + "' contain different numbers of faces."
                    );
                }

                for (int faceIndex1 : firstPatchFaces)
                {
                    const array<double, 3> faceCentre1 =
                        faceCentre(faceIndex1);

                    int matchedFaceIndex = -1;

                    for (int faceIndex2 : secondPatchFaces)
                    {
                        const array<double, 3> faceCentre2 =
                            faceCentre(faceIndex2);

                        bool transverseCoordinatesMatch = true;

                        for (int coordinate = 0;
                            coordinate < 3;
                            ++coordinate)
                        {
                            if (coordinate == normalDirection)
                            {
                                continue;
                            }

                            if (
                                abs(
                                    faceCentre1[coordinate]
                                    - faceCentre2[coordinate]
                                ) > tolerance
                            )
                            {
                                transverseCoordinatesMatch = false;
                                break;
                            }
                        }

                        const bool normalDistanceMatches =
                            abs(
                                abs(
                                    faceCentre1[normalDirection]
                                    - faceCentre2[normalDirection]
                                )
                                - domainLength
                            ) <= tolerance;

                        if (
                            transverseCoordinatesMatch
                            && normalDistanceMatches
                        )
                        {
                            if (matchedFaceIndex != -1)
                            {
                                throw runtime_error(
                                    "FaceAddressedMesh3D: more than one opposite face "
                                    "was found for face "
                                    + to_string(faceIndex1)
                                    + " on patch '"
                                    + firstPatchName
                                    + "'."
                                );
                            }

                            matchedFaceIndex = faceIndex2;
                        }
                    }

                    if (matchedFaceIndex == -1)
                    {
                        throw runtime_error(
                            "FaceAddressedMesh3D: no opposite face was found for face "
                            + to_string(faceIndex1)
                            + " on patch '"
                            + firstPatchName
                            + "'."
                        );
                    }

                    if (
                        oppositeFaceIndex_[faceIndex1] != -1
                        || oppositeFaceIndex_[matchedFaceIndex] != -1
                    )
                    {
                        throw runtime_error(
                            "FaceAddressedMesh3D: a boundary face was assigned "
                            "to more than one opposite face."
                        );
                    }

                    oppositeFaceIndex_[faceIndex1] =
                        matchedFaceIndex;

                    oppositeFaceIndex_[matchedFaceIndex] =
                        faceIndex1;
                }
            };

            // East–west faces differ only in the x direction.
            pairOppositePatches(
                "east",
                "west",
                0,
                dimensionsXYZ_[0]
            );

            // North–south faces differ only in the y direction.
            pairOppositePatches(
                "north",
                "south",
                1,
                dimensionsXYZ_[1]
            );

            // Top–bottom faces differ only in the z direction.
            pairOppositePatches(
                "top",
                "bottom",
                2,
                dimensionsXYZ_[2]
            );
        }

        // Global data

        const int& numPoints() const{
            return numPoints_ ;
        }

        const int& numFaces() const{
            return numFaces_ ;
        }

        const int& numCells() const{
            return numCells_ ;
        }

        const int& numBoundaryPatches() const{
            return numBoundaryPatches_ ;
        }

        const vector<int>& numCellsXYZ() const{
            return numCellXYZ_ ;
        }

        const vector<double>& dimensionXYZ() const{
            return dimensionsXYZ_ ;
        }

        // Geometrical data

        const array<double,3> cellCentre(int cellIndex) const{

            set<int> pointSet ;
            for(size_t i = 0; i < cells_[cellIndex].size(); ++i){
                for(int j: faces_[cells_[cellIndex][i]]){
                    pointSet.insert(j) ;
                }
            }

            int numCellPoints = pointSet.size() ;
            double xCentre=0, yCentre=0, zCentre=0; 

            for(int point: pointSet){
                xCentre += points_[point][0] ;
                yCentre += points_[point][1] ;
                zCentre += points_[point][2] ;
            }

            return {xCentre/numCellPoints, yCentre/numCellPoints, zCentre/numCellPoints} ;
        }

        const array<double,3> faceCentre(int faceIndex) const{

            int numFacePoints = faces_[faceIndex].size() ;
            double xCentre=0, yCentre=0, zCentre=0; 

            for(int point: faces_[faceIndex]){
                xCentre += points_[point][0] ;
                yCentre += points_[point][1] ;
                zCentre += points_[point][2] ;
            }

            return {xCentre/numFacePoints, yCentre/numFacePoints, zCentre/numFacePoints} ;

        }

        double cellVolume(int cellIndex) const{

            array<double,3> cellCentre_ = cellCentre(cellIndex) ;

            double cellVolume = 0 ;

            for(int faceIndex: cells_[cellIndex]){

                vector<array<double,3>> cellPoints ;
                for(int pointIndex: faces_[faceIndex]){
                    cellPoints.push_back(points_[pointIndex]) ;
                }

                cellVolume += calculateVolumeOfTetrahedra(cellCentre_, cellPoints[0], cellPoints[1], cellPoints[2]) +
                              calculateVolumeOfTetrahedra(cellCentre_, cellPoints[0], cellPoints[2], cellPoints[3]) ;
            }

            return cellVolume ;
        }

        const array<double,3> faceAreaVector(int faceIndex) const{
            const array<double,3> faceCentre_ =
                faceCentre(faceIndex);

            vector<array<double,3>> facePoints;

            for (int pointIndex : faces_[faceIndex])
            {
                facePoints.push_back(
                    points_[pointIndex]
                );
            }

            /*
                Calculate the initial normal direction from the face-point
                ordering.
            */
            array<double,3> areaVector =
                returnNormalVector(
                    facePoints[0],
                    facePoints[1],
                    facePoints[2]
                );

            /*
                Calculate the total face area.
            */
            const double area =
                calculateAreaOfTriangle(
                    faceCentre_,
                    facePoints[0],
                    facePoints[1]
                )
                +
                calculateAreaOfTriangle(
                    faceCentre_,
                    facePoints[1],
                    facePoints[2]
                )
                +
                calculateAreaOfTriangle(
                    faceCentre_,
                    facePoints[2],
                    facePoints[3]
                )
                +
                calculateAreaOfTriangle(
                    faceCentre_,
                    facePoints[3],
                    facePoints[0]
                );

            areaVector[0] *= area;
            areaVector[1] *= area;
            areaVector[2] *= area;

            /*
                Orient the area vector consistently.

                Internal face:
                    point from owner cell to neighbour cell.

                Boundary face:
                    point outward from the owner cell towards the face centre.
            */
            const array<int,2>& ownerNeighbour =
                faceOwnerNeighbour_[faceIndex];

            const array<double,3> ownerCentre =
                cellCentre(ownerNeighbour[0]);

            array<double,3> requiredDirection;

            if (ownerNeighbour[1] != -1)
            {
                const array<double,3> neighbourCentre =
                    cellCentre(ownerNeighbour[1]);

                requiredDirection =
                {
                    neighbourCentre[0] - ownerCentre[0],
                    neighbourCentre[1] - ownerCentre[1],
                    neighbourCentre[2] - ownerCentre[2]
                };
            }
            else
            {
                requiredDirection =
                {
                    faceCentre_[0] - ownerCentre[0],
                    faceCentre_[1] - ownerCentre[1],
                    faceCentre_[2] - ownerCentre[2]
                };
            }

            const double orientation =
                areaVector[0] * requiredDirection[0]
                + areaVector[1] * requiredDirection[1]
                + areaVector[2] * requiredDirection[2];

            if (orientation < 0.0)
            {
                areaVector[0] = -areaVector[0];
                areaVector[1] = -areaVector[1];
                areaVector[2] = -areaVector[2];
            }

            return areaVector;
        }

        const array<double,3> faceNormalVector(int faceIndex) const{

            array<double,3> faceVec = faceAreaVector(faceIndex) ;
            double mag = sqrt(faceVec[0]*faceVec[0] + faceVec[1]*faceVec[1] + faceVec[2]*faceVec[2]) ;

            faceVec[0] = faceVec[0]/mag ;
            faceVec[1] = faceVec[1]/mag ;
            faceVec[2] = faceVec[2]/mag ;

            return faceVec ;
        }

        double faceArea(int faceIndex) const{
            array<double,3> faceVec = faceAreaVector(faceIndex) ;
            return sqrt(faceVec[0]*faceVec[0] + faceVec[1]*faceVec[1] + faceVec[2]*faceVec[2]) ;
        }

        const vector<int>& cellFaces(int cellIndex) const{
            return cells_[cellIndex] ;
        }

        // Topological Data

        void boundaryFaceOwner(int patchIndex){
            cout << "Boundary Patch: " << patchIndex << endl ;
            const Patch p = boundaryPatches_[patchIndex] ;
            for(int faceIndex: p.faces()){
                    cout << "Face: " << faceIndex << ": Owner = " << faceOwnerNeighbour_[faceIndex][0] << endl ;
            }
        }

        const vector<int> cellNeighbours(int cellIndex) const{

            vector<int> neighbours ;
            for(int faceIndex: cells_[cellIndex]){
                if(faceOwnerNeighbour_[faceIndex][0]!=cellIndex){
                    neighbours.push_back(faceOwnerNeighbour_[faceIndex][0]) ;
                }
                else{
                    neighbours.push_back(faceOwnerNeighbour_[faceIndex][1]) ;
                }
            }

            return neighbours ;
        }

        const array<int,2> faceOwnerNeighbour(int faceIndex) const{
            return {faceOwnerNeighbour_[faceIndex][0], faceOwnerNeighbour_[faceIndex][1]} ;
        }

        const vector<int>& internalFaces() const{
            return internalFaces_ ;
        }

        const vector<int>& boundaryFaces() const{
            return boundaryFaces_ ;
        }

        const vector<Patch>& boundaryPatches() const{
            return boundaryPatches_ ;
        }

        const vector<int>& oppositeFaceIndex() const{
            return oppositeFaceIndex_ ;
        }

        const array<double,3> findMeshCentre() const{
            return {dimensionsXYZ_[0]/2, dimensionsXYZ_[1]/2, dimensionsXYZ_[2]/2} ;
        }

        // Test Function
        void testFunctions(){
            ofstream out("/Users/amitanshusahoo/Desktop/Project/P1/P1Model/transformationFiles/P1Radiation/Output/TestFileOutputMesh.txt") ;

            createFaceOwnerNeighbourArray() ;
            createInternalAndBoundaryFaceList() ;

            for(size_t i = 0; i < faceOwnerNeighbour_.size(); ++i){
                out << i << ": owner = " << faceOwnerNeighbour_[i][0] << " neighbour = " << faceOwnerNeighbour_[i][1] << endl ;  
            }

            out << endl ;
            out << "Number of Points: " << numPoints() << endl ;
            out << "Number of Cells: " << numCells() << endl ;
            out << "Number of Faces: " << numFaces() << endl ;
            out << "Number of Boundary Patches: " << numBoundaryPatches() << endl ;

            out << endl ;
            array<double,3> cellCentre_ ;
            int numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellCentre_ = cellCentre(i) ;
                out << "Cell Index: " << i << ": Cell Centre = {" << cellCentre_[0] << "," << cellCentre_[1] << "," << cellCentre_[2] << "}" << endl ;
            }

            out << endl ;
            array<double,3> faceCentre_ ;
            int numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceCentre_ = faceCentre(i) ;
                out << "Face Index: " << i << ": Face Centre = {" << faceCentre_[0] << "," << faceCentre_[1] << "," << faceCentre_[2] << "}" << endl ;
            }

            out << endl ;
            double cellVolume_ ;
            numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellVolume_ = cellVolume(i) ;
                out << "Cell Index: " << i << ": Cell Volume = " << cellVolume_ << endl ;
            }

            out << endl ;
            array<double,3> faceVec ;
            numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceVec = faceAreaVector(i) ;
                out << "Face Index: " << i << ": Face Area Vector = {" << faceVec[0] << "," << faceVec[1] << "," << faceVec[2] << "}" << endl ;
            }

            out << endl ;
            int numPatches = numBoundaryPatches() ;
            for(int i=0; i<numPatches; ++i){
                boundaryFaceOwner(i) ;
                out << endl ;
            }

            out << endl ;
            vector<int> cellNeighbours_ ;
            numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellNeighbours_ = cellNeighbours(i) ;
                out << "Cell Index: " << i << ": Cell Neighbours: " ;
                for(int var: cellNeighbours_){
                    out << var << " " ;
                }
                out << endl ;
            }

            out << endl ;
            array<int,2> faceON ;
            numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceON = faceOwnerNeighbour(i) ;
                out << "Face Index: " << i << ": Face {Owner,Neighbour} = {" << faceON[0] << "," << faceON[1] << "}" << endl ;
            }

            out << endl << endl << "Internal Faces: " ;
            for(size_t i=0; i<internalFaces_.size(); ++i){
                out << internalFaces_[i] << " " ;
            }

            out << endl << endl << "Boundary Faces: " ;
            for(size_t i=0; i<boundaryFaces_.size(); ++i){
                out << boundaryFaces_[i] << " " ;
            }

        }


} ;
