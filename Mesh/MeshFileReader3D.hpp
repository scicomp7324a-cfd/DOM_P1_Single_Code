# pragma once

# include<iostream>
# include <fstream>
# include <string>
# include <vector>
# include <array>
# include <sstream>
# include <map>
# include "Patch.hpp"

using namespace std ;

class MeshFileReader3D{

    private:

    vector<array<double,3>> points_ ;
    vector<vector<int>> faces_ ;
    vector<vector<int>> cells_ ;
    vector<Patch> boundaryPatches_ ;

    int numPoints_ ;
    int numFaces_ ;
    int numCells_ ;
    int numBoundaryPatches_ ;

    double lx, ly, lz ;
    int nx, ny, nz ;

    public:

    MeshFileReader3D(vector<string> fileNames){

        readPointsFile(fileNames[0]) ;
        readFacesFile(fileNames[1]) ;
        readCellsFile(fileNames[2]) ;
        readBoundaryPatchesFile(fileNames[3]) ;   
        readDimensionsFile(fileNames[4]) ;   
        //test() ;  

    }

    const vector<array<double,3>>& points() const{
        return points_ ;
    }

    const vector<vector<int>> faces() const{
        return faces_ ;
    }

    const vector<vector<int>> cells() const{
        return cells_ ;
    }

    const vector<Patch> boundaryPatches() const{
        return boundaryPatches_ ;
    }

    const vector<int> numCellsXYZ() const{
        return {nx, ny, nz} ;
    }

    const vector<double> dimensionsXYZ() const{
        return {lx, ly, lz} ;
    }

    int numPoints() const{
        return numPoints_ ;
    }

    int numFaces() const{
        return numFaces_ ;
    }

    int numCells() const{
        return numCells_ ;
    }

    int numBoundaryPatches() const{
        return numBoundaryPatches_ ;
    }

    vector<string> splitWhiteSpace(string s) {
        istringstream iss(s);
        vector<string> out;
        for (string tok; iss >> tok; ) out.push_back(tok);
        return out;
    }

    string stripWhiteSpace(string s){
        auto first = s.find_first_not_of(" \t\r");
        if (first == string::npos) return ""; 
        auto last  = s.find_last_not_of(" \t\r");
        return s.substr(first, last - first + 1);
    }

    void readPointsFile(string pointsFilePath){

        ifstream in(pointsFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){
                line = line.substr(line.find("(")+1, line.find(")")-1) ; 
                strVec = splitWhiteSpace(line) ;
                points_.push_back({stod(strVec[0]), stod(strVec[1]), stod(strVec[2])}) ;

            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numPoints_ = stoi(line) ;
            }
        }
    }

    void readFacesFile(string facesFilePatch){

        ifstream in(facesFilePatch);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){

                int size = stoi(line.substr(0, line.find("("))) ;
                vector<int> arr ;
                line = line.substr(line.find("(")+1, line.find(")")-1) ; 
                strVec = splitWhiteSpace(line) ;
                
                for(int i = 0; i < size; ++i){
                    arr.push_back(stoi(strVec[i])) ;
                }
                faces_.push_back(arr) ;
            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numFaces_ = stoi(line) ;
            }
        }
    }

    void readCellsFile(string cellsFilePath){

        ifstream in(cellsFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){

                int size = stoi(line.substr(0, line.find("("))) ;
                vector<int> arr ;
                line = line.substr(line.find("(")+1, line.find(")")-2) ; 
                strVec = splitWhiteSpace(line) ;

                for(int i = 0; i < size; ++i){
                    arr.push_back(stoi(strVec[i])) ;
                }
                cells_.push_back(arr) ;

            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numCells_ = stoi(line) ;
            }
        }
    }

    void readBoundaryPatchesFile(string boundaryPatchesFilePath){

        ifstream in(boundaryPatchesFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;

        getline(in, line) ;
        line = stripWhiteSpace(line) ;
        numBoundaryPatches_ = stoi(line) ;

        while(getline(in,line)){

            if(line.find("east")!=string::npos || line.find("west")!=string::npos ||
            line.find("north")!=string::npos || line.find("south")!=string::npos ||
            line.find("top")!=string::npos || line.find("bottom")!=string::npos){

                string patchType = stripWhiteSpace(line) ;
                getline(in, line) ; 
                int size = stoi(line) ;
                
                while(getline(in,line)){

                    if(line.find("(")!=string::npos){
                        continue ;
                    }
                    else if(line.find(")")!=string::npos){
                        break ;
                    }
                    else{
                        vector<int> arr ;
                        strVec = splitWhiteSpace(line) ;
                        for(int i = 0; i < size; ++i){
                            arr.push_back(stoi(strVec[i])) ;
                        }
                        boundaryPatches_.push_back(Patch(patchType, boundaryPatches_.size(), arr, size)) ;
                    }
                }
            }
        }
    }

    void readDimensionsFile(string dimensionFileName){

        ifstream in(dimensionFileName) ;
        if(!in){
            cerr << "File could not be open" << endl ;
            throw std::runtime_error("File could not be open") ;
        }

        string line ;
        vector<string> strvec ;
         

        while(getline(in, line)){

            if(line.find("Lx") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                lx = stod(strvec[1]) ;
            }

            if(line.find("Ly") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                ly = stod(strvec[1]) ;
            }

            if(line.find("Lz") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                lz = stod(strvec[1]) ;
            }

            if(line.find("Nx") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                nx = stoi(strvec[1]) ;
            }

            if(line.find("Ny") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                ny = stoi(strvec[1]) ;
            }

            if(line.find("Nz") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;
                nz = stoi(strvec[1]) ;
            }
        }
    }

    void test(){

        ofstream out("/Users/amitanshusahoo/Desktop/Project/P1/P1Model/transformationFiles/P1Radiation/Output/TestFileOutput.txt") ;

        for(size_t i=0; i<points_.size(); ++i){
            out << "Point: {" << points_[i][0] << "," << points_[i][1] << "," << points_[i][2] << "}" << endl ;
        }

        for(vector<int> var: faces_){
            out << "Face: " ;
            for(int val: var){
                out << val << " " ;
            }
            out << endl ;
        }

        for(vector<int> var: cells_){
            out << "Cell: " ;
            for(int val: var){
                out << val << " " ;
            }
            out << endl ;
        }

        for(Patch var: boundaryPatches_){
            out << "Boundary Patch: " << endl ;
            out << "Name: " << var.name() << endl ;
            out << "PatchIndex: " << var.index() << endl ;
            out << "Size: " << var.numFaces() << endl ;
            out << "Faces: " ;
            for(int f: var.faces()){
                out << f << " " ;
            }
            out << endl ;
        }
    }

} ;
