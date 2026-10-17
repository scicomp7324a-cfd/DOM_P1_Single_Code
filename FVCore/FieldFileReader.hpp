# pragma once 

# include<iostream>
# include<array>
# include<vector>
# include<fstream>
# include<sstream>
# include "FaceAddressedMesh3D.hpp"
# include "SparseLinearVector.hpp"
# include "Patch.hpp"
# include "BoundaryCondition.hpp"
# include "FixedValueBC.hpp"
# include "GradientBC.hpp"
# include "MarshakRadiationBC.hpp"
# include "FVMatrix.hpp"
# include "ScalarProperty.hpp"
# include "FaceFluxField.hpp"
# include "Field.hpp"

using namespace std ;

/* *******************************************************************
    FieldFileReader.hpp -  This class contains functions to read field 
    files and set boundary condition pointers depending upon the type
    of condition.  
******************************************************************* */

class FieldFileReader{

    private:
    const FaceAddressedMesh3D& mesh_ ;
    FVMatrix& fvMatrix_ ;
    FaceFluxField& flux_ ;
    ScalarProperty& sp_ ;
    WallProperty wp_ ;
    Field& field_ ;

    SparseLinearVector internal_ ;
    vector<unique_ptr<BoundaryCondition>> boundary_ ;
    
    
    public:
    
    // Constructs the field file reader object using the mesh, fvMatrix, flux, associated scalar properties and field files.
    FieldFileReader(const FaceAddressedMesh3D& mesh, FVMatrix& fvMatrix, FaceFluxField& flux, ScalarProperty& sp, Field& field, string filename) :mesh_(mesh),
     fvMatrix_(fvMatrix),
     flux_(flux),
     sp_(sp),
     wp_(mesh_),
     field_(field),
     internal_(mesh_.numCells())
    {
        boundary_.resize(mesh_.boundaryPatches().size()) ;
        readFixedBCFile(filename) ;
    }

    FieldFileReader(const FaceAddressedMesh3D& mesh, FVMatrix& fvMatrixG, FaceFluxField& flux, ScalarProperty& sp, WallProperty& wp, Field& fieldT, string filename):mesh_(mesh),
                                        fvMatrix_(fvMatrixG),
                                        flux_(flux),
                                        sp_(sp),
                                        wp_(wp),
                                        field_(fieldT),
                                        internal_(mesh_.numCells())               
    {
        boundary_.resize(mesh_.boundaryPatches().size()) ;
        readMarshakRadiationBCFile(filename) ;
    }


    // Non-constant access function for initial internal field values
    SparseLinearVector& internal() {
        return internal_ ;
    }

    // Non-constant access function for boundary field values
    vector<unique_ptr<BoundaryCondition>>& boundary() {
        return boundary_ ;
    }

    // Function to read fixed value or fixed gradient field file 
    void readFixedBCFile(string filename){

        ifstream in(filename);
        vector<string> strvec ;
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        while(getline(in, line)){

            if(line.find("internalField") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;

                if(strvec[1] == "uniform" ){
                    double value = stod(strvec[2]) ;
                    for(int i = 0; i < mesh_.numCells(); ++i){
                        internal_(i) = value ;
                    }
                }
                else if(strvec[1] == "nonuniform"){

                    if(mesh_.numCells() != stoi(strvec[2])){throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of cells in the mesh");}
                    getline(in,line) ;
                    for(int i = 0; i < mesh_.numCells(); ++i){
                        getline(in,line) ;
                        line = stripWhiteSpace(line) ;
                        double value = stod(line) ;
                        internal_(i) = value ;
                    }
                }
            }

            if(line.find("boundaryField") != string::npos){
                while(getline(in,line)){
                    
                    for(Patch p: mesh_.boundaryPatches()){

                        if(line.find(p.name()) != string::npos){

                            while(getline(in, line)){

                                if(line.find("fixedValue") != string::npos){

                                    getline(in, line) ;
                                    line = stripWhiteSpace(line) ;
                                    strvec = splitWhiteSpace(line) ;
                                    
                                    if(strvec[1] == "uniform" ){

                                        map<int, double> faceVal_ ;
                                        double value = stod(strvec[2]) ;
                                        for(int faceIndex: p.faces()){
                                            faceVal_.insert({faceIndex, value}) ;
                                        }
                                        boundary_[p.index()] = make_unique<FixedValueBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ;
                                        
                                    }
                                    else{

                                        if(p.numFaces() != stoi(strvec[2])) {throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of faces in the patch"); }

                                        map<int, double> faceVal_ ;
                                        getline(in, line) ;

                                        for(int faceIndex: p.faces()){
                                            getline(in,line) ;
                                            line = stripWhiteSpace(line) ;
                                            double value = stod(line) ;
                                            faceVal_.insert({faceIndex, value}) ;
                                        }

                                        getline(in,line) ;
                                        boundary_[p.index()] = make_unique<FixedValueBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ;
                                                
                                    }
                                }
                                else if(line.find("fixedGradient") != string::npos || line.find("zeroGradient") != string::npos){

                                    getline(in, line) ;
                                    line = stripWhiteSpace(line) ;
                                    strvec = splitWhiteSpace(line) ;
                                    
                                    if(strvec[1] == "uniform"){
                                        map<int, double> faceVal_ ;
                                        double value = stod(strvec[2]) ;
                                        for(int faceIndex: p.faces()){
                                            faceVal_.insert({faceIndex, value}) ;
                                        }
                                        boundary_[p.index()] = make_unique<GradientBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ; 
                                    }
                                    else{
                                        if(p.numFaces() != stoi(strvec[2])) {throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of faces in the patch"); }

                                        map<int, double> faceVal_ ;
                                        getline(in, line) ;

                                        for(int faceIndex: p.faces()){
                                            getline(in,line) ;
                                            line = stripWhiteSpace(line) ;
                                            double value = stod(line) ;
                                            faceVal_.insert({faceIndex, value}) ;
                                        }

                                        getline(in,line) ;
                                        boundary_[p.index()] = make_unique<GradientBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ;
                                    }
                                }
                                else if(line.find("}") != string::npos){
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        
        }
    }

    // Function to read marshak radiation boundary condition field file
    void readMarshakRadiationBCFile(string filename){
        
        ifstream in(filename);
        vector<string> strvec ;

        if (!in) {
            cerr << "Failed to open file" << endl ;
            throw std::runtime_error("Failed to open file") ;
        }

        string line ;
        while(getline(in, line)){

            if(line.find("internalField") != string::npos){
                line = stripWhiteSpace(line) ;
                strvec = splitWhiteSpace(line) ;

                if(strvec[1] == "uniform" ){
                    double value = stod(strvec[2]) ;
                    for(int i = 0; i < mesh_.numCells(); ++i){
                        internal_(i) = value ;
                    }
                }
                else if(strvec[1] == "nonuniform"){

                    if(mesh_.numCells() != stoi(strvec[2])){throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of cells in the mesh");}
                    getline(in,line) ;
                    for(int i = 0; i < mesh_.numCells(); ++i){
                        getline(in,line) ;
                        line = stripWhiteSpace(line) ;
                        double value = stod(line) ;
                        internal_(i) = value ;
                    }
                }
            }

            if(line.find("boundaryField") != string::npos){
                while(getline(in,line)){
                    
                    for(Patch p: mesh_.boundaryPatches()){

                        if(line.find(p.name()) != string::npos){

                            while(getline(in, line)){

                                if(line.find("marshakRadiation") != string::npos){

                                    getline(in, line) ;
                                    line = stripWhiteSpace(line) ;
                                    strvec = splitWhiteSpace(line) ;
                                    
                                    if(strvec[1] == "uniform" ){

                                        map<int, double> faceVal_ ;
                                        double value = stod(strvec[2]) ;
                                        for(int faceIndex: p.faces()){
                                            faceVal_.insert({faceIndex, value}) ;
                                        }
                                        map<int,double> marshakTemperatureField_ = field_.boundary()[p.index()]->boundaryValues() ;
                                        boundary_[p.index()] = make_unique<MarshakRadiationBC>(fvMatrix_, mesh_, flux_, sp_, wp_(p.index()), faceVal_, marshakTemperatureField_) ;
                                        
                                    }
                                    else{

                                        if(p.numFaces() != stoi(strvec[2])) {throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of faces in the patch"); }

                                        map<int, double> faceVal_ ;
                                        getline(in, line) ;

                                        for(int faceIndex: p.faces()){
                                            getline(in,line) ;
                                            line = stripWhiteSpace(line) ;
                                            double value = stod(line) ;
                                            faceVal_.insert({faceIndex, value}) ;
                                        }
                                        map<int,double> marshakTemperatureField_ = field_.boundary()[p.index()]->boundaryValues() ;
                                        getline(in,line) ;
                                        boundary_[p.index()] = make_unique<MarshakRadiationBC>(fvMatrix_, mesh_, flux_, sp_, wp_(p.index()), faceVal_, marshakTemperatureField_) ;
                                                
                                    }
                                }
                                else if(line.find("fixedGradient") != string::npos || line.find("zeroGradient") != string::npos){

                                                                        getline(in, line) ;
                                    line = stripWhiteSpace(line) ;
                                    strvec = splitWhiteSpace(line) ;
                                    
                                    if(strvec[1] == "uniform"){
                                        map<int, double> faceVal_ ;
                                        double value = stod(strvec[2]) ;
                                        for(int faceIndex: p.faces()){
                                            faceVal_.insert({faceIndex, value}) ;
                                        }
                                        boundary_[p.index()] = make_unique<GradientBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ; 
                                    }
                                    else{
                                        if(p.numFaces() != stoi(strvec[2])) {throw std::runtime_error("Number of values in the non uniform value list is not equal to the number of faces in the patch"); }

                                        map<int, double> faceVal_ ;
                                        getline(in, line) ;

                                        for(int faceIndex: p.faces()){
                                            getline(in,line) ;
                                            line = stripWhiteSpace(line) ;
                                            double value = stod(line) ;
                                            faceVal_.insert({faceIndex, value}) ;
                                        }

                                        getline(in,line) ;
                                        boundary_[p.index()] = make_unique<GradientBC>(fvMatrix_, mesh_, flux_, sp_, faceVal_) ;
                                    }
                                }
                                else if(line.find("}") != string::npos){
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Function to split a string at white spaces
    vector<string> splitWhiteSpace(string s) {
        istringstream iss(s);
        vector<string> out;
        for (string tok; iss >> tok; ) out.push_back(tok);
        return out;
    }

    // Funciton to strip white spaces at the start and end of a string
    string stripWhiteSpace(string s){
        auto first = s.find_first_not_of(" \t\r");
        if (first == string::npos) return ""; 
        auto last  = s.find_last_not_of(" \t\r");
        return s.substr(first, last - first + 1);
    }

} ;
