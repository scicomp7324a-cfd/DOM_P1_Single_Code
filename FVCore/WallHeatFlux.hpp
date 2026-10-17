# pragma once

# include <iostream>
# include<vector>
# include<array>
# include<cmath>

# include "FaceAddressedMesh3D.hpp"
# include "Patch.hpp"
# include "Field.hpp"
# include "FVMatrix.hpp"
# include "ScalarProperty.hpp"

using namespace std;

struct qf{
    int faceIndex_ ;
    double qf_ ;
} ;

class WallHeatFlux{

    private:

        const FaceAddressedMesh3D& mesh_ ;
        FVMatrix& fvMatrixG_ ;
        FVMatrix& fvMatrixT_ ;
        Field& fieldG_ ;
        Field& fieldT_ ;
        const ScalarProperty& gamma_ ;

        map<string,vector<qf>> qw_ ;
        map<string,vector<qf>> Qw_ ;

        map<string,double> avgqw_ ;
        map<string,double> totalQw_ ;

        // Analytical west-east opposite-wall exchange for infinite parallel plates.
        // Sign convention: positive = heat leaving that wall and entering the domain;
        // negative = heat entering/absorbed by that wall from the domain.
        vector<double> radiativeSourceTerm_ ;

        map<string,double> oppositeWallHeatFlux_ ;
        map<string,double> oppositeWallHeatTransfer_ ;

        map<string,double> dimensionlessqw_ ;
        map<string,double> dimensionlessqw2_ ;

        map<string,double> avgWallTemperature_ ;
        map<string,double> avgNearWallTemperature_ ;

        map<string,double> wallTemperatureSlip_ ;
        map<string,double> dimensionlessWallTemperatureSlip_ ;

        map<double,double> centerlineTemperatureInX_ ;
        map<double,double> centerlineGInX_ ;
        map<double,double> centerlineFluxInX_ ;
        map<double,double> dimensionlessCenterlineTemperatureInX_ ;
        map<double,double> dimensionlessCenterlineTemperatureInX2_ ;

        map<double,double> centerlineTemperatureInZ_;
        map<double,double> centerlineGInZ_;
        map<double,double> centerlineFluxInZ_;
        map<double,double> dimensionlessCenterlineTemperatureInZ_;
        map<double,double> centerlineRadiativeSourceTermInX_ ;

        double sigma = 5.670374419e-8  ;
        
    public: 

        WallHeatFlux(const FaceAddressedMesh3D& mesh, FVMatrix& fvMatrixG, FVMatrix& fvMatrixT, Field& fieldG, Field& fieldT, ScalarProperty& gamma)
                                :mesh_(mesh),
                                 fvMatrixG_(fvMatrixG),
                                 fvMatrixT_(fvMatrixT),
                                 fieldG_(fieldG),
                                 fieldT_(fieldT),
                                 gamma_(gamma)
        {}

        map<string,vector<qf>>& qw(){ return qw_ ;}
        const map<string,vector<qf>>& qw() const{ return qw_ ;}

        map<string,double>& avgqw(){ return avgqw_ ;}
        const map<string, double>& avgqw() const{ return avgqw_ ;}

        map<string,double>& totalQw(){ return totalQw_ ;}
        const map<string, double>& totalQw() const{ return totalQw_ ;}

        vector<double>& radiativeSourceTerm(){ return radiativeSourceTerm_ ;}
        const vector<double>& radiativeSourceTerm() const { return radiativeSourceTerm_ ;}

        map<string,double>& oppositeWallHeatFlux(){ return oppositeWallHeatFlux_ ;}
        const map<string,double>& oppositeWallHeatFlux() const{ return oppositeWallHeatFlux_ ;}

        map<string,double>& oppositeWallHeatTransfer(){ return oppositeWallHeatTransfer_ ;}
        const map<string,double>& oppositeWallHeatTransfer() const{ return oppositeWallHeatTransfer_ ;}

        map<string,double>& dimensionlessQw(){ return dimensionlessqw_ ;}
        const map<string,double>& dimensionlessQw() const{ return dimensionlessqw_ ;}

        map<string,double>& dimensionlessQw2(){ return dimensionlessqw2_ ;}
        const map<string,double>& dimensionlessQw2() const{ return dimensionlessqw2_ ;}

        map<string,double>& avgWallTemperature(){ return avgWallTemperature_ ;}
        const map<string, double>& avgWallTemperature() const{ return avgWallTemperature_ ;}

        map<string,double>& avgNearWallTemperature(){ return avgNearWallTemperature_ ;}
        const map<string, double>& avgNearWallTemperature() const{ return avgNearWallTemperature_ ;}

        map<string,double>& wallTemperatureSlip(){ return wallTemperatureSlip_ ;}
        const map<string, double>& wallTemperatureSlip() const{ return wallTemperatureSlip_ ;}

        map<string,double>& dimensionlessWallTemperatureSlip(){ return dimensionlessWallTemperatureSlip_ ;}
        const map<string, double>& dimensionlessWallTemperatureSlip() const{ return dimensionlessWallTemperatureSlip_ ;}
        
        map<double,double>& centerlineTemperature(){return centerlineTemperatureInX_ ;}
        const map<double,double>& centerlineTemperature() const{return centerlineTemperatureInX_ ;} 

        map<double,double>& centerlineG(){return centerlineGInX_ ;}
        const map<double,double>& centerlineG() const{return centerlineGInX_ ;} 

        map<double,double>& centerlineFlux(){return centerlineFluxInX_ ;}
        const map<double,double>& centerlineFlux() const{return centerlineFluxInX_ ;}

        map<double,double>& dimensionlessCenterlineTemperature(){ return dimensionlessCenterlineTemperatureInX_ ;}
        const map<double,double>& dimensionlessCenterlineTemperature() const{ return dimensionlessCenterlineTemperatureInX_ ;}

        map<double,double>& dimensionlessCenterlineTemperature2(){ return dimensionlessCenterlineTemperatureInX2_ ;}
        const map<double,double>& dimensionlessCenterlineTemperature2() const{ return dimensionlessCenterlineTemperatureInX2_ ;}

        map<double,double>& centerlineTemperatureInZ(){ return centerlineTemperatureInZ_;}
        const map<double,double>& centerlineTemperatureInZ() const { return centerlineTemperatureInZ_;}

        map<double,double>& centerlineGInZ(){ return centerlineGInZ_; }
        const map<double,double>& centerlineGInZ() const { return centerlineGInZ_; }

        map<double,double>& centerlineFluxInZ(){ return centerlineFluxInZ_; }
        const map<double,double>& centerlineFluxInZ() const { return centerlineFluxInZ_; }

        map<double,double>& dimensionlessCenterlineTemperatureInZ() { return dimensionlessCenterlineTemperatureInZ_; }
        const map<double,double>& dimensionlessCenterlineTemperatureInZ() const { return dimensionlessCenterlineTemperatureInZ_; }

        map<double,double>& centerlineRadiativeSourceTerm(){ return centerlineRadiativeSourceTermInX_ ;}
        const map<double,double>& centerlineRadiativeSourceTerm() const{ return centerlineRadiativeSourceTermInX_ ;}


        void calculateWallHeatFlux(){

            for(Patch p: mesh_.boundaryPatches()){
                qw_[p.name()] = {} ;
                Qw_[p.name()] = {} ;
                avgqw_[p.name()] = 0 ;
                totalQw_[p.name()] = 0 ;
                double totalArea = 0 ;

                map<int,double> faceVal_ =
                    fieldG_.boundary()[p.index()]->boundaryValues() ;

                for(const auto& [faceIndex,val]: faceVal_){

                    int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                    const double gamma = gamma_(cellOwner) ;
                    const array<double,3> faceArea =
                        mesh_.faceAreaVector(faceIndex) ;

                    double area = sqrt(
                        faceArea[0]*faceArea[0]
                        + faceArea[1]*faceArea[1]
                        + faceArea[2]*faceArea[2]
                    ) ;

                    double q = 0 ;

                    if(fieldG_.boundary()[p.index()]->typeName() == "fixedgradient"){
                        q = gamma*val ;
                    }
                    else{
                        const double Gp = fvMatrixG_.sys().X()(cellOwner) ;
                        const double Gf = val ;
                        const double Dpf = distanceFaceToCell(faceIndex) ;

                        q = gamma*((Gf-Gp)/Dpf) ;
                    }

                    qw_[p.name()].push_back({faceIndex,q}) ;
                    Qw_[p.name()].push_back({faceIndex,q*area}) ;

                    avgqw_[p.name()] += q*area ;
                    totalQw_[p.name()] += q*area ;
                    totalArea += area ;
                }

                avgqw_[p.name()] = avgqw_[p.name()]/totalArea ;
            }
        }

        void calculateDimensionlessQw(){

            vector<string> patchNames ;
            for(Patch p: mesh_.boundaryPatches()){
                patchNames.push_back(p.name()) ;
            }

            for(int i=0 ; i< mesh_.numBoundaryPatches(); ++i){
                if(i%2 == 0){
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i+1]],4) ;
                    if(abs(denom) > 1e-14){dimensionlessqw_[patchNames[i]] = abs(avgqw_[patchNames[i]]/(sigma*(denom))) ;}
                    else{dimensionlessqw_[patchNames[i]] = 0;}
                }
                else{
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i-1]],4) ;
                    if(abs(denom) > 1e-14){dimensionlessqw_[patchNames[i]] = abs(avgqw_[patchNames[i]]/(sigma*(denom))) ;}
                    else{dimensionlessqw_[patchNames[i]] = 0;}
                }
            }
        }

        void calculateDimensionlessQw2(){

            dimensionlessqw2_.clear() ;

            vector<string> patchNames ;
            for(Patch p: mesh_.boundaryPatches()){
                patchNames.push_back(p.name()) ;
            }

            for(int i=0 ; i<mesh_.numBoundaryPatches(); ++i){
                double denom = powl(avgWallTemperature_[patchNames[i]],4) ;
                if(abs(denom) > 1e-14){
                    dimensionlessqw2_[patchNames[i]] = abs(avgqw_[patchNames[i]]/(sigma*denom)) ;
                }
                else{
                    dimensionlessqw2_[patchNames[i]] = 0 ;
                }
            }
        }
       

        void calculateAvgTempartures(){

            for(Patch p: mesh_.boundaryPatches()){
                avgWallTemperature_[p.name()] = 0 ;
                avgNearWallTemperature_[p.name()] = 0 ;
                wallTemperatureSlip_[p.name()] = 0 ;

                double totalTemperature = 0 ;
                double totalNearWallTemperature = 0 ;
                double totalArea = 0 ;

                map<int,double> faceVal_ = fieldT_.boundary()[p.index()]->boundaryValues() ;
                for(const auto& [faceIndex, val]: faceVal_){
                    int cellOwner = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                    const array<double,3> faceArea = mesh_.faceAreaVector(faceIndex) ;
                    double area = sqrt(faceArea[0]*faceArea[0] + faceArea[1]*faceArea[1] + faceArea[2]*faceArea[2]) ;
                    totalTemperature += val*area ;
                    totalNearWallTemperature += fvMatrixT_.sys().X()(cellOwner)*area ;
                    totalArea += area ;
                }
                totalTemperature = totalTemperature/totalArea ;
                totalNearWallTemperature = totalNearWallTemperature/totalArea ;

                avgWallTemperature_[p.name()] = totalTemperature ;
                avgNearWallTemperature_[p.name()] = totalNearWallTemperature ;
                wallTemperatureSlip_[p.name()] = abs(avgWallTemperature_[p.name()] - avgNearWallTemperature_[p.name()]) ;
            }
        }

        void calculateDimensionlessSlipTemperature(){

            vector<string> patchNames ;
            for(Patch p: mesh_.boundaryPatches()){
                patchNames.push_back(p.name()) ;
            }

            for(int i=0 ; i< mesh_.numBoundaryPatches(); ++i){
                if(i%2 == 0){
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i+1]],4) ;
                    if(abs(denom) > 1e-14){ dimensionlessWallTemperatureSlip_[patchNames[i]] = abs((powl(avgWallTemperature_[patchNames[i]],4) - powl(avgNearWallTemperature_[patchNames[i]],4))/(denom)) ;}
                    else{ dimensionlessWallTemperatureSlip_[patchNames[i]] = 0 ;}
                }
                else{
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i-1]],4) ;
                    if(abs(denom) > 1e-14){ dimensionlessWallTemperatureSlip_[patchNames[i]] = abs((powl(avgWallTemperature_[patchNames[i]],4) - powl(avgNearWallTemperature_[patchNames[i]],4))/(denom)) ;}
                    else{ dimensionlessWallTemperatureSlip_[patchNames[i]] = 0 ;}                }
            }
        }

        void calculateCenterlineTemperatureAndG(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int jc = numY/2 ;
            int kc = numZ/2 ;

            for(int i=0; i<mesh_.numCellsXYZ()[0]; ++i){
                int cellIndex = kc*numX*numY + jc*numX + i ;
                double x = mesh_.cellCentre(cellIndex)[0];
                double valueT = fvMatrixT_.sys().X()(cellIndex) ;
                double valueG = fvMatrixG_.sys().X()(cellIndex) ;
                centerlineTemperatureInX_[x] = valueT ;
                centerlineGInX_[x] = valueG ;
            }
        }


        void calculateCenterlineFlux(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int jc = numY/2 ;
            int kc = numZ/2 ;

            centerlineFluxInX_.clear() ;

            if(numX == 1){
                int cellIndex = kc*numX*numY + jc*numX ;
                double x = mesh_.cellCentre(cellIndex)[0] ;
                centerlineFluxInX_[x] = 0.0 ;
                return ;
            }

            for(int i=0; i<numX; ++i){

                int cellIndex = kc*numX*numY + jc*numX + i ;
                double x = mesh_.cellCentre(cellIndex)[0] ;
                double dGdx = 0.0 ;

                if(i == 0){
                    int cellIndexRight = kc*numX*numY + jc*numX + (i+1) ;
                    double xRight = mesh_.cellCentre(cellIndexRight)[0] ;
                    double GRight = fvMatrixG_.sys().X()(cellIndexRight) ;
                    double GCell = fvMatrixG_.sys().X()(cellIndex) ;
                    dGdx = (GRight-GCell)/(xRight-x) ;
                }
                else if(i == numX-1){
                    int cellIndexLeft = kc*numX*numY + jc*numX + (i-1) ;
                    double xLeft = mesh_.cellCentre(cellIndexLeft)[0] ;
                    double GLeft = fvMatrixG_.sys().X()(cellIndexLeft) ;
                    double GCell = fvMatrixG_.sys().X()(cellIndex) ;
                    dGdx = (GCell-GLeft)/(x-xLeft) ;
                }
                else{
                    int cellIndexLeft = kc*numX*numY + jc*numX + (i-1) ;
                    int cellIndexRight = kc*numX*numY + jc*numX + (i+1) ;
                    double xLeft = mesh_.cellCentre(cellIndexLeft)[0] ;
                    double xRight = mesh_.cellCentre(cellIndexRight)[0] ;
                    double GLeft = fvMatrixG_.sys().X()(cellIndexLeft) ;
                    double GRight = fvMatrixG_.sys().X()(cellIndexRight) ;
                    dGdx = (GRight-GLeft)/(xRight-xLeft) ;
                }

                centerlineFluxInX_[x] = -gamma_(cellIndex)*dGdx ;
            }
        }

        void calculateDimensionlessCenterlineTemperature(){

            for(const auto& [key,val]: centerlineTemperatureInX_){
                dimensionlessCenterlineTemperatureInX_[key] = abs((powl(val,4) - powl(avgWallTemperature_["east"],4))/(powl(avgWallTemperature_["west"],4) - powl(avgWallTemperature_["east"],4))) ;
            }

        }

        void calculateDimensionlessCenterlineTemperature2(){
            for(const auto& [key,val]: centerlineTemperatureInX_){
                dimensionlessCenterlineTemperatureInX2_[key] = abs((powl(val,4))/(powl(avgWallTemperature_["east"],4))) ;
            }
        }

        void calculateCenterlineTemperatureAndGInZ(){

            int numX = mesh_.numCellsXYZ()[0];
            int numY = mesh_.numCellsXYZ()[1];
            int numZ = mesh_.numCellsXYZ()[2];

            int ic = numX/2;
            int jc = numY/2;

            centerlineTemperatureInZ_.clear();
            centerlineGInZ_.clear();

            for(int k=0; k<numZ; ++k){

                int cellIndex = k*numX*numY + jc*numX + ic;

                double z = mesh_.cellCentre(cellIndex)[2];

                double valueT = fvMatrixT_.sys().X()(cellIndex);
                double valueG = fvMatrixG_.sys().X()(cellIndex);

                centerlineTemperatureInZ_[z] = valueT;
                centerlineGInZ_[z] = valueG;
            }
        }


        void calculateCenterlineFluxInZ(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int ic = numX/2 ;
            int jc = numY/2 ;

            centerlineFluxInZ_.clear() ;

            if(numZ == 1){
                int cellIndex = jc*numX + ic ;
                double z = mesh_.cellCentre(cellIndex)[2] ;
                centerlineFluxInZ_[z] = 0.0 ;
                return ;
            }

            for(int k=0; k<numZ; ++k){

                int cellIndex = k*numX*numY + jc*numX + ic ;
                double z = mesh_.cellCentre(cellIndex)[2] ;
                double dGdz = 0.0 ;

                if(k == 0){
                    int cellIndexTop = (k+1)*numX*numY + jc*numX + ic ;
                    double zTop = mesh_.cellCentre(cellIndexTop)[2] ;
                    double GTop = fvMatrixG_.sys().X()(cellIndexTop) ;
                    double GCell = fvMatrixG_.sys().X()(cellIndex) ;
                    dGdz = (GTop-GCell)/(zTop-z) ;
                }
                else if(k == numZ-1){
                    int cellIndexBottom = (k-1)*numX*numY + jc*numX + ic ;
                    double zBottom = mesh_.cellCentre(cellIndexBottom)[2] ;
                    double GBottom = fvMatrixG_.sys().X()(cellIndexBottom) ;
                    double GCell = fvMatrixG_.sys().X()(cellIndex) ;
                    dGdz = (GCell-GBottom)/(z-zBottom) ;
                }
                else{
                    int cellIndexBottom = (k-1)*numX*numY + jc*numX + ic ;
                    int cellIndexTop = (k+1)*numX*numY + jc*numX + ic ;
                    double zBottom = mesh_.cellCentre(cellIndexBottom)[2] ;
                    double zTop = mesh_.cellCentre(cellIndexTop)[2] ;
                    double GBottom = fvMatrixG_.sys().X()(cellIndexBottom) ;
                    double GTop = fvMatrixG_.sys().X()(cellIndexTop) ;
                    dGdz = (GTop-GBottom)/(zTop-zBottom) ;
                }

                centerlineFluxInZ_[z] = -gamma_(cellIndex)*dGdz ;
            }
        }

        void calculateDimensionlessCenterlineTemperatureInZ(){

            dimensionlessCenterlineTemperatureInZ_.clear();

            const double Tw = avgWallTemperature_.at("bottom");

            if(std::abs(Tw) < 1e-14){
                throw std::runtime_error(
                    "Cannot calculate dimensionless centerline temperature in Z: avgWallTemperature_[\"bottom\"] is zero."
                );
            }

            for(const auto& [z, T] : centerlineTemperatureInZ_){
                dimensionlessCenterlineTemperatureInZ_[z] = T/Tw;
            }
        }

        void calculateCenterlineRadiativeSourceTerm(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int jc = numY/2 ;
            int kc = numZ/2 ;

            centerlineRadiativeSourceTermInX_.clear() ;

            if(radiativeSourceTerm_.size() != static_cast<size_t>(mesh_.numCells())){
                throw runtime_error("Radiative source term must be calculated before extracting the centerline.") ;
            }

            for(int i=0; i<numX; ++i){
                int cellIndex = kc*numX*numY + jc*numX + i ;
                double x = mesh_.cellCentre(cellIndex)[0] ;
                centerlineRadiativeSourceTermInX_[x] = radiativeSourceTerm_[cellIndex] ;
            }
        }

        double boundaryPatchArea(const string& patchName) const{

            double totalArea = 0.0 ;

            for(Patch p: mesh_.boundaryPatches()){
                if(p.name() == patchName){
                    map<int,double> faceVal_ = fieldG_.boundary()[p.index()]->boundaryValues() ;
                    for(const auto& [faceIndex, val]: faceVal_){
                        const array<double,3> faceArea = mesh_.faceAreaVector(faceIndex) ;
                        const double area = sqrt(faceArea[0]*faceArea[0] + faceArea[1]*faceArea[1] + faceArea[2]*faceArea[2]) ;
                        totalArea += area ;
                    }
                    return totalArea ;
                }
            }

            return 0.0 ;
        }

        void calculateWestEastTransferForParallelPlates(){

            oppositeWallHeatFlux_.clear();
            oppositeWallHeatTransfer_.clear();

            const double qWest = avgqw_.at("west");
            const double qEast = avgqw_.at("east");

            const double QWest = totalQw_.at("west");
            const double QEast = totalQw_.at("east");

            const double qWestToEast = 0.5*(qWest - qEast);
            const double QWestToEast = 0.5*(QWest - QEast);

            oppositeWallHeatFlux_["west"] = qWestToEast;
            oppositeWallHeatFlux_["east"] = -qWestToEast;

            oppositeWallHeatTransfer_["west"] = QWestToEast;
            oppositeWallHeatTransfer_["east"] = -QWestToEast;
        }

        map<string,double> parallelPlateHeatFluxBlackbodyAnalyticalValue(double Twest, double Teast){

            map<string,double> G ;
            double val = sigma*(powl(Twest,4) - powl(Teast,4)) ;
            
            G["west"] = val ;
            G["east"] = -val ;

            return G ; 
        }

        void calculateRadiativeSourceTerm(ScalarProperty& a){
            radiativeSourceTerm_.assign(mesh_.numCells(), 0.0) ;
            for(int i = 0; i < mesh_.numCells(); ++i){
                radiativeSourceTerm_[i] = a(i)*(4*sigma*pow(fvMatrixT_.sys().X()(i),4) - fvMatrixG_.sys().X()(i)) ;
            }
        }
        
        // Function to calculate the face to cell distance inside a control volume
        double distanceFaceToCell(const int& faceIndex){
            array<double,3> faceAreaVec = mesh_.faceAreaVector(faceIndex) ;
            double sfMag = sqrt(faceAreaVec[0]*faceAreaVec[0] + faceAreaVec[1]*faceAreaVec[1] + faceAreaVec[2]*faceAreaVec[2]) ;
            faceAreaVec = {faceAreaVec[0]/sfMag, faceAreaVec[1]/sfMag, faceAreaVec[2]/sfMag};
            array<double,3> p1 = mesh_.faceCentre(faceIndex) ;
            array<double,3> p2 = mesh_.cellCentre(mesh_.faceOwnerNeighbour(faceIndex)[0]) ;
            double d = (p1[0]-p2[0])*faceAreaVec[0] + (p1[1]-p2[1])*faceAreaVec[1] + (p1[2]-p2[2])*faceAreaVec[2] ;
            return abs(d);
        }

    
    
    } ;

