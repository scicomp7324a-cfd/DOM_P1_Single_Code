# pragma once

# include <iostream>
# include <vector>
# include <array>
# include <cmath>
# include <map>
# include <string>
# include <algorithm>
# include <memory>
# include <stdexcept>

# include "FaceAddressedMesh3D.hpp"
# include "DOMatrix.hpp"
# include "Patch.hpp"

struct DomInternalHeatStat{
    int cellIndex = -1 ;
    double x = 0, y = 0, z = 0 ;
};

struct DomWallElementHeatStat{
    int faceIndex = -1 ;
    double x = 0, y = 0, z = 0 ;
    double val = 0;
};

class DiscreteOrdinateWallHeatFlux{

    private:

        const FaceAddressedMesh3D& mesh_ ;
        const vector<unique_ptr<DOMatrix>>& intensity_ ;

        FVMatrix& fvMatrixT_ ;
        Field& fieldT_ ;

        SparseLinearVector& gInternal_ ;
        SparseLinearVector& gFace_ ;

        vector<DomInternalHeatStat> internalHeatFlux_ ;
        vector<DomWallElementHeatStat> wallHeatFlux_ ;

        vector<double> radiativeSourceTerm_ ;
        vector<double> wallNormalHeatFlux_ ;

        vector<double> internalIrradiation_ ;
        vector<double> wallIrradiation_ ;

        map<string,double> averageWallNormalHeatFlux_ ;
        map<string,double> averageWallIrradiation_ ;
        map<string,double> totalWallHeatTransfer_ ;

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

        map<string,double> oppositeWallHeatFlux_ ;
        map<string,double> oppositeWallHeatTransfer_ ;

        map<double,double> centerlineRadiativeSourceTermInX_ ;

        map<string,double> dimensionlessqw_ ;

        map<string,double> dimensionlessqw2_ ; // the denominator in this case is the wall temperature only

        double sigma = 5.670374419e-8  ;
        double pi_ = 3.1415926535897932 ;

    public:

        DiscreteOrdinateWallHeatFlux(const FaceAddressedMesh3D& mesh, const vector<unique_ptr<DOMatrix>>& intensity, FVMatrix& fvMatrixT, Field& fieldT, SparseLinearVector& gInternal, SparseLinearVector& gFace):mesh_(mesh),
                                                                                                         intensity_(intensity),
                                                                                                         fvMatrixT_(fvMatrixT),
                                                                                                         fieldT_(fieldT),
                                                                                                         gInternal_(gInternal),
                                                                                                         gFace_(gFace)
        {}

        vector<DomInternalHeatStat>& internalHeatFlux(){ return internalHeatFlux_ ; }
        const vector<DomInternalHeatStat>& internalHeatFlux() const{ return internalHeatFlux_ ; }

        vector<DomWallElementHeatStat>& wallHeatFlux(){ return wallHeatFlux_ ;}
        const vector<DomWallElementHeatStat>& wallHeatFlux() const{ return wallHeatFlux_ ;}

        vector<double>& wallNormalHeatFlux(){ return wallNormalHeatFlux_ ;}
        const vector<double>& wallNormalHeatFlux() const{ return wallNormalHeatFlux_ ;}

        vector<double>& internalIrradiation(){ return internalIrradiation_ ;}
        const vector<double>& internalIrradiation() const { return internalIrradiation_ ;}

        vector<double>& wallIrradiation(){ return wallIrradiation_ ;}
        const vector<double>& wallIrradiation() const{ return wallIrradiation_ ;}

        vector<double>& radiativeSourceTerm(){ return radiativeSourceTerm_ ;}
        const vector<double>& radiativeSourceTerm() const{ return radiativeSourceTerm_ ;}

        map<string,double>& averageWallNormalHeatFlux() { return averageWallNormalHeatFlux_ ;}
        const map<string,double>& averageWallNormalHeatFlux() const { return averageWallNormalHeatFlux_ ;}

        map<string,double>& totalWallHeatTransfer() { return totalWallHeatTransfer_ ;}
        const map<string,double>& totalWallHeatTransfer() const { return totalWallHeatTransfer_ ;}

        map<string,double>& averageWallIrradiation(){ return averageWallIrradiation_ ;}
        const map<string,double>& averageWallIrradiation() const { return averageWallIrradiation_ ;}

        map<string,double>& avgWallTemperature(){ return avgWallTemperature_ ;}
        const map<string, double>& avgWallTemperature() const{ return avgWallTemperature_ ;}

        map<string,double>& avgNearWallTemperature(){ return avgNearWallTemperature_ ;}
        const map<string, double>& avgNearWallTemperature() const{ return avgNearWallTemperature_ ;}

        map<string,double>& wallTemperatureSlip(){ return wallTemperatureSlip_ ;}
        const map<string, double>& wallTemperatureSlip() const{ return wallTemperatureSlip_ ;}

        map<string,double>& dimensionlessWallTemperatureSlip(){ return dimensionlessWallTemperatureSlip_ ;}
        const map<string, double>& dimensionlessWallTemperatureSlip() const{ return dimensionlessWallTemperatureSlip_ ;}
        
        map<double,double>& centerlineTemperatureInX(){return centerlineTemperatureInX_ ;}
        const map<double,double>& centerlineTemperatureInX() const{return centerlineTemperatureInX_ ;} 

        map<double,double>& centerlineGInX(){return centerlineGInX_ ;}
        const map<double,double>& centerlineGInX() const{return centerlineGInX_ ;} 

        map<double,double>& centerlineFluxInX(){return centerlineFluxInX_ ;}
        const map<double,double>& centerlineFluxInX() const{return centerlineFluxInX_ ;} 

        map<double,double>& dimensionlessCenterlineTemperatureInX(){ return dimensionlessCenterlineTemperatureInX_ ;}
        const map<double,double>& dimensionlessCenterlineTemperatureInX() const{ return dimensionlessCenterlineTemperatureInX_ ;}

        map<double,double>& dimensionlessCenterlineTemperatureInX2(){ return dimensionlessCenterlineTemperatureInX2_ ;}
        const map<double,double>& dimensionlessCenterlineTemperatureInX2() const{ return dimensionlessCenterlineTemperatureInX2_ ;}

        map<double,double>& centerlineTemperatureInZ(){ return centerlineTemperatureInZ_;}
        const map<double,double>& centerlineTemperatureInZ() const { return centerlineTemperatureInZ_;}

        map<double,double>& centerlineGInZ(){ return centerlineGInZ_; }
        const map<double,double>& centerlineGInZ() const { return centerlineGInZ_; }

        map<double,double>& centerlineFluxInZ(){ return centerlineFluxInZ_; }
        const map<double,double>& centerlineFluxInZ() const { return centerlineFluxInZ_; }

        map<double,double>& dimensionlessCenterlineTemperatureInZ() { return dimensionlessCenterlineTemperatureInZ_; }
        const map<double,double>& dimensionlessCenterlineTemperatureInZ() const { return dimensionlessCenterlineTemperatureInZ_; }

        map<string,double>& oppositeWallHeatFlux(){ return oppositeWallHeatFlux_ ;}
        const map<string,double>& oppositeWallHeatFlux() const{ return oppositeWallHeatFlux_ ;}

        map<string,double>& oppositeWallHeatTransfer(){ return oppositeWallHeatTransfer_ ;}
        const map<string,double>& oppositeWallHeatTransfer() const{ return oppositeWallHeatTransfer_ ;}

        map<string,double>& dimensionlessQw(){ return dimensionlessqw_ ;}
        const map<string,double>& dimensionlessQw() const{ return dimensionlessqw_ ;}

        map<string,double>& dimensionlessQw2(){ return dimensionlessqw2_ ; }
        const map<string,double>& dimensionlessQw2() const{ return dimensionlessqw2_ ;}

        map<double,double>& centerlineRadiativeSourceTermInX(){ return centerlineRadiativeSourceTermInX_ ;}
        const map<double,double>& centerlineRadiativeSourceTermInX() const{ return centerlineRadiativeSourceTermInX_ ;}


        void calculateHeatFlux(){

            internalHeatFlux_.assign(mesh_.numCells(), DomInternalHeatStat{}) ;
            wallHeatFlux_.assign(mesh_.numFaces(), DomWallElementHeatStat{}) ;
            wallNormalHeatFlux_.assign(mesh_.numFaces(), 0.0) ;

            for(int i=0; i<mesh_.numCells(); ++i){
                internalHeatFlux_[i].cellIndex = i ; 
                for(size_t j=0; j < intensity_.size(); ++j){
                    const DOMatrix& IF = *intensity_[j];
                    internalHeatFlux_[i].x += IF.I().sys().X()(i)*IF.S().x*IF.S().w ;
                    internalHeatFlux_[i].y += IF.I().sys().X()(i)*IF.S().y*IF.S().w ;
                    internalHeatFlux_[i].z += IF.I().sys().X()(i)*IF.S().z*IF.S().w ;
                }
            }

            for (Patch p : mesh_.boundaryPatches()){

                if (intensity_.front()->iabf().boundary()[p.index()]->typeName() == "periodic"){
                    for(int faceIndex : p.faces()){
                        wallHeatFlux_[faceIndex].faceIndex = faceIndex ;
                        wallHeatFlux_[faceIndex].x = 0.0 ;
                        wallHeatFlux_[faceIndex].y = 0.0 ;
                        wallHeatFlux_[faceIndex].z = 0.0 ;
                        wallNormalHeatFlux_[faceIndex] = 0.0 ;
                    }
                    continue;
                }
                else{
                    for (int faceIndex : p.faces()){
                        
                        wallHeatFlux_[faceIndex].faceIndex =
                            faceIndex;

                        const int ownerCell =
                            mesh_.faceOwnerNeighbour(faceIndex)[0];

                        const array<double, 3> faceNormal =
                            mesh_.faceNormalVector(faceIndex);

                        for (size_t j = 0; j < intensity_.size(); ++j)
                        {
                            const DOMatrix& IF = *intensity_[j];
                            const double directionCosine = IF.S().x * faceNormal[0] + IF.S().y * faceNormal[1] + IF.S().z * faceNormal[2];
                            double faceIntensity = 0.0;

                            if (directionCosine < 0.0)
                            {
                                faceIntensity = IF.iabf().boundary()[p.index()]->boundaryValues()[faceIndex];
                            }
                            else
                            {
                                faceIntensity = IF.I().sys().X()(ownerCell);
                            }

                            wallHeatFlux_[faceIndex].x +=
                                faceIntensity*IF.S().x*IF.S().w;

                            wallHeatFlux_[faceIndex].y +=
                                faceIntensity*IF.S().y*IF.S().w;

                            wallHeatFlux_[faceIndex].z +=
                                faceIntensity*IF.S().z*IF.S().w;
                        }

                        wallNormalHeatFlux_[faceIndex] = -(wallHeatFlux_[faceIndex].x * faceNormal[0] + wallHeatFlux_[faceIndex].y * faceNormal[1] + wallHeatFlux_[faceIndex].z * faceNormal[2]) ;
                        }
                    }
                
            }
        }

        void calculateIrradiation(){

            internalIrradiation_.assign(mesh_.numCells(), 0.0) ;
            wallIrradiation_.assign(mesh_.numFaces(), 0.0) ;

            for(int i=0; i<mesh_.numCells(); ++i){ 
                for(size_t j=0; j < intensity_.size(); ++j){
                    const DOMatrix& IF = *intensity_[j];
                    internalIrradiation_[i] += IF.I().sys().X()(i)*IF.S().w ;
                }
            }

            for(Patch p: mesh_.boundaryPatches()){
                if (intensity_.front()->iabf().boundary()[p.index()]->typeName() == "periodic"){
                    for(int faceIndex: p.faces()){
                        wallIrradiation_[faceIndex] = 0.0 ;
                    }
                    continue ;
                }
                else{
                    for(int faceIndex: p.faces()){
                        array<double,3> n = mesh_.faceNormalVector(faceIndex) ;
                        const int ownerCell = mesh_.faceOwnerNeighbour(faceIndex)[0] ;
                        for(size_t j=0; j < intensity_.size(); ++j){
                            const DOMatrix& IF = *intensity_[j];
                            wallIrradiation_[faceIndex] += IF.I().sys().X()(ownerCell)*IF.S().w*abs(max(IF.S().x*n[0] + IF.S().y*n[1] + IF.S().z*n[2], 0.0)) ;
                        }
                    }
                }
            }
        }

        void calculateAverages(){

            for(Patch p: mesh_.boundaryPatches()){

                double totalHeatTransfer = 0 ;
                double totalIrradiation = 0 ;
                double totalArea = 0 ;

                for(int faceIndex: p.faces()){
                    double area = mesh_.faceArea(faceIndex) ;
                    totalHeatTransfer += wallNormalHeatFlux_[faceIndex]*area ;
                    totalIrradiation += wallIrradiation_[faceIndex]*area ;
                    totalArea += area ;
                }

                if (totalArea <= 0.0){
                    throw runtime_error(
                        "Cannot calculate DOM wall averages for zero-area patch: "
                        + p.name()
                    );
                }

                averageWallNormalHeatFlux_[p.name()] = totalHeatTransfer/totalArea ;
                averageWallIrradiation_[p.name()] = totalIrradiation/totalArea ;
                totalWallHeatTransfer_[p.name()] = totalHeatTransfer ;

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
                    else{ dimensionlessWallTemperatureSlip_[patchNames[i]] = 0 ;}                
                }
            }

        }

        void calculateCenterlineTemperatureAndGInX(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int jc = numY/2 ;
            int kc = numZ/2 ;

            for(int i=0; i<mesh_.numCellsXYZ()[0]; ++i){
                int cellIndex = kc*numX*numY + jc*numX + i ;
                double x = mesh_.cellCentre(cellIndex)[0];
                double valueT = fvMatrixT_.sys().X()(cellIndex) ;
                double valueG = gInternal_(cellIndex) ;
                centerlineTemperatureInX_[x] = valueT ;
                centerlineGInX_[x] = valueG ;
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
                double valueG = gInternal_(cellIndex);

                centerlineTemperatureInZ_[z] = valueT;
                centerlineGInZ_[z] = valueG;
            }
        }

        void calculateCenterlineFluxInX(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int jc = numY/2 ;
            int kc = numZ/2 ;

            centerlineFluxInX_.clear() ;

            for(int i=0; i<numX; ++i){
                int cellIndex = kc*numX*numY + jc*numX + i ;
                double x = mesh_.cellCentre(cellIndex)[0] ;
                centerlineFluxInX_[x] = internalHeatFlux_[cellIndex].x ;
            }
        }

        void calculateCenterlineFluxInZ(){

            int numX = mesh_.numCellsXYZ()[0] ;
            int numY = mesh_.numCellsXYZ()[1] ;
            int numZ = mesh_.numCellsXYZ()[2] ;

            int ic = numX/2 ;
            int jc = numY/2 ;

            centerlineFluxInZ_.clear() ;

            for(int k=0; k<numZ; ++k){
                int cellIndex = k*numX*numY + jc*numX + ic ;
                double z = mesh_.cellCentre(cellIndex)[2] ;
                centerlineFluxInZ_[z] = internalHeatFlux_[cellIndex].z ;
            }
        }

        void calculateDimensionlessCenterlineTemperatureInX(){

            for(const auto& [key,val]: centerlineTemperatureInX_){
                dimensionlessCenterlineTemperatureInX_[key] = abs((powl(val,4) - powl(avgWallTemperature_["east"],4))/(powl(avgWallTemperature_["west"],4) - powl(avgWallTemperature_["east"],4))) ;
            }

        }

        void calculateDimensionlessCenterlineTemperatureInX2(){

            for(const auto& [key,val]: centerlineTemperatureInX_){
                dimensionlessCenterlineTemperatureInX2_[key] = abs((powl(val,4))/(powl(avgWallTemperature_["east"],4))) ;
            }

        }

        void calculateDimensionlessCenterlineTemperatureInZ(){

            dimensionlessCenterlineTemperatureInZ_.clear();
            const double Tw = avgWallTemperature_.at("bottom");

            if(std::abs(Tw) < 1e-14){
                throw std::runtime_error("Cannot calculate dimensionless centerline temperature in Z: avgWallTemperature_[\"bottom\"] is zero.");
            }

            for(const auto& [z, T] : centerlineTemperatureInZ_){
                dimensionlessCenterlineTemperatureInZ_[z] = T/Tw;
            }

        }

        void calculateCenterlineRadiativeSourceTermInX(){

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

        map<string,double> parallelPlateHeatFluxBlackbodyAnalyticalValue(double Twest, double Teast){

            map<string,double> G ;
            double val = sigma*(powl(Twest,4) - powl(Teast,4)) ;
            
            G["west"] = val ;
            G["east"] = -val ;

            return G ; 
        }

        void calculateWestEastTransferForParallelPlates(){

            oppositeWallHeatFlux_.clear();
            oppositeWallHeatTransfer_.clear();

            const double qWest = averageWallNormalHeatFlux_.at("west");
            const double qEast = averageWallNormalHeatFlux_.at("east");

            const double QWest = totalWallHeatTransfer_.at("west");
            const double QEast = totalWallHeatTransfer_.at("east");

            const double qWestToEast = 0.5*(qWest - qEast);
            const double QWestToEast = 0.5*(QWest - QEast);

            oppositeWallHeatFlux_["west"] = qWestToEast;
            oppositeWallHeatFlux_["east"] = -qWestToEast;

            oppositeWallHeatTransfer_["west"] = QWestToEast;
            oppositeWallHeatTransfer_["east"] = -QWestToEast;
        }

        void calculateDimensionlessHeatFlux(){

            vector<string> patchNames ;
            for(Patch p: mesh_.boundaryPatches()){
                patchNames.push_back(p.name()) ;
            }

            for(int i=0 ; i< mesh_.numBoundaryPatches(); ++i){
                if(i%2 == 0){
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i+1]],4) ;
                    if(abs(denom) > 1e-14){dimensionlessqw_[patchNames[i]] = abs(averageWallNormalHeatFlux_[patchNames[i]]/(sigma*(denom))) ;}
                    else{dimensionlessqw_[patchNames[i]] = 0;}
                }
                else{
                    double denom = powl(avgWallTemperature_[patchNames[i]],4) - powl(avgWallTemperature_[patchNames[i-1]],4) ;
                    if(abs(denom) > 1e-14){dimensionlessqw_[patchNames[i]] = abs(averageWallNormalHeatFlux_[patchNames[i]]/(sigma*(denom))) ;}
                    else{dimensionlessqw_[patchNames[i]] = 0;}
                }
            }
        }

        void calculateDimensionlessHeatFlux2(){

            dimensionlessqw2_.clear() ;

            vector<string> patchNames ;
            for(Patch p: mesh_.boundaryPatches()){
                patchNames.push_back(p.name()) ;
            }

            for(int i=0 ; i<mesh_.numBoundaryPatches(); ++i){
                double denom = powl(avgWallTemperature_[patchNames[i]],4) ;
                if(abs(denom) > 1e-14){
                    dimensionlessqw2_[patchNames[i]] = abs(averageWallNormalHeatFlux_[patchNames[i]]/(sigma*denom)) ;
                }
                else{
                    dimensionlessqw2_[patchNames[i]] = 0 ;
                }
            }
        }
        
        void calculateRadiativeSourceTerm(ScalarProperty& a){
            radiativeSourceTerm_.assign(mesh_.numCells(), 0.0) ;
            for(int i = 0; i < mesh_.numCells(); ++i){
                radiativeSourceTerm_[i] = a(i)*(4*sigma*pow(fvMatrixT_.sys().X()(i),4) - gInternal_(i)) ;
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
