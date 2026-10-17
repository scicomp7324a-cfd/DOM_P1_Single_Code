# pragma once

# include "PathConfig.hpp"

# include<iostream>
# include<vector>
# include<array>
# include<map>
# include<chrono>

# include "MeshFileReader3D.hpp"
# include "FaceAddressedMesh3D.hpp"

# include "SparseAddress.hpp"
# include "SparseMatrix.hpp"
# include "SparseLinearVector.hpp"
# include "SparseLinearSystem.hpp"
# include "GaussSeidelSmoother.hpp"
# include "GaussSeidelSolver.hpp"

# include "Convection.hpp"
# include "Diffusion.hpp"
# include "Ddt.hpp"
# include "Source.hpp"
# include "FVMatrix.hpp"
# include "FaceFluxField.hpp"
# include "ScalarProperty.hpp"
# include "BoundaryCondition.hpp"
# include "FixedValueBC.hpp"
# include "GradientBC.hpp"
# include "Patch.hpp"
# include "FieldFileReader.hpp"
# include "Field.hpp"
# include "PressureLaplacian.hpp"
# include "GaussGradient.hpp"
# include "UnderRelaxation.hpp"
# include "WallHeatFlux.hpp"

using namespace std ;


class P1{

    private:

        const vector<string> fieldFiles_ ;

        const SparseAddress addr_ ;
        const FaceAddressedMesh3D mesh_ ;

        FVMatrix G_ ;
        FVMatrix T_ ;

        string n_ ;
        vector<int> numCellsXYZ_ ;
        vector<double> dimensionsXYZ_ ;

        double kVal_ ;
        double aVal_ ;

        ScalarProperty a_ ;
        ScalarProperty k_ ;
        ScalarProperty sigmaS_ ;
        ScalarProperty gamma_ ;
        WallProperty emissivity_ ;

        FaceFluxField flux_ ;

        Field Tinit_ ;
        Field Ginit_ ;

        SparseLinearVector sourceDiagonalG_ ;
        SparseLinearVector sourceTermG_ ;

        SparseLinearVector sourceDiagonalT_ ;
        SparseLinearVector sourceTermT_ ;

        double sigma = 5.670374419e-8  ;
        double relativeTolerance = 0.0000001 ;
        double changeTolerance = 1e-6 ;

        unique_ptr<WallHeatFlux> heatData_;

        string caseName_ ;
        string subCaseName_ ;

        double anisotropyFactor_ = 0 ;

    public:

        P1(vector<string> meshFiles, vector<string> fieldFiles, double a, double sigmaS, double k ,vector<double> emissivity, string dimension, string caseName, string subCaseName)
        :fieldFiles_(fieldFiles),
         addr_(meshFiles),
         mesh_(addr_.mesh()),
         G_(addr_),
         T_(addr_),
         n_(dimension),
         numCellsXYZ_(mesh_.numCellsXYZ()),
         dimensionsXYZ_(mesh_.dimensionXYZ()),
         kVal_(k),
         aVal_(a),
         a_(mesh_, a),
         k_(mesh_, k),
         sigmaS_(mesh_, sigmaS),
         gamma_(mesh_, gamma(a, sigmaS)),
         emissivity_(mesh_, emissivity),
         flux_(addr_),
         Tinit_(),
         Ginit_(),
         sourceDiagonalG_(mesh_.numCells()),
         sourceTermG_(mesh_.numCells()),
         sourceDiagonalT_(mesh_.numCells()),
         sourceTermT_(mesh_.numCells()),
         caseName_(caseName),
         subCaseName_(subCaseName)
         {}

        P1(vector<string> meshFiles, vector<string> fieldFiles, double a, double sigmaS, double k ,vector<double> emissivity, string dimension, string caseName, string subCaseName, double anisotropyFactor)
        :fieldFiles_(fieldFiles),
         addr_(meshFiles),
         mesh_(addr_.mesh()),
         G_(addr_),
         T_(addr_),
         n_(dimension),
         numCellsXYZ_(mesh_.numCellsXYZ()),
         dimensionsXYZ_(mesh_.dimensionXYZ()),
         kVal_(k),
         aVal_(a),
         a_(mesh_, a),
         k_(mesh_, k),
         sigmaS_(mesh_, sigmaS),
         gamma_(mesh_, gammaScattering(a, sigmaS, anisotropyFactor)),
         emissivity_(mesh_, emissivity),
         flux_(addr_),
         Tinit_(),
         Ginit_(),
         sourceDiagonalG_(mesh_.numCells()),
         sourceTermG_(mesh_.numCells()),
         sourceDiagonalT_(mesh_.numCells()),
         sourceTermT_(mesh_.numCells()),
         caseName_(caseName),
         subCaseName_(subCaseName),
         anisotropyFactor_(anisotropyFactor)
         {}

        const unique_ptr<WallHeatFlux>& heatData() const{ return heatData_ ; }
        unique_ptr<WallHeatFlux>& heatData() { return heatData_ ;}

        void runP1(){

            if(kVal_ == 0 && aVal_ ==0){
                runP1PureRadiationNoParticipatingMedium() ;
            }
            else if(kVal_ == 0 && aVal_ != 0){
                runP1PureRadiationParticipatingMedium() ;
            }
            else{
                runP1RadiationConduction() ;
            }

        }

        void runP1PureRadiationNoParticipatingMedium(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;

            double maxFieldChangeG = -1e6 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 1000000000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;

                diffusionG.assemble() ;
                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }

            cout << "West wall heat flux: " << heatData_->avgqw()["west"] << endl;
            cout << "East wall heat flux: " << heatData_->avgqw()["east"] << endl;

            writeVTKRectilinearFileGT(outSimulation) ;
            
        }

        void runP1PureRadiationParticipatingMedium(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;
                T_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;

                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;
                updateT() ;
                setSourceVectorsG() ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;

        }

        void runP1PureRadiationParticipatingMediumTConstant(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;
                T_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;

                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;

                setSourceVectorsG() ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;

        }

        void runP1PureRadiationParticipatingMediumTConstantIsotropicScattering(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;
                T_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;

                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;

                //setSourceVectorsG() ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();

            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTerm() ;
            

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;

        }

        void runP1RadiationConduction(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(n_) ;

            ofstream outInnerG(paths.getOutputDirectory(n_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(n_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outInnerT(paths.getOutputDirectory(n_) / "innerResidualT.csv") ;
            outInnerT << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalT(paths.getOutputDirectory(n_) / "finalRedidualT.csv") ;
            outFinalT << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(n_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(n_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(n_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;
            setSourceVectorsT() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            Diffusion diffusionT(T_, mesh_, Tinit_, k_) ;
            Source sourceT(T_, mesh_, sourceDiagonalT_, sourceTermT_) ;
            
            double maxFieldChangeT = -1e4 ;
            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;
                T_.clear() ;

                SparseLinearVector Tchange_ = T_.sys().X() ;
                SparseLinearVector Gchange_ = G_.sys().X() ;

                diffusionG.assemble() ;
                sourceG.assemble() ;

                diffusionT.assemble() ;
                sourceT.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;
                T_.solve(relativeTolerance, true, outInnerT, outFinalT, count) ;

                setSourceVectorsG() ;
                setSourceVectorsT() ;

                Tchange_ = Tchange_ - T_.sys().X() ;
                Gchange_ = Gchange_ - G_.sys().X() ;

                maxFieldChangeT = Tchange_.maxAbsVal() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;

                fieldsAreChanging = !((maxFieldChangeT < changeTolerance) && (maxFieldChangeG < changeTolerance)) && (count < maxIter);                
                cout << "maxFieldChangeG: " << maxFieldChangeG << " maxFieldChangeT: " << maxFieldChangeT << endl ;
                ++count ;
                
            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();
            heatData_->calculateCenterlineTemperatureAndGInZ() ;
            heatData_->calculateDimensionlessCenterlineTemperatureInZ() ;

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }

            writeVTKRectilinearFileGT(outSimulation) ;
        }

        void runP1RadiationAnnularSourceRectangularEnclosure(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;
            updateAnnularTemperature(0.20, 0.25, 1000) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;
                setSourceVectorsG() ;
                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTerm() ;

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;
            
        }

        void runP1RadiationHotWallColdMediumRectangularEnclosure(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                G_.clear() ;

                SparseLinearVector Gchange_ = G_.sys().X() ;
                setSourceVectorsG() ;
                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTerm() ;

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;
    
        }

        void runP1RadiationScattering(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            FieldFileReader GReader(mesh_, G_, flux_, gamma_, emissivity_, Tinit_, fieldFiles_[1]) ;            
            Ginit_ = Field(GReader.internal(), std::move(GReader.boundary())) ;

            heatData_ = make_unique<WallHeatFlux>(mesh_, G_, T_, Ginit_, Tinit_, gamma_);

            PathConfig paths ;
            paths.createOutputDirectory(caseName_, subCaseName_) ;

            ofstream outInnerG(paths.getOutputDirectory(caseName_, subCaseName_) / "innerResidualG.csv") ;
            outInnerG << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalG(paths.getOutputDirectory(caseName_, subCaseName_) / "finalRedidualG.csv") ;
            outFinalG << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outXG(paths.getOutputDirectory(caseName_, subCaseName_) / "XGFinal.csv") ;
            ofstream outXT(paths.getOutputDirectory(caseName_, subCaseName_) / "XTFinal.csv") ;

            ofstream outSimulation(paths.getOutputDirectory(caseName_, subCaseName_) / "simulation.vtk") ;

            G_.clear() ;
            T_.clear() ;
            G_.initializeSolution(Ginit_.internal()) ;
            T_.initializeSolution(Tinit_.internal()) ;

            setSourceVectorsG() ;

            Diffusion diffusionG(G_, mesh_, Ginit_, gamma_) ;
            Source sourceG(G_, mesh_, sourceDiagonalG_, sourceTermG_) ;

            double maxFieldChangeG = -1e4 ;

            bool fieldsAreChanging = true ;
            int count = 0 ;
            int maxIter = 10000 ;

            while(fieldsAreChanging){

                cout << "Loop count: " << count << endl ;
                SparseLinearVector Gchange_ = G_.sys().X() ;

                G_.clear() ;
                T_.clear() ;

                diffusionG.assemble() ;
                sourceG.assemble() ;

                G_.solve(relativeTolerance, true, outInnerG, outFinalG, count) ;
                updateT() ;
                setSourceVectorsG() ;

                Gchange_ = Gchange_ - G_.sys().X() ;
                maxFieldChangeG = Gchange_.maxAbsVal() ;
                fieldsAreChanging = !((maxFieldChangeG < changeTolerance)) && (count < maxIter);       

                cout << "maxFieldChangeG: " << maxFieldChangeG << endl ;
                ++count ;

            }

            for(Patch p : mesh_.boundaryPatches()){
                Ginit_.boundary()[p.index()]->updateBoundaryValues();
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateWallHeatFlux();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessQw();
            heatData_->calculateDimensionlessQw2();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndG();
            heatData_->calculateDimensionlessCenterlineTemperature();
            heatData_->calculateDimensionlessCenterlineTemperature2();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTerm() ;
            

            for(int i = 0; i < mesh_.numCells(); ++i){
                outXG << i << "," << G_.sys().X()(i) << endl ;
                outXT << i << "," << T_.sys().X()(i) << endl ;
            }
            writeVTKRectilinearFileGT(outSimulation) ;

        }




        // Utility function

        void setSourceVectorsG(){
            for(int i=0; i<mesh_.numCells(); ++i){
                sourceDiagonalG_(i) = a_(i) ;
                sourceTermG_(i) = 4*sigma*a_(i)*powl(T_.sys().X()(i),4) ;
            }
        }

        void setSourceVectorsT(){
            for(int i=0; i<mesh_.numCells(); ++i){
                sourceDiagonalT_(i) = 0 ;
                sourceTermT_(i) = a_(i)*G_.sys().X()(i) - 4*sigma*a_(i)*powl(T_.sys().X()(i),4) ;
            }
        }

        void updateT(){
            for(int i=0; i<mesh_.numCells(); ++i){
                T_.sys().X()(i) = powl((G_.sys().X()(i)/(4*sigma)), 0.25) ;
            }
        }


        void writeVTKRectilinearFileGT(std::ostream& out){

            if(!out){
                throw std::runtime_error("Could not open VTK file for writing") ;
            }

            const int numActiveCells = mesh_.numCells() ;
            const int numX = numCellsXYZ_[0] ;
            const int numY = numCellsXYZ_[1] ;
            const int numZ = numCellsXYZ_[2] ;
            const int numStructuredCells = numX*numY*numZ ;

            if(numX<=0 || numY<=0 || numZ<=0){
                throw std::runtime_error("Invalid structured mesh dimensions while writing P1 VTK file") ;
            }

            const double dx = dimensionsXYZ_[0]/static_cast<double>(numX) ;
            const double dy = dimensionsXYZ_[1]/static_cast<double>(numY) ;
            const double dz = dimensionsXYZ_[2]/static_cast<double>(numZ) ;

            /*
                The VTK file retains the complete Cartesian lattice.

                Active medium cells are mapped to their original Cartesian
                positions. Deleted cells inside the circular pipes remain
                ghost cells and are hidden by ParaView.
            */

            vector<int> activeCellToStructuredCell(numActiveCells, -1) ;
            vector<int> structuredCellToActiveCell(numStructuredCells, -1) ;
            vector<unsigned char> ghostCell(numStructuredCells, 32) ;

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                const array<double,3> centre = mesh_.cellCentre(cellIndex) ;

                int i = static_cast<int>(round(centre[0]/dx-0.5)) ;
                int j = static_cast<int>(round(centre[1]/dy-0.5)) ;
                int k = static_cast<int>(round(centre[2]/dz-0.5)) ;

                if(i<0 || i>=numX || j<0 || j>=numY || k<0 || k>=numZ){
                    throw std::runtime_error(
                        "Active mesh cell could not be mapped to the Cartesian VTK lattice. Cell index: "
                        + to_string(cellIndex)
                    ) ;
                }

                const int structuredCellIndex = k*numX*numY+j*numX+i ;

                if(structuredCellToActiveCell[structuredCellIndex] != -1){
                    throw std::runtime_error(
                        "More than one active mesh cell maps to the same Cartesian VTK cell. Structured cell index: "
                        + to_string(structuredCellIndex)
                    ) ;
                }

                activeCellToStructuredCell[cellIndex] = structuredCellIndex ;
                structuredCellToActiveCell[structuredCellIndex] = cellIndex ;
                ghostCell[structuredCellIndex] = 0 ;
            }

            /*
                Collect the final boundary values of G using global face indices.
            */

            map<int,double> boundaryGValues ;

            for(Patch p: mesh_.boundaryPatches()){

                const map<int,double> patchValues =
                    Ginit_.boundary()[p.index()]->boundaryValues() ;

                for(const auto& [faceIndex,value]: patchValues){
                    boundaryGValues[faceIndex] = value ;
                }
            }

            vector<array<double,3>> activeRadiativeHeatFlux(
                numActiveCells,
                {0.0,0.0,0.0}
            ) ;

            vector<double> activeRadiativeHeatFluxMagnitude(
                numActiveCells,
                0.0
            ) ;

            vector<double> activeRadiativeSourceTerm(
                numActiveCells,
                0.0
            ) ;

            /*
                Calculate the cell-centred gradient using the Gauss theorem.

                This uses the actual mesh connectivity and therefore remains
                valid around the deleted pipe cells.
            */

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                double gradX = 0.0 ;
                double gradY = 0.0 ;
                double gradZ = 0.0 ;

                const double cellVolume = mesh_.cellVolume(cellIndex) ;

                if(cellVolume<=0.0){
                    throw std::runtime_error(
                        "Non-positive cell volume while calculating P1 VTK diagnostics. Cell index: "
                        + to_string(cellIndex)
                    ) ;
                }

                for(int faceIndex: mesh_.cellFaces(cellIndex)){

                    const array<int,2> ownerNeighbour =
                        mesh_.faceOwnerNeighbour(faceIndex) ;

                    array<double,3> faceAreaVector =
                        mesh_.faceAreaVector(faceIndex) ;

                    /*
                        mesh_.faceAreaVector() points from owner to neighbour,
                        or outwards from the owner for a boundary face.

                        Reverse it when the current cell is the neighbour.
                    */

                    if(ownerNeighbour[0] != cellIndex){
                        faceAreaVector[0] = -faceAreaVector[0] ;
                        faceAreaVector[1] = -faceAreaVector[1] ;
                        faceAreaVector[2] = -faceAreaVector[2] ;
                    }

                    double faceG = 0.0 ;

                    if(ownerNeighbour[1] != -1){

                        faceG = 0.5*(
                            G_.sys().X()(ownerNeighbour[0])
                            +
                            G_.sys().X()(ownerNeighbour[1])
                        ) ;
                    }
                    else{

                        const auto boundaryValueIterator =
                            boundaryGValues.find(faceIndex) ;

                        if(boundaryValueIterator == boundaryGValues.end()){
                            throw std::runtime_error(
                                "Boundary G value was not found for face "
                                + to_string(faceIndex)
                            ) ;
                        }

                        faceG = boundaryValueIterator->second ;
                    }

                    gradX += faceG*faceAreaVector[0]/cellVolume ;
                    gradY += faceG*faceAreaVector[1]/cellVolume ;
                    gradZ += faceG*faceAreaVector[2]/cellVolume ;
                }

                const double qx = -gamma_(cellIndex)*gradX ;
                const double qy = -gamma_(cellIndex)*gradY ;
                const double qz = -gamma_(cellIndex)*gradZ ;

                activeRadiativeHeatFlux[cellIndex] = {qx,qy,qz} ;

                activeRadiativeHeatFluxMagnitude[cellIndex] =
                    sqrt(qx*qx+qy*qy+qz*qz) ;

                const double temperature =
                    T_.sys().X()(cellIndex) ;

                const double incidentRadiation =
                    G_.sys().X()(cellIndex) ;

                activeRadiativeSourceTerm[cellIndex] =
                    a_(cellIndex)*(
                        4.0*sigma*pow(temperature,4)
                        -
                        incidentRadiation
                    ) ;
            }

            /*
                Store active-cell values on the complete structured lattice.

                Values belonging to deleted pipe cells remain zero, but those
                cells are hidden using vtkGhostType.
            */

            vector<double> structuredG(numStructuredCells, 0.0) ;
            vector<double> structuredT(numStructuredCells, 0.0) ;

            vector<array<double,3>> structuredRadiativeHeatFlux(
                numStructuredCells,
                {0.0,0.0,0.0}
            ) ;

            vector<double> structuredRadiativeHeatFluxMagnitude(
                numStructuredCells,
                0.0
            ) ;

            vector<double> structuredRadiativeSourceTerm(
                numStructuredCells,
                0.0
            ) ;

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                const int structuredCellIndex =
                    activeCellToStructuredCell[cellIndex] ;

                structuredG[structuredCellIndex] =
                    G_.sys().X()(cellIndex) ;

                structuredT[structuredCellIndex] =
                    T_.sys().X()(cellIndex) ;

                structuredRadiativeHeatFlux[structuredCellIndex] =
                    activeRadiativeHeatFlux[cellIndex] ;

                structuredRadiativeHeatFluxMagnitude[structuredCellIndex] =
                    activeRadiativeHeatFluxMagnitude[cellIndex] ;

                structuredRadiativeSourceTerm[structuredCellIndex] =
                    activeRadiativeSourceTerm[cellIndex] ;
            }

            out << "# vtk DataFile Version 3.0\n" ;
            out << "P1 heat exchanger radiation result\n" ;
            out << "ASCII\n" ;
            out << "DATASET RECTILINEAR_GRID\n" ;

            out << "DIMENSIONS "
                << numX+1 << " "
                << numY+1 << " "
                << numZ+1 << "\n\n" ;

            out << "X_COORDINATES " << numX+1 << " double\n" ;

            for(int i=0; i<=numX; ++i){
                out << setprecision(12)
                    << dx*static_cast<double>(i)
                    << (i<numX ? " " : "\n") ;
            }

            out << "\n" ;

            out << "Y_COORDINATES " << numY+1 << " double\n" ;

            for(int j=0; j<=numY; ++j){
                out << setprecision(12)
                    << dy*static_cast<double>(j)
                    << (j<numY ? " " : "\n") ;
            }

            out << "\n" ;

            out << "Z_COORDINATES " << numZ+1 << " double\n" ;

            for(int k=0; k<=numZ; ++k){
                out << setprecision(12)
                    << dz*static_cast<double>(k)
                    << (k<numZ ? " " : "\n") ;
            }

            out << "\n" ;

            /*
                CELL_DATA must equal the number of cells represented by the
                complete rectilinear grid, not only the number of active cells.
            */

            out << "CELL_DATA " << numStructuredCells << "\n\n" ;

            /*
                VTK hidden-cell flag:
                    0  -> active cell
                    32 -> hidden/deleted cell

                ParaView uses this array to remove the cells lying inside
                the circular pipes from the displayed domain.
            */

            out << "SCALARS vtkGhostType unsigned_char 1\n" ;
            out << "LOOKUP_TABLE default\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << static_cast<int>(
                    ghostCell[structuredCellIndex]
                ) << "\n" ;
            }

            out << "\n" ;

            out << "SCALARS G double 1\n" ;
            out << "LOOKUP_TABLE default\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << setprecision(12)
                    << structuredG[structuredCellIndex]
                    << "\n" ;
            }

            out << "\n" ;

            out << "SCALARS temperature double 1\n" ;
            out << "LOOKUP_TABLE default\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << setprecision(12)
                    << structuredT[structuredCellIndex]
                    << "\n" ;
            }

            out << "\n" ;

            out << "VECTORS radiativeHeatFlux double\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << setprecision(12)
                    << structuredRadiativeHeatFlux[structuredCellIndex][0]
                    << " "
                    << structuredRadiativeHeatFlux[structuredCellIndex][1]
                    << " "
                    << structuredRadiativeHeatFlux[structuredCellIndex][2]
                    << "\n" ;
            }

            out << "\n" ;

            out << "SCALARS radiativeHeatFluxMagnitude double 1\n" ;
            out << "LOOKUP_TABLE default\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << setprecision(12)
                    << structuredRadiativeHeatFluxMagnitude[
                        structuredCellIndex
                    ]
                    << "\n" ;
            }

            out << "\n" ;

            out << "SCALARS radiativeSourceTerm double 1\n" ;
            out << "LOOKUP_TABLE default\n" ;

            for(int structuredCellIndex=0;
                structuredCellIndex<numStructuredCells;
                ++structuredCellIndex){

                out << setprecision(12)
                    << structuredRadiativeSourceTerm[
                        structuredCellIndex
                    ]
                    << "\n" ;
            }

            out << "\n" ;
        }
          

        double gamma(double a, double sigmaS){

            if(a == 0 && sigmaS == 0){
                double L = this->dimensionsXYZ_[0] ;
                return L*100 ;
            }
            else{
                return (double)(1/(3*(a+sigmaS))) ;
            }

        }

        double gammaScattering(double a, double sigmaS, double anisotropyFactor){

            if(abs(a) < 1e-14 && abs(sigmaS) < 1e-14){
                double L = this->dimensionsXYZ_[0] ;
                return L*1000 ;
            }
            else{
                return 1.0/(3.0*(a+sigmaS*(1.0-(anisotropyFactor/3.0)))) ;
            }

        }

        void updateAnnularTemperature(double minRadius, double maxRadius, double temperature){
            array<double,3> centre = mesh_.findMeshCentre() ;
            for(int i = 0; i < mesh_.numCells(); ++i){
                array<double,3> cellCentre = mesh_.cellCentre(i) ;
                double distance = sqrt(pow((cellCentre[0]-centre[0]),2) + pow((cellCentre[1]-centre[1]),2)) ;
                if(distance >= minRadius && distance <= maxRadius){
                    T_.sys().X()(i) = temperature ;
                }
            }
        }
        
} ;

