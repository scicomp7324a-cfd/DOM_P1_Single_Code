# pragma once

# include "PathConfig.hpp"

# include<iostream>
# include<vector>
# include<array>
# include<map>
# include<chrono>
# include<iomanip>
# include<string>
# include<memory>

# include "MeshFileReader3D.hpp"
# include "FaceAddressedMesh3D.hpp"
# include "Patch.hpp"

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

# include "FieldFileReader.hpp"

# include "FixedValueBC.hpp"
# include "GradientBC.hpp"
# include "MarshakRadiationBC.hpp"
# include "PeriodicBC.hpp"
# include "DiscreteOrdinateBC.hpp"

# include "Field.hpp"
# include "PressureLaplacian.hpp"
# include "GaussGradient.hpp"
# include "UnderRelaxation.hpp"

# include "WallHeatFlux.hpp"

# include "OrdinatesConstructor.hpp"
# include "DOMatrix.hpp"
# include "OrdinateField.hpp"

using namespace std ;

class DOM{

    private:

        const vector<string> fieldFiles_ ;

        const SparseAddress addr_ ;
        const FaceAddressedMesh3D mesh_ ;

        FVMatrix T_ ;
        Field Tinit_ ;

        SparseLinearVector sourceDiagonalT_ ;
        SparseLinearVector sourceTermT_ ;

        SparseLinearVector gInternal_ ;
        SparseLinearVector gFace_ ;

        PathConfig paths_ ;

        string meshDimensions_ ;
        vector<int> numCellsXYZ_ ;
        vector<double> dimensionsXYZ_ ;

        double aVal_ ;
        double kVal_ ;
        double sigmaSVal_ ;

        ScalarProperty a_ ;         // Absorption coefficient
        ScalarProperty k_ ;         // Conductivity
        ScalarProperty sigmaS_ ;    // Scattering coefficient
        ScalarProperty gamma_ ;

        WallProperty emissivity_ ;

        int ordinate_ ;             // Choice of the ordinate S2=2, S4=4, S6=6

        vector<SparseLinearVector> simulationXI_ ;  // Simulation vector accumulator
        vector<OrdinateS> ordinateDirections_ ;    // Simulation ordinate directions

        vector<unique_ptr<DOMatrix>> intensities_ ;
        FaceFluxField flux_ ;

        double sigma = 5.670374419e-8  ;
        double pi = 3.1415926535897932 ;
        double relativeTolerance = 0.000000001 ;
        double changeTolerance = 1e-7 ;

        unique_ptr<DiscreteOrdinateWallHeatFlux> heatData_ ;

        string caseName_ ;
        string subCaseName_ ;

        double anisotropyFactor_ = 0 ;

    public:

        DOM(vector<string> meshFiles, vector<string> fieldFiles, double a, double k, double sigmaS, vector<double> emissivity, string dimension, string caseName, string subCaseName, int ordinate)
        :fieldFiles_(fieldFiles),
         addr_(meshFiles),
         mesh_(addr_.mesh()),
         T_(addr_),
         Tinit_(),
         sourceDiagonalT_(mesh_.numCells()),
         sourceTermT_(mesh_.numCells()),
         gInternal_(mesh_.numCells()),
         gFace_(mesh_.numFaces()),
         meshDimensions_(dimension),
         numCellsXYZ_(mesh_.numCellsXYZ()),
         dimensionsXYZ_(mesh_.dimensionXYZ()),
         aVal_(a),
         kVal_(k),
         sigmaSVal_(sigmaS),
         a_(mesh_, a),
         k_(mesh_, k),
         sigmaS_(mesh_, sigmaS),
         gamma_(mesh_, gamma(a, sigmaS)),
         emissivity_(mesh_, emissivity),
         ordinate_(ordinate),
         flux_(addr_),
         caseName_(caseName),
         subCaseName_(subCaseName)
        {}

        DOM(vector<string> meshFiles, vector<string> fieldFiles, double a, double k, double sigmaS, vector<double> emissivity, string dimension, string caseName, string subCaseName, int ordinate, double anisotropyFactor)
        :fieldFiles_(fieldFiles),
         addr_(meshFiles),
         mesh_(addr_.mesh()),
         T_(addr_),
         Tinit_(),
         sourceDiagonalT_(mesh_.numCells()),
         sourceTermT_(mesh_.numCells()),
         gInternal_(mesh_.numCells()),
         gFace_(mesh_.numFaces()),
         meshDimensions_(dimension),
         numCellsXYZ_(mesh_.numCellsXYZ()),
         dimensionsXYZ_(mesh_.dimensionXYZ()),
         aVal_(a),
         kVal_(k),
         sigmaSVal_(sigmaS),
         a_(mesh_, a),
         k_(mesh_, k),
         sigmaS_(mesh_, sigmaS),
         gamma_(mesh_, gamma(a, sigmaS)),
         emissivity_(mesh_, emissivity),
         ordinate_(ordinate),
         flux_(addr_),
         caseName_(caseName),
         subCaseName_(subCaseName),
         anisotropyFactor_(anisotropyFactor)
        {}

        unique_ptr<DiscreteOrdinateWallHeatFlux>& heatData(){ return heatData_ ;}
        const unique_ptr<DiscreteOrdinateWallHeatFlux>& heatData() const{ return heatData_ ;}

        void runDOMPureRadiationNoParticipatingMedium(){
            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]);
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary()));

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());

            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.

                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Update the source vector for the next outer iteration.
                    */
                    setSourceVectorsI(IF);

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                if(count%10==0){
                    writeVTKRectilinearIntensityFile(outSimulation);
                }
                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */
            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }
            updateGInternal() ;
            /*
                Calculate:
                    - internal radiative heat flux;
                    - wall radiative heat flux;
                    - wall-normal radiative heat flux;
                    - internal irradiation;
                    - wall irradiation;
                    - patch-averaged wall quantities.
            */
            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_,intensities_, T_, Tinit_, gInternal_, gFace_);
            
            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();

            const vector<double>& wallNormalHeatFlux = heatData_->wallNormalHeatFlux();
            const vector<double>& wallIrradiation = heatData_->wallIrradiation();
            const map<string, double>& averageWallNormalHeatFlux = heatData_->averageWallNormalHeatFlux();
            const map<string, double>& averageWallIrradiation = heatData_->averageWallIrradiation();

            /*
                Create the post-processing output file.
            */
            const filesystem::path postProcessingFilePath = outputDirectory / "DOM_PostProcessing.txt";
            ofstream postProcessingFile( postProcessingFilePath );

            if (!postProcessingFile){
                throw runtime_error("DOM: unable to create post-processing file: " + postProcessingFilePath.string());
            }

            postProcessingFile
                << scientific
                << setprecision(10);

            postProcessingFile
                << "============================================================\n"
                << "             DISCRETE ORDINATES POST-PROCESSING\n"
                << "============================================================\n\n";

            postProcessingFile
                << "Case name                     : "
                << caseName_
                << "\n";

            postProcessingFile
                << "Subcase name                  : "
                << subCaseName_
                << "\n";

            postProcessingFile
                << "Number of ordinates           : "
                << intensities_.size()
                << "\n";

            postProcessingFile
                << "Number of outer iterations    : "
                << count
                << "\n";

            postProcessingFile
                << "Final maximum intensity change: "
                << maxFieldChangeI
                << "\n";

            postProcessingFile
                << "Outer convergence tolerance   : "
                << changeTolerance
                << "\n";

            postProcessingFile
                << "Linear relative tolerance     : "
                << relativeTolerance
                << "\n\n";

            postProcessingFile
                << "------------------------------------------------------------\n"
                << "Physical wall results\n"
                << "------------------------------------------------------------\n";

            double globalIntegratedWallHeatRate = 0.0;
            double totalAbsoluteWallHeatRate = 0.0;

            /*
                Patch is again copied because Patch::index() is non-const.
            */

            for (Patch p : mesh_.boundaryPatches())
            {
                const string boundaryType = intensities_.front()->iabf().boundary()[p.index()]->typeName();

                /*
                    Periodic boundaries are not physical radiating walls.
                */
                if (boundaryType == "periodic"){
                    continue;
                }

                double patchArea = 0.0;
                double integratedNormalHeatRate = 0.0;
                double integratedIncidentIrradiation = 0.0;
                double minimumNormalHeatFlux = numeric_limits<double>::max();
                double maximumNormalHeatFlux = numeric_limits<double>::lowest();
                double minimumIrradiation = numeric_limits<double>::max();
                double maximumIrradiation = numeric_limits<double>::lowest();

                for (int faceIndex : p.faces())
                {
                    const double faceArea = mesh_.faceArea(faceIndex);
                    const double normalHeatFlux = wallNormalHeatFlux[faceIndex];
                    const double irradiation = wallIrradiation[faceIndex];

                    patchArea += faceArea;

                    integratedNormalHeatRate += normalHeatFlux * faceArea;
                    integratedIncidentIrradiation += irradiation * faceArea;

                    minimumNormalHeatFlux = min(minimumNormalHeatFlux, normalHeatFlux);
                    maximumNormalHeatFlux = max(maximumNormalHeatFlux,normalHeatFlux);
                    minimumIrradiation = min(minimumIrradiation,irradiation);
                    maximumIrradiation = max(maximumIrradiation,irradiation);
                }

                globalIntegratedWallHeatRate += integratedNormalHeatRate;
                totalAbsoluteWallHeatRate += abs(integratedNormalHeatRate);

                postProcessingFile
                    << "\nPatch name                          : "
                    << p.name()
                    << "\n";

                postProcessingFile
                    << "Boundary-condition type              : "
                    << boundaryType
                    << "\n";

                postProcessingFile
                    << "Number of faces                      : "
                    << p.faces().size()
                    << "\n";

                postProcessingFile
                    << "Patch area [m^2]                     : "
                    << patchArea
                    << "\n";

                const auto averageHeatFluxIterator = averageWallNormalHeatFlux.find(p.name());

                if (averageHeatFluxIterator != averageWallNormalHeatFlux.end())
                {
                    postProcessingFile
                        << "Average wall-normal heat flux [W/m^2]: "
                        << averageHeatFluxIterator->second
                        << "\n";
                }

                const auto averageIrradiationIterator = averageWallIrradiation.find(p.name());

                if (averageIrradiationIterator != averageWallIrradiation.end())
                {
                    postProcessingFile
                        << "Average wall irradiation [W/m^2]    : "
                        << averageIrradiationIterator->second
                        << "\n";
                }

                postProcessingFile
                    << "Integrated wall-normal heat rate [W] : "
                    << integratedNormalHeatRate
                    << "\n";

                postProcessingFile
                    << "Integrated incident irradiation [W]  : "
                    << integratedIncidentIrradiation
                    << "\n";

                postProcessingFile
                    << "Minimum wall-normal flux [W/m^2]     : "
                    << minimumNormalHeatFlux
                    << "\n";

                postProcessingFile
                    << "Maximum wall-normal flux [W/m^2]     : "
                    << maximumNormalHeatFlux
                    << "\n";

                postProcessingFile
                    << "Minimum wall irradiation [W/m^2]     : "
                    << minimumIrradiation
                    << "\n";

                postProcessingFile
                    << "Maximum wall irradiation [W/m^2]     : "
                    << maximumIrradiation
                    << "\n";
            }

            /*
                Calculate the normalised global energy-balance error.
            */
            double relativeGlobalBalanceError = 0.0;

            if (totalAbsoluteWallHeatRate > numeric_limits<double>::epsilon()){
                relativeGlobalBalanceError = abs(globalIntegratedWallHeatRate) / totalAbsoluteWallHeatRate;
            }

            postProcessingFile
                << "\n------------------------------------------------------------\n"
                << "Global radiative balance\n"
                << "------------------------------------------------------------\n";

            postProcessingFile
                << "Sum of integrated wall-normal heat rates [W]: "
                << globalIntegratedWallHeatRate
                << "\n";

            postProcessingFile
                << "Sum of absolute wall heat rates [W]         : "
                << totalAbsoluteWallHeatRate
                << "\n";

            postProcessingFile
                << "Relative global balance error               : "
                << relativeGlobalBalanceError
                << "\n";

            postProcessingFile
                << "\nSign convention:\n"
                << "A positive wall-normal heat flux indicates radiative energy\n"
                << "travelling in the direction of the outward wall normal.\n";

            postProcessingFile
                << "============================================================\n";

            postProcessingFile.close();

            cout
                << "\nDOM simulation converged after "
                << count
                << " outer iterations.\n";

            cout
                << "DOM post-processing data written to:\n"
                << postProcessingFilePath
                << "\n";
        }

        void runDOMPureRadiationParticipatingMedium(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());
            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());

                updateGInternal() ;
                updateTInternal() ;

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsI(*intensityPtr);
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }



                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */
            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;
            updateTInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();

        }

        void runDOMPureRadiationParticipatingMediumConstantT(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());
            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());

                updateGInternal() ;
                //updateTInternal() ;

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsI(*intensityPtr);
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */
            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();

        }

        void runDOMPureRadiationParticipatingMediumConstantTIsotropicScattering(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());
            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsIIsotropicScattering(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());

                updateGInternal() ;
                //updateTInternal() ;

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsIIsotropicScattering(*intensityPtr);
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */
            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTermInX() ;

        }

        void runDOMRadiationConduction(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());

            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;
            double maxFieldChangeT = -1e4 ;
            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            ofstream outInnerT(paths_.getOutputDirectory(meshDimensions_) / "innerResidualT.csv") ;
            outInnerT << "iteration,initial_residual,final_residual,ratio" << endl ;

            ofstream outFinalT(paths_.getOutputDirectory(meshDimensions_) / "finalRedidualT.csv") ;
            outFinalT << "outer_iteration,initial_residual,final_residual,ratio" << endl ;

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i){
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            Diffusion diffusionT(T_, mesh_, Tinit_, k_) ;
            Source sourceT(T_, mesh_, sourceDiagonalT_, sourceTermT_) ;

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);
                ofstream outSimulation(simulationFilePath);

                if (!outSimulation){
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());
                SparseLinearVector Tchange_ = T_.sys().X() ;

                for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_){
                    DOMatrix& intensity_ = *intensityPtr_;
                    for (Patch p : mesh_.boundaryPatches()){
                        intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                    }
                }

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsI(*intensityPtr);
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                updateGInternal() ;
                setSourceVectorsT() ;

                T_.clear() ;
                diffusionT.assemble() ;
                sourceT.assemble() ;
                T_.solve(relativeTolerance, true, outInnerT, outFinalT, count) ;

                Tchange_ = Tchange_ - T_.sys().X() ;
                maxFieldChangeT = Tchange_.maxAbsVal() ;
                maxFieldChangeI = iterationMaximumFieldChange;

                writeVTKRectilinearIntensityFile(outSimulation);
                outSimulation.close();

                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = !(converged && (maxFieldChangeT < changeTolerance)) && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */

            if (!(converged && (maxFieldChangeT < changeTolerance)))
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;

            /*
                Update the final temperature boundary values.
            */

            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();
            heatData_->calculateCenterlineTemperatureAndGInZ() ;
            heatData_->calculateDimensionlessCenterlineTemperatureInZ() ;

        }

        void runDOMRadiationAnnularSourceRectangularEnclosure(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());
            updateAnnularTemperature(0.20, 0.25, 1000) ;

            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();

            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */

            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_(mesh_.numCells());

                //updateGInternal() ;
                //updateTInternal() ;

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsI(*intensityPtr);
                }

                for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_){
                    DOMatrix& intensity_ = *intensityPtr_ ;
                    for (Patch p : mesh_.boundaryPatches()){
                        intensity_.iabf().boundary()[p.index()]->updateBoundaryValues() ;
                    }
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();
                    setSourceVectorsI(IF) ;
                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange >= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTermInX() ;

        }

        void runDOMRadiationHotWallColdMediumRectangularEnclosure(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());

            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();

            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */

            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
                setSourceVectorsI(intensity_);
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging){
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_(mesh_.numCells());

                for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_){
                    DOMatrix& intensity_ = *intensityPtr_ ;
                    for (Patch p : mesh_.boundaryPatches()){
                        intensity_.iabf().boundary()[p.index()]->updateBoundaryValues() ;
                    }
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();
                    setSourceVectorsI(IF) ;
                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange >= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();
            heatData_->calculateRadiativeSourceTerm(a_) ;
            heatData_->calculateCenterlineRadiativeSourceTermInX() ;

        }

        void runDOMRadiationScattering(){

            FieldFileReader TReader(mesh_, T_, flux_, k_, Tinit_, fieldFiles_[0]) ;
            Tinit_ = Field(TReader.internal(), std::move(TReader.boundary())) ;

            paths_.createOutputDirectory(caseName_, subCaseName_);

            T_.clear();
            T_.initializeSolution(Tinit_.internal());
            heatData_ = make_unique<DiscreteOrdinateWallHeatFlux>(mesh_, intensities_, T_, Tinit_, gInternal_, gFace_);

            bool fieldsAreChanging = true;
            bool converged = false;

            int count = 0;
            const int maxIter = 10000;
            double maxFieldChangeI = numeric_limits<double>::max();
            const filesystem::path outputDirectory = paths_.getOutputDirectory(caseName_,subCaseName_);
            const filesystem::path simulationDirectory = outputDirectory / "SimulationFolder";

            filesystem::create_directories(simulationDirectory);

            /*
                Construct the ordinate directions and intensity equations.
                intensities_ is owned by the DOM object and is populated
                directly by OrdinatesConstructor.
            */

            OrdinatesConstructor ordinatesConstructor( addr_, mesh_, paths_, emissivity_, Tinit_, intensities_, ordinate_);

            /*
                Allocate one stored solution vector for every ordinate.
            */
            for (size_t i = 0; i < intensities_.size(); ++i)
            {
                simulationXI_.emplace_back(mesh_.numCells());
            }

            vector<Convection> convectionI_;
            vector<Source> sourceI_;

            convectionI_.reserve(intensities_.size());
            sourceI_.reserve(intensities_.size());

            /*
                Construct the convection and source discretisation objects
                for every ordinate.
            */

            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                convectionI_.emplace_back(intensity_.I(), intensity_.F(), intensity_.iabf());
                sourceI_.emplace_back(intensity_.I(), mesh_, intensity_.sourceDiagonalI(), intensity_.sourceTermI());
                intensity_.openFile(caseName_, subCaseName_);
                ordinateDirections_.push_back(intensity_.S());
            }

            /*
                Outer angular iteration.
            */
            while (fieldsAreChanging)
            {
                const filesystem::path simulationFilePath = simulationDirectory / makeSimulationFileName(count);

                ofstream outSimulation(simulationFilePath);

                if (!outSimulation)
                {
                    throw runtime_error("DOM: unable to create simulation output file: " + simulationFilePath.string());
                }

                bool notConvergedInTotal = false;
                double iterationMaximumFieldChange = 0.0;
                SparseLinearVector Ichange_( mesh_.numCells());

                updateGInternal() ;
                updateTInternal() ;

                for (auto& intensityPtr : intensities_){
                    setSourceVectorsIAnisotropicScatteringSingleDirection(*intensityPtr);
                }

                for (size_t i = 0; i < intensities_.size(); ++i)
                {
                    DOMatrix& IF = *intensities_[i];

                    /*
                        Save the intensity solution from the previous
                        outer iteration.
                    */

                    Ichange_ = IF.I().sys().X();

                    /*
                        Clear and reassemble the current ordinate equation.
                    */
                    IF.I().clear();

                    convectionI_[i].assemble();
                    sourceI_[i].assemble();
                    IF.I().solve(relativeTolerance, true, IF.outInnerI(), IF.outFinalI(), count);

                    /*
                        Store the current intensity solution for VTK output.
                    */
                    simulationXI_[i] = IF.I().sys().X();

                    /*
                        Calculate the maximum change in the current
                        ordinate intensity field.
                    */
                    Ichange_ = Ichange_ - IF.I().sys().X();
                    const double currentFieldChange = Ichange_.maxAbsVal();

                    iterationMaximumFieldChange = max(iterationMaximumFieldChange, currentFieldChange);

                    if (currentFieldChange>= changeTolerance){
                        notConvergedInTotal = true;
                    }
                }

                maxFieldChangeI = iterationMaximumFieldChange;
                writeVTKRectilinearIntensityFile(outSimulation);

                outSimulation.close();
                ++count;
                converged = !notConvergedInTotal;
                fieldsAreChanging = notConvergedInTotal && (count < maxIter);
            }

            /*
                Prevent post-processing of a non-converged solution.
            */
            if (!converged)
            {
                throw runtime_error(
                    "DOM outer iteration failed to converge within "
                    + to_string(maxIter)
                    + " iterations. Final maximum intensity change = "
                    + to_string(maxFieldChangeI)
                );
            }

            /*
                Update the final intensity boundary values before calculating
                wall heat flux and irradiation.

                Patch is intentionally copied here because Patch::index()
                is currently not declared const.
            */
            for (unique_ptr<DOMatrix>& intensityPtr_ : intensities_)
            {
                DOMatrix& intensity_ = *intensityPtr_;
                for (Patch p : mesh_.boundaryPatches()){
                    intensity_.iabf().boundary()[p.index()]->updateBoundaryValues();
                }
            }

            updateGInternal() ;
            updateTInternal() ;

            /*
                Update the final temperature boundary values.
            */
            for (Patch p : mesh_.boundaryPatches()){
                Tinit_.boundary()[p.index()]->updateBoundaryValues();
            }

            heatData_->calculateHeatFlux();
            heatData_->calculateIrradiation();
            heatData_->calculateAverages();
            heatData_->calculateAvgTempartures();
            heatData_->calculateDimensionlessHeatFlux();
            heatData_->calculateDimensionlessHeatFlux2();
            heatData_->calculateDimensionlessSlipTemperature();
            heatData_->calculateCenterlineTemperatureAndGInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX();
            heatData_->calculateDimensionlessCenterlineTemperatureInX2();
            heatData_->calculateRadiativeSourceTerm(a_);
            heatData_->calculateCenterlineRadiativeSourceTermInX() ;

        }



        // Utility functions

        string makeSimulationFileName(int step){
            std::ostringstream name;

            name << "simulation"
                << std::setw(6) << std::setfill('0') << step
                << ".vtk";

            return name.str();
        }

        void writeVTKRectilinearIntensityFile(ostream& out){

            if(!out){
                throw runtime_error("Could not open VTK file for writing") ;
            }

            const int numActiveCells = mesh_.numCells() ;
            const int numX = numCellsXYZ_[0] ;
            const int numY = numCellsXYZ_[1] ;
            const int numZ = numCellsXYZ_[2] ;
            const int numStructuredCells = numX*numY*numZ ;

            if(numX<=0 || numY<=0 || numZ<=0){
                throw runtime_error("Invalid structured mesh dimensions while writing DOM VTK file") ;
            }

            if(simulationXI_.size()!=ordinateDirections_.size()){
                throw runtime_error("The number of stored DOM intensity fields does not match the number of ordinate directions") ;
            }

            for(size_t ord=0; ord<simulationXI_.size(); ++ord){
                if(simulationXI_[ord].size()!=numActiveCells){
                    throw runtime_error(
                        "Stored DOM intensity field size does not match the number of active mesh cells for ordinate "
                        + to_string(ord)
                    ) ;
                }
            }

            const double dx = dimensionsXYZ_[0]/static_cast<double>(numX) ;
            const double dy = dimensionsXYZ_[1]/static_cast<double>(numY) ;
            const double dz = dimensionsXYZ_[2]/static_cast<double>(numZ) ;

            /*
                Map every active medium cell to its original Cartesian
                position. Deleted cells inside the pipes remain hidden
                cells in the VTK rectilinear grid.
            */

            vector<int> activeCellToStructuredCell(numActiveCells,-1) ;
            vector<int> structuredCellToActiveCell(numStructuredCells,-1) ;
            vector<unsigned char> ghostCell(numStructuredCells,32) ;

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                const array<double,3> centre = mesh_.cellCentre(cellIndex) ;

                int i = static_cast<int>(round(centre[0]/dx-0.5)) ;
                int j = static_cast<int>(round(centre[1]/dy-0.5)) ;
                int k = static_cast<int>(round(centre[2]/dz-0.5)) ;

                if(i<0 || i>=numX || j<0 || j>=numY || k<0 || k>=numZ){
                    throw runtime_error(
                        "Active mesh cell could not be mapped to the Cartesian DOM VTK lattice. Cell index: "
                        + to_string(cellIndex)
                    ) ;
                }

                const int structuredCellIndex = k*numX*numY+j*numX+i ;

                if(structuredCellToActiveCell[structuredCellIndex]!=-1){
                    throw runtime_error(
                        "More than one active mesh cell maps to the same Cartesian DOM VTK cell. Structured cell index: "
                        + to_string(structuredCellIndex)
                    ) ;
                }

                activeCellToStructuredCell[cellIndex] = structuredCellIndex ;
                structuredCellToActiveCell[structuredCellIndex] = cellIndex ;
                ghostCell[structuredCellIndex] = 0 ;
            }

            /*
                Calculate DOM diagnostics only for the active medium cells.
            */

            vector<double> activeIrradiation(numActiveCells,0.0) ;
            vector<double> activeHeatFluxX(numActiveCells,0.0) ;
            vector<double> activeHeatFluxY(numActiveCells,0.0) ;
            vector<double> activeHeatFluxZ(numActiveCells,0.0) ;
            vector<double> activeHeatFluxMagnitude(numActiveCells,0.0) ;
            vector<double> activeRadiativeSourceTerm(numActiveCells,0.0) ;

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                for(size_t ord=0; ord<simulationXI_.size(); ++ord){

                    const double intensity = simulationXI_[ord](cellIndex) ;
                    const OrdinateS& S = ordinateDirections_[ord] ;

                    activeIrradiation[cellIndex] += intensity*S.w ;

                    activeHeatFluxX[cellIndex] += intensity*S.x*S.w ;
                    activeHeatFluxY[cellIndex] += intensity*S.y*S.w ;
                    activeHeatFluxZ[cellIndex] += intensity*S.z*S.w ;
                }

                activeHeatFluxMagnitude[cellIndex] = sqrt(
                    activeHeatFluxX[cellIndex]*activeHeatFluxX[cellIndex]
                    +
                    activeHeatFluxY[cellIndex]*activeHeatFluxY[cellIndex]
                    +
                    activeHeatFluxZ[cellIndex]*activeHeatFluxZ[cellIndex]
                ) ;

                activeRadiativeSourceTerm[cellIndex] = a_(cellIndex)*(
                    4.0*sigma*pow(T_.sys().X()(cellIndex),4)
                    -
                    activeIrradiation[cellIndex]
                ) ;
            }

            /*
                Transfer all active-cell values to the complete structured
                lattice. Deleted pipe cells remain zero and are hidden using
                vtkGhostType.
            */

            vector<double> structuredIrradiation(numStructuredCells,0.0) ;
            vector<double> structuredTemperature(numStructuredCells,0.0) ;
            vector<double> structuredHeatFluxX(numStructuredCells,0.0) ;
            vector<double> structuredHeatFluxY(numStructuredCells,0.0) ;
            vector<double> structuredHeatFluxZ(numStructuredCells,0.0) ;
            vector<double> structuredHeatFluxMagnitude(numStructuredCells,0.0) ;
            vector<double> structuredRadiativeSourceTerm(numStructuredCells,0.0) ;

            vector<vector<array<double,3>>> structuredOrdinateIntensity(
                simulationXI_.size(),
                vector<array<double,3>>(
                    numStructuredCells,
                    {0.0,0.0,0.0}
                )
            ) ;

            for(int cellIndex=0; cellIndex<numActiveCells; ++cellIndex){

                const int structuredCellIndex =
                    activeCellToStructuredCell[cellIndex] ;

                structuredIrradiation[structuredCellIndex] =
                    activeIrradiation[cellIndex] ;

                structuredTemperature[structuredCellIndex] =
                    T_.sys().X()(cellIndex) ;

                structuredHeatFluxX[structuredCellIndex] =
                    activeHeatFluxX[cellIndex] ;

                structuredHeatFluxY[structuredCellIndex] =
                    activeHeatFluxY[cellIndex] ;

                structuredHeatFluxZ[structuredCellIndex] =
                    activeHeatFluxZ[cellIndex] ;

                structuredHeatFluxMagnitude[structuredCellIndex] =
                    activeHeatFluxMagnitude[cellIndex] ;

                structuredRadiativeSourceTerm[structuredCellIndex] =
                    activeRadiativeSourceTerm[cellIndex] ;

                for(size_t ord=0; ord<simulationXI_.size(); ++ord){

                    const double intensity =
                        simulationXI_[ord](cellIndex) ;

                    const OrdinateS& S =
                        ordinateDirections_[ord] ;

                    structuredOrdinateIntensity[ord][structuredCellIndex] = {
                        intensity*S.x,
                        intensity*S.y,
                        intensity*S.z
                    } ;
                }
            }

            out << "# vtk DataFile Version 3.0" << endl ;
            out << "DOM heat exchanger radiation result" << endl ;
            out << "ASCII" << endl ;
            out << "DATASET RECTILINEAR_GRID" << endl ;

            out << "DIMENSIONS "
                << numX+1 << " "
                << numY+1 << " "
                << numZ+1 << endl << endl ;

            out << "X_COORDINATES " << numX+1 << " double" << endl ;

            for(int i=0; i<=numX; ++i){
                out << setprecision(12)
                    << dx*static_cast<double>(i)
                    << (i<numX ? " " : "\n") ;
            }

            out << endl ;

            out << "Y_COORDINATES " << numY+1 << " double" << endl ;

            for(int j=0; j<=numY; ++j){
                out << setprecision(12)
                    << dy*static_cast<double>(j)
                    << (j<numY ? " " : "\n") ;
            }

            out << endl ;

            out << "Z_COORDINATES " << numZ+1 << " double" << endl ;

            for(int k=0; k<=numZ; ++k){
                out << setprecision(12)
                    << dz*static_cast<double>(k)
                    << (k<numZ ? " " : "\n") ;
            }

            out << endl ;

            /*
                The complete structured cell count must be written here.
                Active cells contain the DOM solution. Deleted cells inside
                the pipes are represented as hidden VTK cells.
            */

            out << "CELL_DATA " << numStructuredCells << endl << endl ;

            out << "SCALARS vtkGhostType unsigned_char 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){
                out << static_cast<int>(ghostCell[structuredCellIndex]) << endl ;
            }

            out << endl ;

            out << "SCALARS irradiation_G double 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){
                out << setprecision(12)
                    << structuredIrradiation[structuredCellIndex] << endl ;
            }

            out << endl ;

            out << "SCALARS temperature_T double 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){
                out << setprecision(12)
                    << structuredTemperature[structuredCellIndex] << endl ;
            }

            out << endl ;

            out << "VECTORS radiativeHeatFlux double" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){

                out << setprecision(12)
                    << structuredHeatFluxX[structuredCellIndex] << " "
                    << structuredHeatFluxY[structuredCellIndex] << " "
                    << structuredHeatFluxZ[structuredCellIndex] << endl ;
            }

            out << endl ;

            out << "SCALARS radiativeHeatFluxMagnitude double 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){
                out << setprecision(12)
                    << structuredHeatFluxMagnitude[structuredCellIndex] << endl ;
            }

            out << endl ;

            out << "SCALARS radiativeSourceTerm double 1" << endl ;
            out << "LOOKUP_TABLE default" << endl ;

            for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){
                out << setprecision(12)
                    << structuredRadiativeSourceTerm[structuredCellIndex] << endl ;
            }

            out << endl ;

            /*
                Store every ordinate intensity as a vector aligned with its
                discrete ordinate direction.
            */

            for(size_t ord=0; ord<simulationXI_.size(); ++ord){

                out << "VECTORS intensity_ordinate_" << ord << " double" << endl ;

                for(int structuredCellIndex=0; structuredCellIndex<numStructuredCells; ++structuredCellIndex){

                    out << setprecision(12)
                        << structuredOrdinateIntensity[ord][structuredCellIndex][0] << " "
                        << structuredOrdinateIntensity[ord][structuredCellIndex][1] << " "
                        << structuredOrdinateIntensity[ord][structuredCellIndex][2] << endl ;
                }

                out << endl ;
            }
        }


        void setSourceVectorsI(DOMatrix& doMatrix){
            
            for(int i=0; i < mesh_.numCells(); ++i){

                double S ;
                if(a_(i)==0 && sigmaS_(i)==0){
                    S = 0 ;
                }
                else{
                    S = (1 - (sigmaS_(i)/(a_(i)+sigmaS_(i))))*(sigma*pow(T_.sys().X()(i),4)/pi) ;
                }
                doMatrix.sourceDiagonalI()(i) = (a_(i)+sigmaS_(i)) ;
                doMatrix.sourceTermI()(i) = (a_(i)+sigmaS_(i))*S ;
            }
        }

        void setSourceVectorsIIsotropicScattering(DOMatrix& doMatrix){
            
            for(int i=0; i < mesh_.numCells(); ++i){
                double S ;
                double omega = (sigmaS_(i)/(a_(i)+sigmaS_(i))) ;
                S = (1 - omega)*(sigma*pow(T_.sys().X()(i),4)/pi) + (omega/(4*pi))*gInternal_(i) ;
                doMatrix.sourceDiagonalI()(i) = (a_(i)+sigmaS_(i)) ;
                doMatrix.sourceTermI()(i) = (a_(i)+sigmaS_(i))*S ;
            }
        }

        void setSourceVectorsIAnisotropicScatteringSingleDirection(DOMatrix& doMatrix){
            
            for(int i=0; i < mesh_.numCells(); ++i){
                double S ;
                double omega = (sigmaS_(i)/(a_(i)+sigmaS_(i))) ;
                S = (1 - omega)*(sigma*pow(T_.sys().X()(i),4)/pi) + (omega/(4*pi))*gInternal_(i) ;

                for(unique_ptr<DOMatrix>& intensityPtr_ : intensities_){
                    DOMatrix& intensity = *intensityPtr_;
                    S += (omega/(4*pi))*(anisotropyFactor_*(doMatrix.S().x*intensity.S().x)*intensity.S().w*intensity.I().sys().X()(i)) ;
                }

                doMatrix.sourceDiagonalI()(i) = (a_(i)+sigmaS_(i)) ;
                doMatrix.sourceTermI()(i) = (a_(i)+sigmaS_(i))*S ;
            }
        }

        // used to get the irradiation from the intensity field 
        void updateGInternal(){
            gInternal_.clear() ;
            for(int i = 0; i < mesh_.numCells(); ++i){
                for(unique_ptr<DOMatrix>& intensity: intensities_){
                    DOMatrix& IF = *intensity ;
                    gInternal_(i) += IF.I().sys().X()(i)*IF.S().w ; 
                }
            }
        }

        void updateTInternal(){
            for(int i=0; i<mesh_.numCells(); ++i){
                T_.sys().X()(i) = powl((gInternal_(i)/(4*sigma)), 0.25) ;
            }
        }

        void setSourceVectorsT(){
            for(int i=0; i<mesh_.numCells(); ++i){
                sourceDiagonalT_(i) = 0.0 ;
                sourceTermT_(i) = -a_(i)*(4*sigma*pow(T_.sys().X()(i),4) - gInternal_(i)) ;
            }
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

        
    
};

