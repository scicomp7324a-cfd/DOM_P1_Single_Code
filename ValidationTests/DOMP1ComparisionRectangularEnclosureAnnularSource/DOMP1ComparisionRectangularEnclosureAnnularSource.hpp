# pragma once

# include "PathConfig.hpp"

# include <iostream>
# include <vector>
# include <map>
# include <filesystem>
# include <string>
# include <fstream>
# include <cstdlib>
# include <stdexcept>
# include <cmath>
# include <chrono>
# include <iterator>

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
# include "P1.hpp"

# include "DOMatrix.hpp"
# include "OrdinateField.hpp"
# include "OrdinatesConstructor.hpp"
# include "PeriodicBC.hpp"
# include "DiscreteOrdinateBC.hpp"
# include "DiscreteOrdinateWallHeatFlux.hpp"
# include "DOM.hpp"

using namespace std;

class DOMP1ComparisionRectangularEnclosureAnnularSource{

    public:

    DOMP1ComparisionRectangularEnclosureAnnularSource(){}

    void test(){

        P1RadiativeEnergyAndFluxProfileTest() ;
        DOMRadiativeEnergyAndFluxProfileTest(2) ;
        DOMRadiativeEnergyAndFluxProfileTest(4) ;
        DOMRadiativeEnergyAndFluxProfileTest(6) ;
        DOMRadiativeEnergyAndFluxProfileTest(8) ;
        P1MeshConvergenceTest() ;
        DOMMeshConvergenceTest(2) ;
        DOMMeshConvergenceTest(4) ;
        DOMMeshConvergenceTest(6) ;
        DOMMeshConvergenceTest(8) ;
        DOMAngularConvergenceTest() ;

    }

    void P1RadiativeEnergyAndFluxProfileTest(){
                
        // Creates and validates all required input/output paths before starting the simulations.
        PathConfig paths ;
            
        // Stops execution if any required path or output directory setup fails.
        try{
            paths.validate() ;
            paths.ensureOutputDir() ;
        }
        catch(const exception& e){
            cerr << "Path error: " << e.what() << endl ;
            throw runtime_error(e.what()) ;
        }

        string parentFolder = "DOMP1ComparisionRectangularEnclosureAnnularSource" ;
        string caseName = "P1RadiativeEnergyAndFluxProfileTest" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
                
        // Lists the mesh cases used in the study.
        vector<string> subCaseNames = {"320x320x1"} ;
            
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 100.0 ;

        // Non-scattering medium.
        double sigmaS = 0.0 ;

        // Zero conductivity.
        double k = 0 ;

        // Cold black walls in x and y and periodic boundaries in z.
        vector<double> emissivity = {1.0,1.0,1.0,1.0,0.0,0.0} ; 

        // Physical configuration.
        double domainLength = 1.0 ;
        double backgroundTemperature = 0.0 ;
        double sourceTemperature = 1000.0 ;
        double innerRadius = 0.20 ;
        double outerRadius = 0.25 ;

        for(string subCaseName: subCaseNames){

            paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the field files for the initial G and T data.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the P1 solver for the current case and executes the full simulation.
            string p1CaseName = parentFolder + "/" + caseName ;
            P1 p1(meshFiles, fieldFiles, absorptionCoefficient, sigmaS, k, emissivity, subCaseName, p1CaseName, subCaseName) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            p1.runP1RadiationAnnularSourceRectangularEnclosure() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            p1.heatData()->calculateCenterlineTemperatureAndG() ;
            p1.heatData()->calculateCenterlineFlux() ;
            p1.heatData()->calculateCenterlineRadiativeSourceTerm() ;

            auto firstCenterlineCoordinate = p1.heatData()->centerlineG().begin() ;
            auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
            double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
            double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

            ofstream centerline_G(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineG.csv") ;
            ofstream centerline_flux(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineFlux.csv") ;
            ofstream centerline_source(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineRadiativeSourceTerm.csv") ;
            ofstream symmetry_info(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "symmetry_checks.csv") ;

            centerline_G << "background_temperature," << backgroundTemperature << endl ;
            centerline_G << "source_temperature," << sourceTemperature << endl ;
            centerline_G << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_G << "scattering_coefficient," << sigmaS << endl ;
            centerline_G << "domain_length," << domainLength << endl ;
            centerline_G << "inner_radius," << innerRadius << endl ;
            centerline_G << "outer_radius," << outerRadius << endl ;
            centerline_G << "runtime," << elapsed.count() << endl ;
            centerline_G << "sub_case," << subCaseName << endl << endl ;
            centerline_G << "x,x_over_L,G_P1" << endl ;

            centerline_flux << "background_temperature," << backgroundTemperature << endl ;
            centerline_flux << "source_temperature," << sourceTemperature << endl ;
            centerline_flux << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_flux << "scattering_coefficient," << sigmaS << endl ;
            centerline_flux << "domain_length," << domainLength << endl ;
            centerline_flux << "inner_radius," << innerRadius << endl ;
            centerline_flux << "outer_radius," << outerRadius << endl ;
            centerline_flux << "runtime," << elapsed.count() << endl ;
            centerline_flux << "sub_case," << subCaseName << endl << endl ;
            centerline_flux << "x,x_over_L,radiative_flux_x_P1,radiative_flux_magnitude_P1" << endl ;

            centerline_source << "background_temperature," << backgroundTemperature << endl ;
            centerline_source << "source_temperature," << sourceTemperature << endl ;
            centerline_source << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_source << "scattering_coefficient," << sigmaS << endl ;
            centerline_source << "domain_length," << domainLength << endl ;
            centerline_source << "inner_radius," << innerRadius << endl ;
            centerline_source << "outer_radius," << outerRadius << endl ;
            centerline_source << "runtime," << elapsed.count() << endl ;
            centerline_source << "sub_case," << subCaseName << endl << endl ;
            centerline_source << "x,x_over_L,radiative_source_term_P1" << endl ;

            for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << (key-westWallLocation)/domainLength << "," << val << endl ; }
            for(auto& [key,val]: p1.heatData()->centerlineFlux()){ centerline_flux << key << "," << (key-westWallLocation)/domainLength << "," << val << "," << abs(val) << endl ; }
            for(auto& [key,val]: p1.heatData()->centerlineRadiativeSourceTerm()){ centerline_source << key << "," << (key-westWallLocation)/domainLength << "," << val << endl ; }

            vector<double> centerlineGValues ;
            vector<double> centerlineFluxValues ;
            vector<double> centerlineSourceValues ;

            for(auto& [key,val]: p1.heatData()->centerlineG()){ centerlineGValues.push_back(val) ; }
            for(auto& [key,val]: p1.heatData()->centerlineFlux()){ centerlineFluxValues.push_back(val) ; }
            for(auto& [key,val]: p1.heatData()->centerlineRadiativeSourceTerm()){ centerlineSourceValues.push_back(val) ; }

            double maximumG = 0 ;
            double maximumFlux = 0 ;
            double maximumSource = 0 ;

            double maximumGSymmetryError = 0 ;
            double maximumFluxAntisymmetryError = 0 ;
            double maximumFluxMagnitudeSymmetryError = 0 ;
            double maximumSourceSymmetryError = 0 ;

            for(size_t i=0; i<centerlineGValues.size(); ++i){

                size_t mirrorIndex = centerlineGValues.size()-1-i ;

                maximumG = max(maximumG,abs(centerlineGValues[i])) ;
                maximumFlux = max(maximumFlux,abs(centerlineFluxValues[i])) ;
                maximumSource = max(maximumSource,abs(centerlineSourceValues[i])) ;

                maximumGSymmetryError = max(maximumGSymmetryError,abs(centerlineGValues[i]-centerlineGValues[mirrorIndex])) ;
                maximumFluxAntisymmetryError = max(maximumFluxAntisymmetryError,abs(centerlineFluxValues[i]+centerlineFluxValues[mirrorIndex])) ;
                maximumFluxMagnitudeSymmetryError = max(maximumFluxMagnitudeSymmetryError,abs(abs(centerlineFluxValues[i])-abs(centerlineFluxValues[mirrorIndex]))) ;
                maximumSourceSymmetryError = max(maximumSourceSymmetryError,abs(centerlineSourceValues[i]-centerlineSourceValues[mirrorIndex])) ;
            }

            double relativeGSymmetryError = maximumG > numeric_limits<double>::epsilon() ? maximumGSymmetryError/maximumG : 0 ;
            double relativeFluxAntisymmetryError = maximumFlux > numeric_limits<double>::epsilon() ? maximumFluxAntisymmetryError/maximumFlux : 0 ;
            double relativeFluxMagnitudeSymmetryError = maximumFlux > numeric_limits<double>::epsilon() ? maximumFluxMagnitudeSymmetryError/maximumFlux : 0 ;
            double relativeSourceSymmetryError = maximumSource > numeric_limits<double>::epsilon() ? maximumSourceSymmetryError/maximumSource : 0 ;

            symmetry_info << "quantity,symmetry_type,maximum_absolute_error,maximum_relative_error" << endl ;
            symmetry_info << "G,even," << maximumGSymmetryError << "," << relativeGSymmetryError << endl ;
            symmetry_info << "radiative_flux_x,odd," << maximumFluxAntisymmetryError << "," << relativeFluxAntisymmetryError << endl ;
            symmetry_info << "radiative_flux_magnitude,even," << maximumFluxMagnitudeSymmetryError << "," << relativeFluxMagnitudeSymmetryError << endl ;
            symmetry_info << "radiative_source_term,even," << maximumSourceSymmetryError << "," << relativeSourceSymmetryError << endl ;
                    
            centerline_G << endl ;
            centerline_flux << endl ;
            centerline_source << endl ;
        }
    }

    void DOMRadiativeEnergyAndFluxProfileTest(int ordinate){
                
        // Creates and validates all required input/output paths before starting the simulations.
        PathConfig paths ;
            
        // Stops execution if any required path or output directory setup fails.
        try{
            paths.validate() ;
            paths.ensureOutputDir() ;
        }
        catch(const exception& e){
            cerr << "Path error: " << e.what() << endl ;
            throw runtime_error(e.what()) ;
        }

        string parentFolder = "DOMP1ComparisionRectangularEnclosureAnnularSource" ;
        string caseName = "DOMRadiativeEnergyAndFluxProfileTest" ;
        string outputCaseName = caseName + "/ordinate" + to_string(ordinate) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
                
        // Lists the mesh cases used in the study.
        vector<string> subCaseNames = {"320x320x1"} ;
            
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 100.0 ;

        // Non-scattering medium.
        double sigmaS = 0.0 ;

        // Zero conductivity.
        double k = 0 ;

        // Cold black walls in x and y and periodic boundaries in z.
        vector<double> emissivity = {1.0,1.0,1.0,1.0,0.0,0.0} ; 

        // Physical configuration.
        double domainLength = 1.0 ;
        double backgroundTemperature = 0.0 ;
        double sourceTemperature = 1000.0 ;
        double innerRadius = 0.20 ;
        double outerRadius = 0.25 ;

        paths.createOutputDirectory(parentFolder, outputCaseName, "") ;

        for(string subCaseName: subCaseNames){

            paths.createOutputDirectory(parentFolder, outputCaseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the field files for the initial G and T data.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the DOM solver for the current case and executes the full simulation.
            string domCaseName = parentFolder + "/" + outputCaseName ;
            DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, subCaseName, domCaseName, subCaseName, ordinate) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            dom.runDOMRadiationAnnularSourceRectangularEnclosure() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            dom.heatData()->calculateCenterlineTemperatureAndGInX() ;
            dom.heatData()->calculateCenterlineFluxInX() ;
            dom.heatData()->calculateCenterlineRadiativeSourceTermInX() ;

            auto firstCenterlineCoordinate = dom.heatData()->centerlineGInX().begin() ;
            auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
            double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
            double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

            ofstream centerline_G(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "centerlineG.csv") ;
            ofstream centerline_flux(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "centerlineFlux.csv") ;
            ofstream centerline_source(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "centerlineRadiativeSourceTerm.csv") ;
            ofstream symmetry_info(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "symmetry_checks.csv") ;

            centerline_G << "background_temperature," << backgroundTemperature << endl ;
            centerline_G << "source_temperature," << sourceTemperature << endl ;
            centerline_G << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_G << "scattering_coefficient," << sigmaS << endl ;
            centerline_G << "domain_length," << domainLength << endl ;
            centerline_G << "inner_radius," << innerRadius << endl ;
            centerline_G << "outer_radius," << outerRadius << endl ;
            centerline_G << "ordinate_order," << ordinate << endl ;
            centerline_G << "runtime," << elapsed.count() << endl ;
            centerline_G << "sub_case," << subCaseName << endl << endl ;
            centerline_G << "x,x_over_L,G_DOM" << endl ;

            centerline_flux << "background_temperature," << backgroundTemperature << endl ;
            centerline_flux << "source_temperature," << sourceTemperature << endl ;
            centerline_flux << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_flux << "scattering_coefficient," << sigmaS << endl ;
            centerline_flux << "domain_length," << domainLength << endl ;
            centerline_flux << "inner_radius," << innerRadius << endl ;
            centerline_flux << "outer_radius," << outerRadius << endl ;
            centerline_flux << "ordinate_order," << ordinate << endl ;
            centerline_flux << "runtime," << elapsed.count() << endl ;
            centerline_flux << "sub_case," << subCaseName << endl << endl ;
            centerline_flux << "x,x_over_L,radiative_flux_x_DOM,radiative_flux_magnitude_DOM" << endl ;

            centerline_source << "background_temperature," << backgroundTemperature << endl ;
            centerline_source << "source_temperature," << sourceTemperature << endl ;
            centerline_source << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_source << "scattering_coefficient," << sigmaS << endl ;
            centerline_source << "domain_length," << domainLength << endl ;
            centerline_source << "inner_radius," << innerRadius << endl ;
            centerline_source << "outer_radius," << outerRadius << endl ;
            centerline_source << "ordinate_order," << ordinate << endl ;
            centerline_source << "runtime," << elapsed.count() << endl ;
            centerline_source << "sub_case," << subCaseName << endl << endl ;
            centerline_source << "x,x_over_L,radiative_source_term_DOM" << endl ;

            for(auto& [key,val]: dom.heatData()->centerlineGInX()){ centerline_G << key << "," << (key-westWallLocation)/domainLength << "," << val << endl ; }
            for(auto& [key,val]: dom.heatData()->centerlineFluxInX()){ centerline_flux << key << "," << (key-westWallLocation)/domainLength << "," << val << "," << abs(val) << endl ; }
            for(auto& [key,val]: dom.heatData()->centerlineRadiativeSourceTermInX()){ centerline_source << key << "," << (key-westWallLocation)/domainLength << "," << val << endl ; }

            vector<double> centerlineGValues ;
            vector<double> centerlineFluxValues ;
            vector<double> centerlineSourceValues ;

            for(auto& [key,val]: dom.heatData()->centerlineGInX()){ centerlineGValues.push_back(val) ; }
            for(auto& [key,val]: dom.heatData()->centerlineFluxInX()){ centerlineFluxValues.push_back(val) ; }
            for(auto& [key,val]: dom.heatData()->centerlineRadiativeSourceTermInX()){ centerlineSourceValues.push_back(val) ; }

            double maximumG = 0 ;
            double maximumFlux = 0 ;
            double maximumSource = 0 ;

            double maximumGSymmetryError = 0 ;
            double maximumFluxAntisymmetryError = 0 ;
            double maximumFluxMagnitudeSymmetryError = 0 ;
            double maximumSourceSymmetryError = 0 ;

            for(size_t i=0; i<centerlineGValues.size(); ++i){

                size_t mirrorIndex = centerlineGValues.size()-1-i ;

                maximumG = max(maximumG,abs(centerlineGValues[i])) ;
                maximumFlux = max(maximumFlux,abs(centerlineFluxValues[i])) ;
                maximumSource = max(maximumSource,abs(centerlineSourceValues[i])) ;

                maximumGSymmetryError = max(maximumGSymmetryError,abs(centerlineGValues[i]-centerlineGValues[mirrorIndex])) ;
                maximumFluxAntisymmetryError = max(maximumFluxAntisymmetryError,abs(centerlineFluxValues[i]+centerlineFluxValues[mirrorIndex])) ;
                maximumFluxMagnitudeSymmetryError = max(maximumFluxMagnitudeSymmetryError,abs(abs(centerlineFluxValues[i])-abs(centerlineFluxValues[mirrorIndex]))) ;
                maximumSourceSymmetryError = max(maximumSourceSymmetryError,abs(centerlineSourceValues[i]-centerlineSourceValues[mirrorIndex])) ;
            }

            double relativeGSymmetryError = maximumG > numeric_limits<double>::epsilon() ? maximumGSymmetryError/maximumG : 0 ;
            double relativeFluxAntisymmetryError = maximumFlux > numeric_limits<double>::epsilon() ? maximumFluxAntisymmetryError/maximumFlux : 0 ;
            double relativeFluxMagnitudeSymmetryError = maximumFlux > numeric_limits<double>::epsilon() ? maximumFluxMagnitudeSymmetryError/maximumFlux : 0 ;
            double relativeSourceSymmetryError = maximumSource > numeric_limits<double>::epsilon() ? maximumSourceSymmetryError/maximumSource : 0 ;

            symmetry_info << "quantity,symmetry_type,maximum_absolute_error,maximum_relative_error" << endl ;
            symmetry_info << "G,even," << maximumGSymmetryError << "," << relativeGSymmetryError << endl ;
            symmetry_info << "radiative_flux_x,odd," << maximumFluxAntisymmetryError << "," << relativeFluxAntisymmetryError << endl ;
            symmetry_info << "radiative_flux_magnitude,even," << maximumFluxMagnitudeSymmetryError << "," << relativeFluxMagnitudeSymmetryError << endl ;
            symmetry_info << "radiative_source_term,even," << maximumSourceSymmetryError << "," << relativeSourceSymmetryError << endl ;
                    
            centerline_G << endl ;
            centerline_flux << endl ;
            centerline_source << endl ;
        }
    }

    void P1MeshConvergenceTest(){
                
        // Creates and validates all required input/output paths before starting the simulations.
        PathConfig paths ;
            
        // Stops execution if any required path or output directory setup fails.
        try{
            paths.validate() ;
            paths.ensureOutputDir() ;
        }
        catch(const exception& e){
            cerr << "Path error: " << e.what() << endl ;
            throw runtime_error(e.what()) ;
        }

        string parentFolder = "DOMP1ComparisionRectangularEnclosureAnnularSource" ;
        string caseName = "P1MeshConvergenceTest" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
                
        // Lists the square meshes used in the spatial convergence study.
        vector<int> meshSizes = {20,40,80,160,320} ;

        // Physical dimensions of the square enclosure.
        double domainLength = 1.0 ;
            
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 100.0 ;

        // Non-scattering medium.
        double sigmaS = 0.0 ;

        // No conduction.
        double k = 0.0 ;

        // Cold black physical walls in x and y, with periodic boundaries in z.
        vector<double> emissivity = {1.0,1.0,1.0,1.0,0.0,0.0} ;

        // Stores the spatial convergence data.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, caseName) / "mesh_convergence.csv") ;
        convergence_info << "cells_x,"
                        << "cells_y,"
                        << "total_cells_xy,"
                        << "cell_size,"
                        << "optical_cell_thickness,"
                        << "successive_grid_relative_L2_change,"
                        << "observed_order,"
                        << "runtime" << endl ;

        double previousError = 0 ;
        double previousCellSize = 0 ;
        map<double,double> previousCenterlineG ;
        bool firstMesh = true ;

        for(int meshSize: meshSizes){

            string subCaseName = to_string(meshSize) + "x" + to_string(meshSize) + "x1" ;

            paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the prescribed temperature and initial irradiation fields.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs and runs the P1 annular-source simulation.
            string p1CaseName = parentFolder + "/" + caseName ;
            P1 p1(meshFiles, fieldFiles, absorptionCoefficient, sigmaS, k, emissivity, subCaseName, p1CaseName, subCaseName) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            p1.runP1RadiationAnnularSourceRectangularEnclosure() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            p1.heatData()->calculateCenterlineTemperatureAndG() ;

            double cellSize = domainLength/meshSize ;
            double opticalCellThickness = absorptionCoefficient*cellSize ;

            convergence_info << meshSize << ","
                            << meshSize << ","
                            << meshSize*meshSize << ","
                            << cellSize << ","
                            << opticalCellThickness << "," ;

            if(firstMesh){
                convergence_info << ",," ;
                firstMesh = false ;
            }
            else{
                double relativeL2Error = calculateCenterlineRelativeL2Error(previousCenterlineG,p1.heatData()->centerlineG()) ;
                convergence_info << relativeL2Error << "," ;

                if(previousError == 0){
                    convergence_info << "," ;
                }
                else{
                    double observedOrder = log(previousError/relativeL2Error)/log(previousCellSize/cellSize) ;
                    convergence_info << observedOrder << "," ;
                }

                previousError = relativeL2Error ;
            }

            convergence_info << elapsed.count() << endl ;

            previousCenterlineG = p1.heatData()->centerlineG() ;
            previousCellSize = cellSize ;
        }
    }

    void DOMMeshConvergenceTest(int ordinate){
                
        // Creates and validates all required input/output paths before starting the simulations.
        PathConfig paths ;
            
        // Stops execution if any required path or output directory setup fails.
        try{
            paths.validate() ;
            paths.ensureOutputDir() ;
        }
        catch(const exception& e){
            cerr << "Path error: " << e.what() << endl ;
            throw runtime_error(e.what()) ;
        }

        string parentFolder = "DOMP1ComparisionRectangularEnclosureAnnularSource" ;
        string caseName = "DOMMeshConvergenceTest" ;
        string outputCaseName = caseName + "/ordinate" + to_string(ordinate) ;

        paths.createOutputDirectory(parentFolder, outputCaseName, "") ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
                
        // Lists the square mesh sizes used in the convergence study.
        vector<int> meshSizes = {20,40,80,160,320} ;

        // Physical length of the square enclosure.
        double domainLength = 1.0 ;
            
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 100.0 ;

        // Non-scattering medium.
        double sigmaS = 0.0 ;

        // Zero conductivity.
        double k = 0.0 ;

        // Cold black walls in x and y and periodic boundaries in z.
        vector<double> emissivity = {1.0,1.0,1.0,1.0,0.0,0.0} ; 

        // Stores the convergence data in a directly plottable format.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, outputCaseName, "") / "mesh_convergence.csv") ;
        convergence_info << "ordinate_order,"
                        << "cells_x,"
                        << "cells_y,"
                        << "total_cells_xy,"
                        << "cell_size,"
                        << "optical_cell_thickness,"
                        << "successive_grid_relative_L2_change,"
                        << "observed_order,"
                        << "runtime" << endl ;

        double previousError = 0 ;
        double previousCellSize = 0 ;
        map<double,double> previousCenterlineG ;
        bool firstMesh = true ;

        for(int meshSize: meshSizes){

            string subCaseName = to_string(meshSize) + "x" + to_string(meshSize) + "x1" ;

            paths.createOutputDirectory(parentFolder, outputCaseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the prescribed temperature and initial irradiation fields.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the DOM solver for the current case and executes the annular-source simulation.
            string domCaseName = parentFolder + "/" + outputCaseName ;
            DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, subCaseName, domCaseName, subCaseName, ordinate) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            dom.runDOMRadiationAnnularSourceRectangularEnclosure() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            double cellSize = domainLength/meshSize ;
            double opticalCellThickness = absorptionCoefficient*cellSize ;

            convergence_info << ordinate << ","
                            << meshSize << ","
                            << meshSize << ","
                            << meshSize*meshSize << ","
                            << cellSize << ","
                            << opticalCellThickness << "," ;

            if(firstMesh){
                convergence_info << ",," ;
                firstMesh = false ;
            }
            else{
                double relativeL2Error = calculateCenterlineRelativeL2Error(previousCenterlineG,dom.heatData()->centerlineGInX()) ;
                convergence_info << relativeL2Error << "," ;

                if(previousError == 0){
                    convergence_info << "," ;
                }
                else{
                    double observedOrder = log(previousError/relativeL2Error)/log(previousCellSize/cellSize) ;
                    convergence_info << observedOrder << "," ;
                }

                previousError = relativeL2Error ;
            }

            convergence_info << elapsed.count() << endl ;

            previousCenterlineG = dom.heatData()->centerlineGInX() ;
            previousCellSize = cellSize ;
        }
    }

    void DOMAngularConvergenceTest(){
                
        // Creates and validates all required input/output paths before starting the simulations.
        PathConfig paths ;
            
        // Stops execution if any required path or output directory setup fails.
        try{
            paths.validate() ;
            paths.ensureOutputDir() ;
        }
        catch(const exception& e){
            cerr << "Path error: " << e.what() << endl ;
            throw runtime_error(e.what()) ;
        }

        string parentFolder = "DOMP1ComparisionRectangularEnclosureAnnularSource" ;
        string caseName = "DOMAngularConvergenceTest" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;

        // Lists the ordinate orders used in the angular convergence study.
        vector<int> ordinateOrders = {2,4,6,8} ;

        // Fixed spatial mesh used for every ordinate order.
        string meshSubCaseName = "320x320x1" ;

        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 100.0 ;

        // Non-scattering medium.
        double sigmaS = 0.0 ;

        // Zero conductivity.
        double k = 0.0 ;

        // Cold black walls in x and y and periodic boundaries in z.
        vector<double> emissivity = {1.0,1.0,1.0,1.0,0.0,0.0} ;

        // Stores the angular convergence data in a directly plottable format.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, caseName) / "angular_convergence.csv") ;
        convergence_info << "ordinate_order,"
                        << "cells_x,"
                        << "cells_y,"
                        << "total_cells_xy,"
                        << "successive_ordinate_relative_L2_change,"
                        << "observed_angular_order,"
                        << "runtime" << endl ;

        double previousError = 0 ;
        double previousAngularSpacing = 0 ;
        map<double,double> previousCenterlineG ;
        bool firstOrdinate = true ;

        for(int ordinate: ordinateOrders){

            string subCaseName = "ordinate" + to_string(ordinate) ;
            string outputCaseName = caseName + "/" + subCaseName ;

            paths.createOutputDirectory(parentFolder, outputCaseName, meshSubCaseName) ;

            // Retrieves the fixed 320x320x1 mesh files.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, meshSubCaseName) ;

            // Defines the prescribed temperature and initial irradiation fields.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the DOM solver and runs the annular-source simulation for the current ordinate order.
            string domCaseName = parentFolder + "/" + outputCaseName ;
            DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, meshSubCaseName, domCaseName, meshSubCaseName, ordinate) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            dom.runDOMRadiationAnnularSourceRectangularEnclosure() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            double angularSpacing = 1.0/ordinate ;

            convergence_info << ordinate << ","
                            << 320 << ","
                            << 320 << ","
                            << 320*320 << "," ;

            if(firstOrdinate){
                convergence_info << ",," ;
                firstOrdinate = false ;
            }
            else{
                double relativeL2Error = calculateCenterlineRelativeL2Error(previousCenterlineG,dom.heatData()->centerlineGInX()) ;
                convergence_info << relativeL2Error << "," ;

                if(previousError == 0){
                    convergence_info << "," ;
                }
                else{
                    double observedOrder = log(previousError/relativeL2Error)/log(previousAngularSpacing/angularSpacing) ;
                    convergence_info << observedOrder << "," ;
                }

                previousError = relativeL2Error ;
            }

            convergence_info << elapsed.count() << endl ;

            previousCenterlineG = dom.heatData()->centerlineGInX() ;
            previousAngularSpacing = angularSpacing ;

            ofstream centerline_G(paths.getOutputDirectory(parentFolder, outputCaseName, meshSubCaseName) / "centerlineG.csv") ;
            centerline_G << "ordinate_order," << ordinate << endl ;
            centerline_G << "mesh," << meshSubCaseName << endl ;
            centerline_G << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_G << "scattering_coefficient," << sigmaS << endl << endl ;
            centerline_G << "x,G_DOM" << endl ;

            for(auto& [coordinate,value]: dom.heatData()->centerlineGInX()){
                centerline_G << coordinate << "," << value << endl ;
            }
        }
    }

    double calculateCenterlineRelativeL2Error(const map<double,double>& coarseCenterline, const map<double,double>& fineCenterline){

        if(coarseCenterline.empty() || fineCenterline.size() < 2){
            throw runtime_error("Insufficient centerline data for relative L2 error calculation.") ;
        }

        double errorSum = 0 ;
        double fineSum = 0 ;

        for(const auto& [coordinate,coarseValue]: coarseCenterline){

            auto upperValue = fineCenterline.lower_bound(coordinate) ;

            double interpolatedFineValue ;

            if(upperValue == fineCenterline.begin()){
                interpolatedFineValue = upperValue->second ;
            }
            else if(upperValue == fineCenterline.end()){
                interpolatedFineValue = prev(upperValue)->second ;
            }
            else{
                auto lowerValue = prev(upperValue) ;
                double interpolationFactor = (coordinate-lowerValue->first)/(upperValue->first-lowerValue->first) ;
                interpolatedFineValue = lowerValue->second+interpolationFactor*(upperValue->second-lowerValue->second) ;
            }

            errorSum += powl(coarseValue-interpolatedFineValue,2) ;
            fineSum += powl(interpolatedFineValue,2) ;
        }

        if(fineSum == 0){
            throw runtime_error("Fine-grid centerline norm is zero.") ;
        }

        return sqrt(errorSum/fineSum) ;
    }

} ;