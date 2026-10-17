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

class DOMP1ComparisonAnisotropicAndIsotropicScattering{

    public:

    void test(){ 
        P1Scattering() ;
        DOMScattering() ;
        DOMAngularConvergence() ;
    }

    void DOMAngularConvergence(){

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

        string parentFolder = "DOMP1ComparisonAnisotropicAndIsotropicScattering" ;
        string caseName = "DOMAngularConvergence" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;

        // Lists the ordinate orders used in the angular convergence study.
        vector<int> ordinateOrders = {2,4,6,8} ;

        // Fixed spatial mesh used for every ordinate order.
        string meshSubCaseName = "300x1x1" ;

        // Absorption coefficient is zero because the medium is purely scattering.
        double absorptionCoefficient = 0.0 ;

        // Since absorption is zero, the scattering coefficient is equal to the extinction coefficient.
        vector<double> sigmaSValues = {0.1,1.0,10.0} ;

        // Isotropic and forward-scattering cases used in the angular convergence study.
        vector<double> anisotropyFactors = {0.0,1.0} ;

        // Zero conductivity.
        double k = 0.0 ;

        // Black west and east walls with periodic transverse boundaries.
        vector<double> emissivity = {1.0,1.0,0.0,0.0,0.0,0.0} ;

        // Stores the angular convergence data in a directly plottable format.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, caseName) / "angular_convergence.csv") ;
        convergence_info << "anisotropy_factor,"
                        << "extinction_coefficient,"
                        << "scattering_coefficient,"
                        << "ordinate_order,"
                        << "mesh,"
                        << "east_wall_heat_flux,"
                        << "dimensionless_east_wall_heat_flux,"
                        << "runtime" << endl ;

        for(double anisotropyFactor: anisotropyFactors){
            for(double sigmaS: sigmaSValues){
                for(int ordinate: ordinateOrders){

                    string anisotropySubCaseName = "anisotropyFactor" + to_string(anisotropyFactor) ;
                    string extinctionCoefficientSubCaseName = "extinctionCoefficient" + to_string(sigmaS) + "_sigmaS" + to_string(sigmaS) ;
                    string ordinateSubCaseName = "ordinate" + to_string(ordinate) ;

                    string outputCaseName = caseName + "/"
                                        + anisotropySubCaseName + "/"
                                        + extinctionCoefficientSubCaseName ;

                    paths.createOutputDirectory(parentFolder, outputCaseName, ordinateSubCaseName) ;

                    // Retrieves the fixed 100x1x1 mesh files.
                    meshFiles = paths.getMeshFilenames(parentFolder, caseName, meshSubCaseName) ;

                    // Defines the prescribed temperature and initial irradiation fields.
                    fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

                    // Constructs the DOM solver and runs the scattering simulation for the current ordinate order.
                    string domCaseName = parentFolder + "/" + outputCaseName ;
                    DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, meshSubCaseName, domCaseName, ordinateSubCaseName, ordinate, anisotropyFactor) ;

                    auto start = std::chrono::high_resolution_clock::now() ;
                    dom.runDOMRadiationScattering() ;
                    auto end = std::chrono::high_resolution_clock::now() ;
                    std::chrono::duration<double> elapsed = end-start ;

                    // Calculates the dimensionless wall heat flux using the wall temperature as the reference temperature.
                    dom.heatData()->calculateDimensionlessHeatFlux2() ;

                    double eastWallHeatFlux = dom.heatData()->averageWallNormalHeatFlux()["east"] ;
                    double dimensionlessEastWallHeatFlux = dom.heatData()->dimensionlessQw2()["east"] ;

                    convergence_info << anisotropyFactor << ","
                                    << sigmaS << ","
                                    << sigmaS << ","
                                    << ordinate << ","
                                    << meshSubCaseName << ","
                                    << eastWallHeatFlux << ","
                                    << dimensionlessEastWallHeatFlux << ","
                                    << elapsed.count() << endl ;

                    ofstream angular_result(paths.getOutputDirectory(parentFolder, outputCaseName, ordinateSubCaseName) / "angular_convergence_result.csv") ;
                    angular_result << "anisotropy_factor,"
                                << "extinction_coefficient,"
                                << "scattering_coefficient,"
                                << "ordinate_order,"
                                << "mesh,"
                                << "east_wall_heat_flux,"
                                << "dimensionless_east_wall_heat_flux,"
                                << "runtime" << endl ;

                    angular_result << anisotropyFactor << ","
                                << sigmaS << ","
                                << sigmaS << ","
                                << ordinate << ","
                                << meshSubCaseName << ","
                                << eastWallHeatFlux << ","
                                << dimensionlessEastWallHeatFlux << ","
                                << elapsed.count() << endl ;
                }
            }
        }
    }


    void P1Scattering(){
        
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

        string parentFolder = "DOMP1ComparisonAnisotropicAndIsotropicScattering" ;
        string caseName = "P1Scattering" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;

        // Fixed spatial mesh used for every extinction coefficient and anisotropy factor.
        string meshSubCaseName = "300x1x1" ;

        // Physical plate separation.
        double plateSeparation = 1.0 ;

        // Absorption coefficient is zero because the medium is purely scattering.
        double absorptionCoefficient = 0.0 ;

        // Since absorption is zero, the scattering coefficient is equal to the extinction coefficient.
        vector<double> sigmaSValues = {0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10.0} ;

        // Anisotropy factors used for the dimensionless wall heat-flux comparison.
        vector<double> anisotropyFactors = {1.0,0.5,0.0,-0.5,-1.0} ;

        // Zero conductivity.
        double k = 0.0 ;

        // Black west and east walls with periodic transverse boundaries.
        vector<double> emissivity = {1.0,1.0,0.0,0.0,0.0,0.0} ;

        // Stores the dimensionless wall heat-flux results required for Test 2.
        ofstream heat_flux_info(paths.getOutputDirectory(parentFolder, caseName) / "dimensionless_heat_flux.csv") ;
        heat_flux_info << "anisotropy_factor,"
                    << "asymmetry_factor,"
                    << "extinction_coefficient,"
                    << "scattering_coefficient,"
                    << "transport_scattering_coefficient,"
                    << "mesh,"
                    << "east_wall_heat_flux,"
                    << "dimensionless_east_wall_heat_flux,"
                    << "runtime" << endl ;

        for(double anisotropyFactor: anisotropyFactors){
            for(double sigmaS: sigmaSValues){

                // Converts the anisotropy factor used by Singh and Mishra to the P1 asymmetry factor.
                double asymmetryFactor = anisotropyFactor/3.0 ;

                // P1 accounts for anisotropic scattering through the transport scattering coefficient.
                double transportScatteringCoefficient = sigmaS*(1.0-asymmetryFactor) ;

                string anisotropySubCaseName = "anisotropyFactor" + to_string(anisotropyFactor) ;
                string extinctionCoefficientSubCaseName = "extinctionCoefficient" + to_string(sigmaS) + "_sigmaS" + to_string(sigmaS) ;

                string outputCaseName = caseName + "/"
                                    + anisotropySubCaseName + "/"
                                    + extinctionCoefficientSubCaseName ;

                paths.createOutputDirectory(parentFolder, outputCaseName, meshSubCaseName) ;

                // Retrieves the fixed 100x1x1 mesh files.
                meshFiles = paths.getMeshFilenames(parentFolder, caseName, meshSubCaseName) ;

                // Defines the prescribed temperature and initial irradiation fields.
                fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

                // Constructs the P1 solver and runs the scattering simulation.
                string p1CaseName = parentFolder + "/" + outputCaseName ;
                P1 p1(meshFiles, fieldFiles, absorptionCoefficient, sigmaS, k, emissivity, meshSubCaseName, p1CaseName, meshSubCaseName, anisotropyFactor) ;

                auto start = std::chrono::high_resolution_clock::now() ;
                p1.runP1RadiationScattering() ;
                auto end = std::chrono::high_resolution_clock::now() ;
                std::chrono::duration<double> elapsed = end-start ;

                // Calculates the wall and centreline diagnostics required for Tests 2 and 3.
                //heatData_->calculateWallHeatFlux();
                //heatData_->calculateAvgTempartures();
                //heatData_->calculateDimensionlessQw();
                //heatData_->calculateDimensionlessSlipTemperature();
                //heatData_->calculateCenterlineTemperatureAndG();
                //heatData_->calculateDimensionlessCenterlineTemperature();

                double eastWallHeatFlux = p1.heatData()->avgqw()["east"] ;
                double dimensionlessEastWallHeatFlux = p1.heatData()->dimensionlessQw()["east"] ;

                heat_flux_info << anisotropyFactor << ","
                            << asymmetryFactor << ","
                            << sigmaS << ","
                            << sigmaS << ","
                            << transportScatteringCoefficient << ","
                            << meshSubCaseName << ","
                            << eastWallHeatFlux << ","
                            << dimensionlessEastWallHeatFlux << ","
                            << elapsed.count() << endl ;

                ofstream heat_flux_result(paths.getOutputDirectory(parentFolder, outputCaseName, meshSubCaseName) / "dimensionless_heat_flux.csv") ;
                heat_flux_result << "anisotropy_factor,"
                                << "asymmetry_factor,"
                                << "extinction_coefficient,"
                                << "scattering_coefficient,"
                                << "transport_scattering_coefficient,"
                                << "mesh,"
                                << "east_wall_heat_flux,"
                                << "dimensionless_east_wall_heat_flux,"
                                << "runtime" << endl ;

                heat_flux_result << anisotropyFactor << ","
                                << asymmetryFactor << ","
                                << sigmaS << ","
                                << sigmaS << ","
                                << transportScatteringCoefficient << ","
                                << meshSubCaseName << ","
                                << eastWallHeatFlux << ","
                                << dimensionlessEastWallHeatFlux << ","
                                << elapsed.count() << endl ;

                // Test 3 requires emissive-power profiles only for a = 0, +1 and -1 and beta = 0.1, 1, 3, 5 and 10.
                bool emissivePowerAnisotropy = anisotropyFactor == 0.0 || anisotropyFactor == 1.0 || anisotropyFactor == -1.0 ;
                bool emissivePowerExtinctionCoefficient = sigmaS == 0.1 || sigmaS == 1.0 || sigmaS == 3.0 || sigmaS == 5.0 || sigmaS == 10.0 ;

                if(emissivePowerAnisotropy && emissivePowerExtinctionCoefficient){

                    auto firstCenterlineCoordinate = p1.heatData()->dimensionlessCenterlineTemperature2().begin() ;
                    auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
                    double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
                    double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

                    ofstream emissive_power(paths.getOutputDirectory(parentFolder, outputCaseName, meshSubCaseName) / "dimensionless_emissive_power.csv") ;
                    emissive_power << "anisotropy_factor," << anisotropyFactor << endl ;
                    emissive_power << "asymmetry_factor," << asymmetryFactor << endl ;
                    emissive_power << "extinction_coefficient," << sigmaS << endl ;
                    emissive_power << "scattering_coefficient," << sigmaS << endl ;
                    emissive_power << "transport_scattering_coefficient," << transportScatteringCoefficient << endl ;
                    emissive_power << "plate_separation," << plateSeparation << endl ;
                    emissive_power << "mesh," << meshSubCaseName << endl << endl ;
                    emissive_power << "x,x_over_X,dimensionless_emissive_power_P1" << endl ;

                    for(auto& [coordinate,value]: p1.heatData()->dimensionlessCenterlineTemperature2()){
                        emissive_power << coordinate << ","
                                    << (coordinate-westWallLocation)/plateSeparation << ","
                                    << value << endl ;
                    }
                }
            }
        }
    }


    void DOMScattering(){
    
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

    string parentFolder = "DOMP1ComparisonAnisotropicAndIsotropicScattering" ;
    string caseName = "DOMScattering" ;

    paths.createOutputDirectory(parentFolder, caseName) ;

    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;

    // Lists the ordinate orders used for the DOM scattering study.
    vector<int> ordinateOrders = {2,4,6,8} ;

    // Fixed spatial mesh used for every ordinate order, extinction coefficient and anisotropy factor.
    string meshSubCaseName = "300x1x1" ;

    // Physical plate separation.
    double plateSeparation = 1.0 ;

    // Absorption coefficient is zero because the medium is purely scattering.
    double absorptionCoefficient = 0.0 ;

    // Since absorption is zero, the scattering coefficient is equal to the extinction coefficient.
    vector<double> sigmaSValues = {0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,6.0,7.0,8.0,9.0,10.0} ;

    // Anisotropy factors used for the dimensionless wall heat-flux comparison.
    vector<double> anisotropyFactors = {1.0,0.5,0.0,-0.5,-1.0} ;

    // Zero conductivity.
    double k = 0.0 ;

    // Black west and east walls with periodic transverse boundaries.
    vector<double> emissivity = {1.0,1.0,0.0,0.0,0.0,0.0} ;

    // Stores the dimensionless wall heat-flux results required for Test 2.
    ofstream heat_flux_info(paths.getOutputDirectory(parentFolder, caseName) / "dimensionless_heat_flux.csv") ;
    heat_flux_info << "anisotropy_factor,"
                   << "extinction_coefficient,"
                   << "scattering_coefficient,"
                   << "ordinate_order,"
                   << "mesh,"
                   << "east_wall_heat_flux,"
                   << "dimensionless_east_wall_heat_flux,"
                   << "runtime" << endl ;

    for(double anisotropyFactor: anisotropyFactors){
        for(double sigmaS: sigmaSValues){
            for(int ordinate: ordinateOrders){

                string anisotropySubCaseName = "anisotropyFactor" + to_string(anisotropyFactor) ;
                string extinctionCoefficientSubCaseName = "extinctionCoefficient" + to_string(sigmaS) + "_sigmaS" + to_string(sigmaS) ;
                string ordinateSubCaseName = "ordinate" + to_string(ordinate) ;

                string outputCaseName = caseName + "/"
                                      + anisotropySubCaseName + "/"
                                      + extinctionCoefficientSubCaseName ;

                paths.createOutputDirectory(parentFolder, outputCaseName, ordinateSubCaseName) ;

                // Retrieves the fixed 100x1x1 mesh files.
                meshFiles = paths.getMeshFilenames(parentFolder, caseName, meshSubCaseName) ;

                // Defines the prescribed temperature and initial irradiation fields.
                fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

                // Constructs the DOM solver and runs the scattering simulation for the current ordinate order.
                string domCaseName = parentFolder + "/" + outputCaseName ;
                DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, meshSubCaseName, domCaseName, ordinateSubCaseName, ordinate, anisotropyFactor) ;

                auto start = std::chrono::high_resolution_clock::now() ;
                dom.runDOMRadiationScattering() ;
                auto end = std::chrono::high_resolution_clock::now() ;
                std::chrono::duration<double> elapsed = end-start ;

                // Calculates the wall and centreline diagnostics required for Tests 2 and 3.
                dom.heatData()->calculateHeatFlux() ;
                dom.heatData()->calculateIrradiation() ;
                dom.heatData()->calculateAverages() ;
                dom.heatData()->calculateAvgTempartures() ;
                dom.heatData()->calculateCenterlineTemperatureAndGInX() ;
                dom.heatData()->calculateDimensionlessHeatFlux2() ;
                dom.heatData()->calculateDimensionlessCenterlineTemperatureInX2() ;

                double eastWallHeatFlux = dom.heatData()->averageWallNormalHeatFlux()["east"] ;
                double dimensionlessEastWallHeatFlux = dom.heatData()->dimensionlessQw2()["east"] ;

                heat_flux_info << anisotropyFactor << ","
                               << sigmaS << ","
                               << sigmaS << ","
                               << ordinate << ","
                               << meshSubCaseName << ","
                               << eastWallHeatFlux << ","
                               << dimensionlessEastWallHeatFlux << ","
                               << elapsed.count() << endl ;

                ofstream heat_flux_result(paths.getOutputDirectory(parentFolder, outputCaseName, ordinateSubCaseName) / "dimensionless_heat_flux.csv") ;
                heat_flux_result << "anisotropy_factor,"
                                 << "extinction_coefficient,"
                                 << "scattering_coefficient,"
                                 << "ordinate_order,"
                                 << "mesh,"
                                 << "east_wall_heat_flux,"
                                 << "dimensionless_east_wall_heat_flux,"
                                 << "runtime" << endl ;

                heat_flux_result << anisotropyFactor << ","
                                 << sigmaS << ","
                                 << sigmaS << ","
                                 << ordinate << ","
                                 << meshSubCaseName << ","
                                 << eastWallHeatFlux << ","
                                 << dimensionlessEastWallHeatFlux << ","
                                 << elapsed.count() << endl ;

                // Test 3 requires emissive-power profiles only for a = 0, +1 and -1 and beta = 0.1, 1, 3, 5 and 10.
                bool emissivePowerAnisotropy = anisotropyFactor == 0.0 || anisotropyFactor == 1.0 || anisotropyFactor == -1.0 ;
                bool emissivePowerExtinctionCoefficient = sigmaS == 0.1 || sigmaS == 1.0 || sigmaS == 3.0 || sigmaS == 5.0 || sigmaS == 10.0 ;

                if(emissivePowerAnisotropy && emissivePowerExtinctionCoefficient){

                    auto firstCenterlineCoordinate = dom.heatData()->dimensionlessCenterlineTemperatureInX2().begin() ;
                    auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
                    double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
                    double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

                    ofstream emissive_power(paths.getOutputDirectory(parentFolder, outputCaseName, ordinateSubCaseName) / "dimensionless_emissive_power.csv") ;
                    emissive_power << "anisotropy_factor," << anisotropyFactor << endl ;
                    emissive_power << "extinction_coefficient," << sigmaS << endl ;
                    emissive_power << "scattering_coefficient," << sigmaS << endl ;
                    emissive_power << "ordinate_order," << ordinate << endl ;
                    emissive_power << "plate_separation," << plateSeparation << endl ;
                    emissive_power << "mesh," << meshSubCaseName << endl << endl ;
                    emissive_power << "x,x_over_X,dimensionless_emissive_power_DOM" << endl ;

                    for(auto& [coordinate,value]: dom.heatData()->dimensionlessCenterlineTemperatureInX2()){
                        emissive_power << coordinate << ","
                                       << (coordinate-westWallLocation)/plateSeparation << ","
                                       << value << endl ;
                    }
                }
            }
        }
    }
}

} ;