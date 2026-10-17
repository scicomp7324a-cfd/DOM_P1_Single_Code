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

class P1PreliminaryValidation{

    public:

    void test(){
        p1ParallelPlateSimulationPureRadiationNoParticipatingMediumClosed() ;
        p1ParallelPlateSimulationPureRadiationParticipatingMedium() ;
        p1SquarePlateSimulationRadiationConduction() ;
    }

    // Parallel plate simulation of pure radiation with no participating medium
    void p1ParallelPlateSimulationPureRadiationNoParticipatingMediumClosed(){

    /*
        This function deals with heat transfer between to parallel plates conditioned such that they behave as infinite parallel plates. In this case the west wall temperature is taken to be 1000 K and the east wall temperature is varied and the diagnostic is accumulate for all the cases. The analytical calculation of heat flux is done in 
    */
    
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

    string parentFolder = "P1PreliminaryValidation" ;
    string caseName = "ParallelPlateNoMediumVacuumClosed" ;
    vector<string> subCaseNames = {"300","400","500","600","700","800","900","1000"} ; // East Wall Temperature is taken here as sub case name as it is varying
    //vector<string> subCaseNames = {"300"} ;
    paths.createOutputDirectory(parentFolder, caseName) ;
    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;
    
    // Lists the mesh sizes used in the study.
    string dimensions = "20x5x5" ;

    // Absorption coefficient value. No participating medium
    double a = 0.0 ;

    // Zero scattering
    double sigmaS = 0 ;

    // Zero conductivity
    double k = 0 ;

    // Emissivity values along the walls
    vector<double> emissivity = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0} ; 

    // East Wall temperatures 
    
    double Twest = 1000 ;

    // Diagnostics 
    ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder,caseName) / "heat_diagnostics.csv") ;
    heat_flux_transfer_info << "west_wall_temperature,"
                            << "east_wall_temperature,"
                            << "runtime,"
                            << "east_wall_heat_flux,"
                            << "east_wall_analytical_heat_flux,"
                            << "west_wall_heat_flux,"
                            << "west_wall_analytical_heat_flux,"
                            << "north_wall_heat_flux,"
                            << "south_wall_heat_flux,"
                            << "top_wall_heat_flux,"
                            << "bottom_wall_heat_flux,"
                            << "heat_flux_sum,"
                            << "heat_flux_east_west_wall_sum,"
                            << "west_east_transfer_flux_seen_at_east,"
                            << "west_east_transfer_flux_seen_at_west,"
                            << "east_wall_heat_transfer,"
                            << "west_wall_heat_transfer,"
                            << "north_wall_heat_transfer,"
                            << "south_wall_heat_transfer,"
                            << "top_wall_heat_transfer,"
                            << "bottom_wall_heat_transfer,"
                            << "heat_transfer_sum,"
                            << "heat_transfer_east_west_wall_sum," 
                            << "opposite_wall_heat_transfer_east,"
                            << "opposite_wall_heat_transfer_west" << endl ;

    
    for(string subCaseName: subCaseNames){

        paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

        // Retrieves the mesh files for the current mesh configuration.
        meshFiles = paths.getMeshFilenames(parentFolder, caseName, dimensions) ;

        // Defines the field files for the initial G and T data.
        fieldFiles = {paths.getTFieldFile(parentFolder, caseName, subCaseName), paths.getGFieldFile(parentFolder, caseName, subCaseName)} ;

        // Constructs the SIMPLE solver for the current case and executes the full simulation.
        string p1CaseName = parentFolder + "/" + caseName ;
        P1 p1(meshFiles, fieldFiles, a, sigmaS, k, emissivity, dimensions, p1CaseName, subCaseName) ;
        auto start = std::chrono::high_resolution_clock::now();
        p1.runP1PureRadiationNoParticipatingMedium() ;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        p1.heatData()->calculateWestEastTransferForParallelPlates() ;
        map<string,double> GAnalytical = p1.heatData()->parallelPlateHeatFluxBlackbodyAnalyticalValue(Twest,stod(subCaseName)) ;

        heat_flux_transfer_info << Twest << ","
                                << subCaseName << ","
                                << elapsed.count() << ","
                                << p1.heatData()->avgqw()["east"] << ","
                                << GAnalytical["east"] << ","
                                << p1.heatData()->avgqw()["west"] << ","
                                << GAnalytical["west"] << ","
                                << p1.heatData()->avgqw()["north"] << ","
                                << p1.heatData()->avgqw()["south"] << ","
                                << p1.heatData()->avgqw()["top"] << ","
                                << p1.heatData()->avgqw()["bottom"] << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"] + p1.heatData()->avgqw()["north"] + p1.heatData()->avgqw()["south"] + p1.heatData()->avgqw()["top"] + p1.heatData()->avgqw()["bottom"]) << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatFlux()["east"] << ","
                                << p1.heatData()->oppositeWallHeatFlux()["west"] << ","
                                << p1.heatData()->totalQw()["east"] << ","
                                << p1.heatData()->totalQw()["west"] << ","
                                << p1.heatData()->totalQw()["north"] << ","
                                << p1.heatData()->totalQw()["south"] << ","
                                << p1.heatData()->totalQw()["top"] << ","
                                << p1.heatData()->totalQw()["bottom"] << ","
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"] + p1.heatData()->totalQw()["north"] + p1.heatData()->totalQw()["south"] + p1.heatData()->totalQw()["top"] + p1.heatData()->totalQw()["bottom"]) << "," 
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["east"] << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["west"] << endl ;

        ofstream centerline_G(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineG.csv") ;
        centerline_G << "west_wall_temperature," << Twest << endl ;
        centerline_G << "east_wall_temperature," << subCaseName << endl << endl ;
        centerline_G << "x_centerline_coordinate, G" << endl ;
        for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << val << endl ; }
        centerline_G << endl ;

    }
}

    // Parallel plate simulation of pure radiation with participating medium
    void p1ParallelPlateSimulationPureRadiationParticipatingMedium(){

    /*
        This function deals with heat transfer between to parallel plates conditioned such that they behave as infinite parallel plates. In this case the west wall temperature is taken to be 1000 K and the east wall temperature is kept at 500 K and absorption coefficient is varied. The analytical calculation of heat flux in included for parallel plates. 
    */
    
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

    string parentFolder = "P1PreliminaryValidation" ;
    string caseName = "ParallelPlateParticipatingMediumClosed" ;

    paths.createOutputDirectory(parentFolder, caseName) ;

    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;
    
    // Lists the mesh sizes used in the study.
    string dimensions = "20x5x5" ;

    // Absorption coefficient value. No participating medium
    vector<string> subCaseNames = {"0.1","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","5.0","7.0","9.0","10.0"} ;

    // Zero scattering
    double sigmaS = 0 ;

    // Zero conductivity
    double k = 0 ;

    // Emissivity values along the walls
    vector<double> emissivity = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0} ; 

    // East Wall temperatures 
    
    double Twest = 1000 ;
    double Teast = 500 ;

    // Diagnostics 
    ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder,caseName) / "heat_diagnostics.csv") ;
    heat_flux_transfer_info << "west_wall_temperature,"
                            << "east_wall_temperature,"
                            << "runtime,"
                            << "absorption_coefficient,"
                            << "east_wall_heat_flux,"
                            << "east_wall_analytical_heat_flux,"
                            << "west_wall_heat_flux,"
                            << "west_wall_analytical_heat_flux,"
                            << "north_wall_heat_flux,"
                            << "south_wall_heat_flux,"
                            << "top_wall_heat_flux,"
                            << "bottom_wall_heat_flux,"
                            << "heat_flux_sum,"
                            << "heat_flux_east_west_wall_sum,"
                            << "west_east_transfer_flux_seen_at_east,"
                            << "west_east_transfer_flux_seen_at_west,"
                            << "east_wall_heat_transfer,"
                            << "west_wall_heat_transfer,"
                            << "north_wall_heat_transfer,"
                            << "south_wall_heat_transfer,"
                            << "top_wall_heat_transfer,"
                            << "bottom_wall_heat_transfer,"
                            << "heat_transfer_sum,"
                            << "heat_transfer_east_west_wall_sum," 
                            << "opposite_wall_heat_transfer_east,"
                            << "opposite_wall_heat_transfer_west,"
                            << "dimensionless_heat_flux_east," 
                            << "dimensionless_heat_flux_west," 
                            << "dimensionless_heat_flux_north," 
                            << "dimensionless_heat_flux_south," 
                            << "dimensionless_heat_flux_top," 
                            << "dimensionless_heat_flux_bottom" << endl ;

    // Approximate digitised analytical data from Fig. 4(b) taken from the siegel howell
    ofstream analytical_dimensionless_heat_flux(paths.getOutputDirectory(parentFolder,caseName) / "analytical_dimensionless_heat_flux.csv") ;
    analytical_dimensionless_heat_flux << "absorption_coefficient,dimension_heat_flux_east" << endl ;
    vector<double> tauL_heatFlux = { 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00 };
    vector<double> psi_heatFlux = { 0.81, 0.70, 0.61, 0.55, 0.50, 0.45, 0.42, 0.39, 0.36, 0.34, 0.32, 0.30 };
    for(size_t i = 0; i < tauL_heatFlux.size(); ++i){ analytical_dimensionless_heat_flux << tauL_heatFlux[i] << "," << psi_heatFlux[i] << endl ;}

    ofstream analytical_dimensionless_centerline_temperature(paths.getOutputDirectory(parentFolder,caseName) / "analytical_dimensionless_centerline_temperature.csv") ;
    analytical_dimensionless_centerline_temperature << "x/l,0.1,1.0,2.0,10.0" << endl ;
    vector<double> xOverL = { 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0 };
    vector<double> phi_tau_0p1 = { 0.58, 0.55, 0.54, 0.54, 0.49, 0.46, 0.43 };
    vector<double> phi_tau_1p0 = { 0.80, 0.75, 0.69, 0.58, 0.47, 0.36, 0.24};
    vector<double> phi_tau_2p0 = { 0.88, 0.82, 0.75, 0.60, 0.42, 0.29, 0.15 };
    vector<double> phi_tau_10p0 = { 0.99, 0.91, 0.81, 0.58, 0.33, 0.16, 0.06 };
    for(size_t i = 0; i < xOverL.size(); ++i){analytical_dimensionless_centerline_temperature << xOverL[i] << "," << phi_tau_0p1[i] << "," << phi_tau_1p0[i] << "," << phi_tau_2p0[i] << "," << phi_tau_10p0[i] << endl ;}

    for(string subCaseName: subCaseNames){

        paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

        // Retrieves the mesh files for the current mesh configuration.
        meshFiles = paths.getMeshFilenames(parentFolder, caseName, dimensions) ;

        // Defines the field files for the initial G and T data.
        fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

        // Constructs the SIMPLE solver for the current case and executes the full simulation.
        string p1CaseName = parentFolder + "/" + caseName ;
        P1 p1(meshFiles, fieldFiles, stod(subCaseName), sigmaS, k, emissivity, dimensions, p1CaseName, subCaseName) ;
        auto start = std::chrono::high_resolution_clock::now();
        p1.runP1PureRadiationParticipatingMedium() ;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        p1.heatData()->calculateWestEastTransferForParallelPlates() ;
        map<string,double> GAnalytical = p1.heatData()->parallelPlateHeatFluxBlackbodyAnalyticalValue(Twest,Teast) ;

        heat_flux_transfer_info << Twest << ","
                                << Teast << ","
                                << elapsed.count() << ","
                                << subCaseName << ","
                                << p1.heatData()->avgqw()["east"] << ","
                                << GAnalytical["east"] << ","
                                << p1.heatData()->avgqw()["west"] << ","
                                << GAnalytical["west"] << ","
                                << p1.heatData()->avgqw()["north"] << ","
                                << p1.heatData()->avgqw()["south"] << ","
                                << p1.heatData()->avgqw()["top"] << ","
                                << p1.heatData()->avgqw()["bottom"] << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"] + p1.heatData()->avgqw()["north"] + p1.heatData()->avgqw()["south"] + p1.heatData()->avgqw()["top"] + p1.heatData()->avgqw()["bottom"]) << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatFlux()["east"] << ","
                                << p1.heatData()->oppositeWallHeatFlux()["west"] << ","
                                << p1.heatData()->totalQw()["east"] << ","
                                << p1.heatData()->totalQw()["west"] << ","
                                << p1.heatData()->totalQw()["north"] << ","
                                << p1.heatData()->totalQw()["south"] << ","
                                << p1.heatData()->totalQw()["top"] << ","
                                << p1.heatData()->totalQw()["bottom"] << ","
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"] + p1.heatData()->totalQw()["north"] + p1.heatData()->totalQw()["south"] + p1.heatData()->totalQw()["top"] + p1.heatData()->totalQw()["bottom"]) << "," 
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["east"] << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["west"] << ","
                                << p1.heatData()->dimensionlessQw()["east"] << ","
                                << p1.heatData()->dimensionlessQw()["west"] << ","
                                << p1.heatData()->dimensionlessQw()["north"] << ","
                                << p1.heatData()->dimensionlessQw()["south"] << ","
                                << p1.heatData()->dimensionlessQw()["top"] << ","
                                << p1.heatData()->dimensionlessQw()["bottom"] << endl ;

        ofstream centerline_G(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineG.csv") ;
        ofstream centerline_T(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineT.csv") ;
        ofstream dimensionless_centerline_T(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "dimensionlessCenterlineT.csv") ;

        centerline_G << "west_wall_temperature," << Twest << endl ;
        centerline_G << "east_wall_temperature," << Teast << endl ;
        centerline_G << "absorption_coefficient,"<< subCaseName << endl << endl ;
        centerline_G << "x_over_L, G" << endl ;

        centerline_T << "west_wall_temperature," << Twest << endl ;
        centerline_T << "east_wall_temperature," << Teast << endl ;
        centerline_T << "absorption_coefficient,"<< subCaseName << endl << endl ;
        centerline_T << "x_over_L, T" << endl ;

        dimensionless_centerline_T << "west_wall_temperature," << Twest << endl ;
        dimensionless_centerline_T << "east_wall_temperature," << Teast << endl ;
        dimensionless_centerline_T << "absorption_coefficient,"<< subCaseName << endl << endl ;
        dimensionless_centerline_T << "x_over_L, dimensionless_T" << endl ;

        for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << val << endl ; }
        for(auto& [key,val]: p1.heatData()->centerlineTemperature()){ centerline_T << key << "," << val << endl ;}
        for(auto& [key,val]: p1.heatData()->dimensionlessCenterlineTemperature()){ dimensionless_centerline_T << key << "," << val << endl ;}
        
        centerline_G << endl ;
        centerline_T << endl ;
        dimensionless_centerline_T << endl ;

    }
    
}

    // Square plate simulation of pure radiation with participating medium
    void p1SquarePlateSimulationRadiationConduction(){

    /*
        This function deals with heat transfer between to parallel plates conditioned such that they behave as infinite parallel plates. In this case the west wall temperature is taken to be 1000 K and the east wall temperature is kept at 500 K and absorption coefficient is varied. The analytical calculation of heat flux in included for parallel plates. 
    */
    
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

    string parentFolder = "P1PreliminaryValidation" ;
    string caseName = "SquarePlateRadiationConductionClosed" ;

    paths.createOutputDirectory(parentFolder, caseName) ;
    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;
    
    // Lists the mesh sizes used in the study.
    string dimensions = "40x5x40" ;

    // Absorption coefficient value. No participating medium
    double a = 1 ;

    // Zero scattering
    double sigmaS = 0 ;
    double sigma = 5.670374419e-8 ;

    // Emissivity values along the walls
    vector<double> emissivity = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0} ; 

    // East Wall temperatures 
    
    double Tbottom = 1000 ;
    double Ttop = Tbottom/2 ;
    double Twest = Tbottom/2 ;
    double Teast = Tbottom/2 ;
    
    // Stark number 
    vector<double> N = {1, 0.1, 0.05} ;
    // Conductivity is the subcase in this case 
    double sc = 4*sigma*(powl(Tbottom,3) - powl(Ttop,3))/a ;
    vector<string> subCaseNames = {to_string(sc*N[0]), to_string(sc*N[1]), to_string(sc*N[2])} ;

    // Diagnostics 
    ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder,caseName) / "heat_diagnostics.csv") ;
    heat_flux_transfer_info << "bottom_wall_temperature,"
                            << "top_wall_temperature,"
                            << "west_wall_temperature,"
                            << "east_wall_temperature,"
                            << "runtime,"
                            << "absorption_coefficient,"
                            << "conductivity,"
                            << "east_wall_heat_flux,"
                            << "west_wall_heat_flux,"
                            << "north_wall_heat_flux,"
                            << "south_wall_heat_flux,"
                            << "top_wall_heat_flux,"
                            << "bottom_wall_heat_flux,"
                            << "heat_flux_sum,"
                            << "heat_flux_east_west_wall_sum,"
                            << "west_east_transfer_flux_seen_at_east,"
                            << "west_east_transfer_flux_seen_at_west,"
                            << "east_wall_heat_transfer,"
                            << "west_wall_heat_transfer,"
                            << "north_wall_heat_transfer,"
                            << "south_wall_heat_transfer,"
                            << "top_wall_heat_transfer,"
                            << "bottom_wall_heat_transfer,"
                            << "heat_transfer_sum,"
                            << "heat_transfer_east_west_wall_sum," 
                            << "opposite_wall_heat_transfer_east,"
                            << "opposite_wall_heat_transfer_west,"
                            << "dimensionless_heat_flux_east," 
                            << "dimensionless_heat_flux_west," 
                            << "dimensionless_heat_flux_north," 
                            << "dimensionless_heat_flux_south," 
                            << "dimensionless_heat_flux_top," 
                            << "dimensionless_heat_flux_bottom" << endl ;

    ofstream analytical_dimensionless_centerline_temperature(paths.getOutputDirectory(parentFolder, caseName) / "analytical_dimensionless_centerline_temperature.csv") ;
    analytical_dimensionless_centerline_temperature << "y/l,N0.05,N0.1,N1" << endl ;
    // Approximate digitised Rousse 2000 validation data from Fig. 7

    vector<double> yOverLy_fig7 = { 0.00, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00 };
    vector<double> TbyTw_N_1_fig7 = { 1.000, 0.951, 0.901, 0.815, 0.737, 0.677, 0.629, 0.592, 0.562, 0.539, 0.519, 0.500 };
    vector<double> TbyTw_N_0p1_fig7 = { 1.000, 0.952, 0.905, 0.822, 0.757, 0.709, 0.664, 0.624, 0.593, 0.561, 0.534, 0.500};
    vector<double> TbyTw_N_0p05_fig7 = { 1.000, 0.953, 0.907, 0.828, 0.768, 0.726, 0.685, 0.649, 0.613, 0.578, 0.545, 0.500 };

    for(size_t i = 0; i < yOverLy_fig7.size(); ++i){analytical_dimensionless_centerline_temperature << yOverLy_fig7[i] << "," <<  TbyTw_N_0p05_fig7[i] << "," << TbyTw_N_0p1_fig7[i] << "," << TbyTw_N_1_fig7[i] << endl ;}

    
    for(string subCaseName: subCaseNames){

        paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

        // Retrieves the mesh files for the current mesh configuration.
        meshFiles = paths.getMeshFilenames(parentFolder, caseName, dimensions) ;

        // Defines the field files for the initial G and T data.
        fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

        // Constructs the SIMPLE solver for the current case and executes the full simulation.
        string p1CaseName = parentFolder + "/" + caseName ;
        P1 p1(meshFiles, fieldFiles, a, sigmaS, stod(subCaseName), emissivity, dimensions, p1CaseName, subCaseName) ;
        auto start = std::chrono::high_resolution_clock::now();
        p1.runP1RadiationConduction() ;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        p1.heatData()->calculateWestEastTransferForParallelPlates() ;

        heat_flux_transfer_info << Tbottom << ","
                                << Ttop << ","
                                << Twest << ","
                                << Teast << ","
                                << elapsed.count() << ","
                                << a << ","
                                << subCaseName << ","
                                << p1.heatData()->avgqw()["east"] << ","
                                << p1.heatData()->avgqw()["west"] << ","
                                << p1.heatData()->avgqw()["north"] << ","
                                << p1.heatData()->avgqw()["south"] << ","
                                << p1.heatData()->avgqw()["top"] << ","
                                << p1.heatData()->avgqw()["bottom"] << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"] + p1.heatData()->avgqw()["north"] + p1.heatData()->avgqw()["south"] + p1.heatData()->avgqw()["top"] + p1.heatData()->avgqw()["bottom"]) << ","
                                << (p1.heatData()->avgqw()["east"] + p1.heatData()->avgqw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatFlux()["east"] << ","
                                << p1.heatData()->oppositeWallHeatFlux()["west"] << ","
                                << p1.heatData()->totalQw()["east"] << ","
                                << p1.heatData()->totalQw()["west"] << ","
                                << p1.heatData()->totalQw()["north"] << ","
                                << p1.heatData()->totalQw()["south"] << ","
                                << p1.heatData()->totalQw()["top"] << ","
                                << p1.heatData()->totalQw()["bottom"] << ","
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"] + p1.heatData()->totalQw()["north"] + p1.heatData()->totalQw()["south"] + p1.heatData()->totalQw()["top"] + p1.heatData()->totalQw()["bottom"]) << "," 
                                << (p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"]) << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["east"] << ","
                                << p1.heatData()->oppositeWallHeatTransfer()["west"] << ","
                                << p1.heatData()->dimensionlessQw()["east"] << ","
                                << p1.heatData()->dimensionlessQw()["west"] << ","
                                << p1.heatData()->dimensionlessQw()["north"] << ","
                                << p1.heatData()->dimensionlessQw()["south"] << ","
                                << p1.heatData()->dimensionlessQw()["top"] << ","
                                << p1.heatData()->dimensionlessQw()["bottom"] << endl ;

        ofstream centerline_G(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineGInX.csv") ;
        ofstream centerline_T(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineTInX.csv") ;
        ofstream dimensionless_centerline_T(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "dimensionlessCenterlineTInX.csv") ;

        ofstream centerline_G_In_Z(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineGInZ.csv") ;
        ofstream centerline_T_In_Z(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineTInZ.csv") ;
        ofstream dimensionless_centerline_T_In_Z(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "dimensionlessCenterlineTInZ.csv") ;

        centerline_G << "west_wall_temperature," << Twest << endl ;
        centerline_G << "east_wall_temperature," << Teast << endl ;
        centerline_G << "top_wall_temperature," << Ttop << endl ;
        centerline_G << "bottom_wall_temperature," << Tbottom << endl ;
        centerline_G << "conductivity,"<< subCaseName << endl << endl ;
        centerline_G << "x_over_L, G" << endl ;

        centerline_T << "west_wall_temperature," << Twest << endl ;
        centerline_T << "east_wall_temperature," << Teast << endl ;
        centerline_T << "top_wall_temperature," << Ttop << endl ;
        centerline_T << "bottom_wall_temperature," << Tbottom << endl ;
        centerline_T << "conductivity,"<< subCaseName << endl << endl ;
        centerline_T << "x_over_L, T" << endl ;

        dimensionless_centerline_T << "west_wall_temperature," << Twest << endl ;
        dimensionless_centerline_T << "east_wall_temperature," << Teast << endl ;
        dimensionless_centerline_T << "top_wall_temperature," << Ttop << endl ;
        dimensionless_centerline_T << "bottom_wall_temperature," << Tbottom << endl ;
        dimensionless_centerline_T << "conductivity,"<< subCaseName << endl << endl ;
        dimensionless_centerline_T << "x_over_L, dimensionless_T" << endl ;

        centerline_G_In_Z << "west_wall_temperature," << Twest << endl ;
        centerline_G_In_Z << "east_wall_temperature," << Teast << endl ;
        centerline_G_In_Z << "top_wall_temperature," << Ttop << endl ;
        centerline_G_In_Z << "bottom_wall_temperature," << Tbottom << endl ;
        centerline_G_In_Z << "conductivity,"<< subCaseName << endl << endl ;
        centerline_G_In_Z << "z_over_L, G" << endl ;

        centerline_T_In_Z << "west_wall_temperature," << Twest << endl ;
        centerline_T_In_Z << "east_wall_temperature," << Teast << endl ;
        centerline_T_In_Z << "top_wall_temperature," << Ttop << endl ;
        centerline_T_In_Z << "bottom_wall_temperature," << Tbottom << endl ;
        centerline_T_In_Z << "conductivity,"<< subCaseName << endl << endl ;
        centerline_T_In_Z << "z_over_L, T" << endl ;

        dimensionless_centerline_T_In_Z << "west_wall_temperature," << Twest << endl ;
        dimensionless_centerline_T_In_Z << "east_wall_temperature," << Teast << endl ;
        dimensionless_centerline_T_In_Z << "top_wall_temperature," << Ttop << endl ;
        dimensionless_centerline_T_In_Z << "bottom_wall_temperature," << Tbottom << endl ;
        dimensionless_centerline_T_In_Z << "conductivity,"<< subCaseName << endl << endl ;
        dimensionless_centerline_T_In_Z << "z_over_L, dimensionless_T" << endl ;

        for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << val << endl ; }
        for(auto& [key,val]: p1.heatData()->centerlineTemperature()){ centerline_T << key << "," << val << endl ;}
        for(auto& [key,val]: p1.heatData()->dimensionlessCenterlineTemperature()){ dimensionless_centerline_T << key << "," << val << endl ;}

        for(auto& [key,val]: p1.heatData()->centerlineGInZ()){ centerline_G_In_Z << key << "," << val << endl ; }
        for(auto& [key,val]: p1.heatData()->centerlineTemperatureInZ()){ centerline_T_In_Z << key << "," << val << endl ;}
        for(auto& [key,val]: p1.heatData()->dimensionlessCenterlineTemperatureInZ()){ dimensionless_centerline_T_In_Z << key << "," << val << endl ;}
        
        centerline_G << endl ;
        centerline_T << endl ;
        dimensionless_centerline_T << endl ;

    }
}

} ;