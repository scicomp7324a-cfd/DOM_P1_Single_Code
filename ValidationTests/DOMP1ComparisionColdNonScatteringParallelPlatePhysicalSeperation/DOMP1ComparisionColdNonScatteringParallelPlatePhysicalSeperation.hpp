# pragma once

# include "PathConfig.hpp"

# include <iostream>
# include <vector>
# include <array>
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

struct ParallelPlateAnalyticalData{
    map<double,double> centerlineG ;
    map<double,double> centerlineFlux ;
    map<string,double> wallHeatFlux ;
};

class DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation{

    public:

    DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation(){}

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

        string parentFolder = "DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation" ;
        string caseName = "P1RadiativeEnergyAndFluxProfileTest" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
            
        // Lists the physical plate separation cases used in the study.
        vector<string> subCaseNames = {"300x1x1_1", "300x1x1_10"} ;
        
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 2.0 ;

        // Zero scattering.
        double sigmaS = 0 ;

        // Zero conductivity.
        double k = 0 ;

        // Emissivity values along the walls.
        vector<double> emissivity = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0} ; 

        // Both plates are maintained at 500 K.
        double Twest = 500 ;
        double Teast = 500 ;

        // Diagnostics.
        ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder, caseName) / "heat_diagnostics.csv") ;
        heat_flux_transfer_info << "west_wall_temperature,"
                                << "east_wall_temperature,"
                                << "runtime,"
                                << "absorption_coefficient,"
                                << "plate_separation,"
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

        for(string subCaseName: subCaseNames){

            double plateSeparation = 0 ;

            if(subCaseName == "300x1x1_1"){
                plateSeparation = 1.0 ;
            }
            else if(subCaseName == "300x1x1_10"){
                plateSeparation = 10.0 ;
            }
            else{
                throw runtime_error("Unknown physical plate separation subcase.") ;
            }

            paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the field files for the initial G and T data.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the P1 solver for the current case and executes the full simulation.
            string p1CaseName = parentFolder + "/" + caseName ;
            P1 p1(meshFiles, fieldFiles, absorptionCoefficient, sigmaS, k, emissivity, subCaseName, p1CaseName, subCaseName) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            p1.runP1PureRadiationParticipatingMediumTConstant() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            p1.heatData()->calculateWestEastTransferForParallelPlates() ;
            p1.heatData()->calculateCenterlineTemperatureAndG() ;
            p1.heatData()->calculateCenterlineFlux() ;

            ParallelPlateAnalyticalData analyticalData = calculateParallelPlateColdMediumAnalyticalValues(p1.heatData()->centerlineG(),plateSeparation,Twest,Teast,absorptionCoefficient) ;
            auto firstCenterlineCoordinate = p1.heatData()->centerlineG().begin() ;
            auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
            double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
            double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

            heat_flux_transfer_info << Twest << ","
                                    << Teast << ","
                                    << elapsed.count() << ","
                                    << absorptionCoefficient << ","
                                    << plateSeparation << ","
                                    << p1.heatData()->avgqw()["east"] << ","
                                    << analyticalData.wallHeatFlux["east"] << ","
                                    << p1.heatData()->avgqw()["west"] << ","
                                    << analyticalData.wallHeatFlux["west"] << ","
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
            ofstream centerline_flux(paths.getOutputDirectory(parentFolder, caseName, subCaseName) / "centerlineFlux.csv") ;

            centerline_G << "west_wall_temperature," << Twest << endl ;
            centerline_G << "east_wall_temperature," << Teast << endl ;
            centerline_G << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_G << "plate_separation," << plateSeparation << endl ;
            centerline_G << "sub_case," << subCaseName << endl << endl ;
            centerline_G << "x,x_over_L,G_P1,G_analytical" << endl ;

            centerline_flux << "west_wall_temperature," << Twest << endl ;
            centerline_flux << "east_wall_temperature," << Teast << endl ;
            centerline_flux << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_flux << "plate_separation," << plateSeparation << endl ;
            centerline_flux << "sub_case," << subCaseName << endl << endl ;
            centerline_flux << "x,x_over_L,radiative_flux_x_P1,radiative_flux_x_analytical" << endl ;

            for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << (key-westWallLocation)/plateSeparation << "," << val << "," << analyticalData.centerlineG.at(key) << endl ; }
            for(auto& [key,val]: p1.heatData()->centerlineFlux()){ centerline_flux << key << "," << (key-westWallLocation)/plateSeparation << "," << val << "," << analyticalData.centerlineFlux.at(key) << endl ; }
                
            centerline_G << endl ;
            centerline_flux << endl ;
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

        string parentFolder = "DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation" ;
        string caseName = "DOMRadiativeEnergyAndFluxProfileTest" ;
        string outputCaseName = caseName + "/ordinate" + to_string(ordinate) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
            
        // Lists the physical plate separation cases used in the study.
        vector<string> subCaseNames = {"300x1x1_1", "300x1x1_10"} ;
        
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 2.0 ;

        // Zero scattering.
        double sigmaS = 0 ;

        // Zero conductivity.
        double k = 0 ;

        // Emissivity values along the walls.
        vector<double> emissivity = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0} ; 

        // Both plates are maintained at 500 K.
        double Twest = 500 ;
        double Teast = 500 ;

        paths.createOutputDirectory(parentFolder, outputCaseName, "") ;

        // Diagnostics.
        ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder, outputCaseName, "") / "heat_diagnostics.csv") ;
        heat_flux_transfer_info << "west_wall_temperature,"
                                << "east_wall_temperature,"
                                << "runtime,"
                                << "absorption_coefficient,"
                                << "plate_separation,"
                                << "number_of_ordinates,"
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

        for(string subCaseName: subCaseNames){

            double plateSeparation = 0 ;

            if(subCaseName == "300x1x1_1"){
                plateSeparation = 1.0 ;
            }
            else if(subCaseName == "300x1x1_10"){
                plateSeparation = 10.0 ;
            }
            else{
                throw runtime_error("Unknown physical plate separation subcase.") ;
            }

            paths.createOutputDirectory(parentFolder, outputCaseName, subCaseName) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

            // Defines the field files for the initial G and T data.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

            // Constructs the DOM solver for the current case and executes the full simulation.
            string domCaseName = parentFolder + "/" + outputCaseName ;
            DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, subCaseName, domCaseName, subCaseName, ordinate) ;

            auto start = std::chrono::high_resolution_clock::now() ;
            dom.runDOMPureRadiationParticipatingMediumConstantT() ;
            auto end = std::chrono::high_resolution_clock::now() ;
            std::chrono::duration<double> elapsed = end-start ;

            dom.heatData()->calculateWestEastTransferForParallelPlates() ;
            dom.heatData()->calculateCenterlineFluxInX() ;

            ParallelPlateAnalyticalData analyticalData = calculateParallelPlateColdMediumAnalyticalValues(dom.heatData()->centerlineGInX(),plateSeparation,Twest,Teast,absorptionCoefficient) ;

            auto firstCenterlineCoordinate = dom.heatData()->centerlineGInX().begin() ;
            auto secondCenterlineCoordinate = next(firstCenterlineCoordinate) ;
            double centerlineCellSpacing = secondCenterlineCoordinate->first-firstCenterlineCoordinate->first ;
            double westWallLocation = firstCenterlineCoordinate->first-0.5*centerlineCellSpacing ;

            heat_flux_transfer_info << Twest << ","
                                    << Teast << ","
                                    << elapsed.count() << ","
                                    << absorptionCoefficient << ","
                                    << plateSeparation << ","
                                    << ordinate << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["east"] << ","
                                    << analyticalData.wallHeatFlux["east"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["west"] << ","
                                    << analyticalData.wallHeatFlux["west"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["north"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["south"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["top"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["bottom"] << ","
                                    << (dom.heatData()->averageWallNormalHeatFlux()["east"] + dom.heatData()->averageWallNormalHeatFlux()["west"] + dom.heatData()->averageWallNormalHeatFlux()["north"] + dom.heatData()->averageWallNormalHeatFlux()["south"] + dom.heatData()->averageWallNormalHeatFlux()["top"] + dom.heatData()->averageWallNormalHeatFlux()["bottom"]) << ","
                                    << (dom.heatData()->averageWallNormalHeatFlux()["east"] + dom.heatData()->averageWallNormalHeatFlux()["west"]) << ","
                                    << dom.heatData()->oppositeWallHeatFlux()["east"] << ","
                                    << dom.heatData()->oppositeWallHeatFlux()["west"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["east"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["west"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["north"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["south"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["top"] << ","
                                    << dom.heatData()->totalWallHeatTransfer()["bottom"] << ","
                                    << (dom.heatData()->totalWallHeatTransfer()["east"] + dom.heatData()->totalWallHeatTransfer()["west"] + dom.heatData()->totalWallHeatTransfer()["north"] + dom.heatData()->totalWallHeatTransfer()["south"] + dom.heatData()->totalWallHeatTransfer()["top"] + dom.heatData()->totalWallHeatTransfer()["bottom"]) << "," 
                                    << (dom.heatData()->totalWallHeatTransfer()["east"] + dom.heatData()->totalWallHeatTransfer()["west"]) << ","
                                    << dom.heatData()->oppositeWallHeatTransfer()["east"] << ","
                                    << dom.heatData()->oppositeWallHeatTransfer()["west"] << ","
                                    << dom.heatData()->dimensionlessQw()["east"] << ","
                                    << dom.heatData()->dimensionlessQw()["west"] << ","
                                    << dom.heatData()->dimensionlessQw()["north"] << ","
                                    << dom.heatData()->dimensionlessQw()["south"] << ","
                                    << dom.heatData()->dimensionlessQw()["top"] << ","
                                    << dom.heatData()->dimensionlessQw()["bottom"] << endl ;

            ofstream centerline_G(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "centerlineG.csv") ;
            ofstream centerline_flux(paths.getOutputDirectory(parentFolder, outputCaseName, subCaseName) / "centerlineFlux.csv") ;

            centerline_G << "west_wall_temperature," << Twest << endl ;
            centerline_G << "east_wall_temperature," << Teast << endl ;
            centerline_G << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_G << "plate_separation," << plateSeparation << endl ;
            centerline_G << "number_of_ordinates," << ordinate << endl ;
            centerline_G << "sub_case," << subCaseName << endl << endl ;
            centerline_G << "x,x_over_L,G_DOM,G_analytical" << endl ;

            centerline_flux << "west_wall_temperature," << Twest << endl ;
            centerline_flux << "east_wall_temperature," << Teast << endl ;
            centerline_flux << "absorption_coefficient," << absorptionCoefficient << endl ;
            centerline_flux << "plate_separation," << plateSeparation << endl ;
            centerline_flux << "number_of_ordinates," << ordinate << endl ;
            centerline_flux << "sub_case," << subCaseName << endl << endl ;
            centerline_flux << "x,x_over_L,radiative_flux_x_DOM,radiative_flux_x_analytical" << endl ;

            for(auto& [key,val]: dom.heatData()->centerlineGInX()){ centerline_G << key << "," << (key-westWallLocation)/plateSeparation << "," << val << "," << analyticalData.centerlineG.at(key) << endl ; }
            for(auto& [key,val]: dom.heatData()->centerlineFluxInX()){ centerline_flux << key << "," << (key-westWallLocation)/plateSeparation << "," << val << "," << analyticalData.centerlineFlux.at(key) << endl ; }
                
            centerline_G << endl ;
            centerline_flux << endl ;
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

        string parentFolder = "DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation" ;
        string caseName = "P1MeshConvergenceTest" ;

        paths.createOutputDirectory(parentFolder, caseName) ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
            
        // Lists the mesh sizes used in the convergence study.
        vector<int> meshSizes = {10,20,30,40,50,75,100,150,200,250,300} ;

        // Lists the physical plate separations used in the study.
        vector<double> plateSeparations = {1.0,10.0} ;
        
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 2.0 ;

        // Zero scattering.
        double sigmaS = 0 ;

        // Zero conductivity.
        double k = 0 ;

        // Emissivity values along the walls.
        vector<double> emissivity = {1.0,1.0,0.0,0.0,0.0,0.0} ; 

        // Both plates are maintained at 500 K.
        double Twest = 500 ;
        double Teast = 500 ;

        // Stores the convergence data in a directly plottable format.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, caseName) / "mesh_convergence.csv") ;
        convergence_info << "plate_separation,"
                        << "number_of_cells,"
                        << "cell_size,"
                        << "optical_cell_thickness,"
                        << "relative_L2_irradiation_error,"
                        << "observed_order,"
                        << "runtime" << endl ;

        for(double plateSeparation: plateSeparations){

            double previousError = 0 ;
            double previousCellSize = 0 ;
            bool firstMesh = true ;

            for(int meshSize: meshSizes){

                string distanceName ;

                if(plateSeparation == 1.0){
                    distanceName = "1" ;
                }
                else if(plateSeparation == 10.0){
                    distanceName = "10" ;
                }
                else{
                    throw runtime_error("Unknown physical plate separation.") ;
                }

                string subCaseName = to_string(meshSize) + "x1x1_" + distanceName ;

                paths.createOutputDirectory(parentFolder, caseName, subCaseName) ;

                // Retrieves the mesh files for the current mesh configuration.
                meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

                // Defines the field files for the initial G and T data.
                fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

                // Constructs the P1 solver for the current case and executes the full simulation.
                string p1CaseName = parentFolder + "/" + caseName ;
                P1 p1(meshFiles, fieldFiles, absorptionCoefficient, sigmaS, k, emissivity, subCaseName, p1CaseName, subCaseName) ;

                auto start = std::chrono::high_resolution_clock::now() ;
                p1.runP1PureRadiationParticipatingMediumTConstant() ;
                auto end = std::chrono::high_resolution_clock::now() ;
                std::chrono::duration<double> elapsed = end-start ;

                p1.heatData()->calculateCenterlineTemperatureAndG() ;

                ParallelPlateAnalyticalData analyticalData = calculateParallelPlateColdMediumAnalyticalValues(p1.heatData()->centerlineG(),plateSeparation,Twest,Teast,absorptionCoefficient) ;

                double relativeL2Error = calculateCenterlineRelativeL2Error(p1.heatData()->centerlineG(),analyticalData.centerlineG) ;
                double cellSize = plateSeparation/meshSize ;
                double opticalCellThickness = absorptionCoefficient*cellSize ;

                convergence_info << plateSeparation << ","
                                << meshSize << ","
                                << cellSize << ","
                                << opticalCellThickness << ","
                                << relativeL2Error << "," ;

                if(firstMesh){
                    convergence_info << "," ;
                    firstMesh = false ;
                }
                else{
                    double observedOrder = log(previousError/relativeL2Error)/log(previousCellSize/cellSize) ;
                    convergence_info << observedOrder << "," ;
                }

                convergence_info << elapsed.count() << endl ;

                previousError = relativeL2Error ;
                previousCellSize = cellSize ;
            }
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

        string parentFolder = "DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation" ;
        string caseName = "DOMMeshConvergenceTest" ;
        string outputCaseName = caseName + "/ordinate" + to_string(ordinate) ;

        paths.createOutputDirectory(parentFolder, outputCaseName, "") ;

        // Stores the mesh and field file lists used to initialise each simulation case.
        vector<string> meshFiles, fieldFiles ;
            
        // Lists the mesh sizes used in the convergence study.
        vector<int> meshSizes = {10,20,30,40,50,75,100,150,200,250,300} ;

        // Lists the physical plate separations used in the study.
        vector<double> plateSeparations = {1.0,10.0} ;
        
        // Absorption coefficient of the participating medium.
        double absorptionCoefficient = 2.0 ;

        // Zero scattering.
        double sigmaS = 0 ;

        // Zero conductivity.
        double k = 0 ;

        // Emissivity values along the walls.
        vector<double> emissivity = {1.0,1.0,0.0,0.0,0.0,0.0} ; 

        // Both plates are maintained at 500 K.
        double Twest = 500 ;
        double Teast = 500 ;

        // Stores the convergence data in a directly plottable format.
        ofstream convergence_info(paths.getOutputDirectory(parentFolder, outputCaseName, "") / "mesh_convergence.csv") ;
        convergence_info << "number_of_ordinates,"
                        << "plate_separation,"
                        << "number_of_cells,"
                        << "cell_size,"
                        << "optical_cell_thickness,"
                        << "relative_L2_irradiation_error,"
                        << "observed_order,"
                        << "runtime" << endl ;

        for(double plateSeparation: plateSeparations){

            double previousError = 0 ;
            double previousCellSize = 0 ;
            bool firstMesh = true ;

            for(int meshSize: meshSizes){

                string distanceName ;

                if(plateSeparation == 1.0){
                    distanceName = "1" ;
                }
                else if(plateSeparation == 10.0){
                    distanceName = "10" ;
                }
                else{
                    throw runtime_error("Unknown physical plate separation.") ;
                }

                string subCaseName = to_string(meshSize) + "x1x1_" + distanceName ;

                paths.createOutputDirectory(parentFolder, outputCaseName, subCaseName) ;

                // Retrieves the mesh files for the current mesh configuration.
                meshFiles = paths.getMeshFilenames(parentFolder, caseName, subCaseName) ;

                // Defines the field files for the initial G and T data.
                fieldFiles = {paths.getTFieldFile(parentFolder, caseName), paths.getGFieldFile(parentFolder, caseName)} ;

                // Constructs the DOM solver for the current case and executes the full simulation.
                string domCaseName = parentFolder + "/" + outputCaseName ;
                DOM dom(meshFiles, fieldFiles, absorptionCoefficient, k, sigmaS, emissivity, subCaseName, domCaseName, subCaseName, ordinate) ;

                auto start = std::chrono::high_resolution_clock::now() ;
                dom.runDOMPureRadiationParticipatingMediumConstantT() ;
                auto end = std::chrono::high_resolution_clock::now() ;
                std::chrono::duration<double> elapsed = end-start ;

                ParallelPlateAnalyticalData analyticalData = calculateParallelPlateColdMediumAnalyticalValues(dom.heatData()->centerlineGInX(),plateSeparation,Twest,Teast,absorptionCoefficient) ;

                double relativeL2Error = calculateCenterlineRelativeL2Error(dom.heatData()->centerlineGInX(),analyticalData.centerlineG) ;
                double cellSize = plateSeparation/meshSize ;
                double opticalCellThickness = absorptionCoefficient*cellSize ;

                convergence_info << ordinate << ","
                                << plateSeparation << ","
                                << meshSize << ","
                                << cellSize << ","
                                << opticalCellThickness << ","
                                << relativeL2Error << "," ;

                if(firstMesh){
                    convergence_info << "," ;
                    firstMesh = false ;
                }
                else{
                    double observedOrder = log(previousError/relativeL2Error)/log(previousCellSize/cellSize) ;
                    convergence_info << observedOrder << "," ;
                }

                convergence_info << elapsed.count() << endl ;

                previousError = relativeL2Error ;
                previousCellSize = cellSize ;
            }
        }
    }

    ParallelPlateAnalyticalData calculateParallelPlateColdMediumAnalyticalValues(const map<double,double>& centerlineCoordinates, double plateSeparation, double Twest, double Teast, double absorptionCoefficient){

        const array<double,20> gaussPoints = {
            -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
            -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154196,
            -0.2277858511416451, -0.0765265211334973, 0.0765265211334973, 0.2277858511416451,
            0.3737060887154196, 0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
            0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949
        } ;

        const array<double,20> gaussWeights = {
            0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767048,
            0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
            0.1491729864726037, 0.1527533871307258, 0.1527533871307258, 0.1491729864726037,
            0.1420961093183821, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
            0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
        } ;

        ParallelPlateAnalyticalData analyticalData ;

        if(centerlineCoordinates.size() < 2){ throw runtime_error("At least two centreline coordinates are required for the analytical solution.") ;}

        double sigma = 5.670374419e-8 ;
        double westWallIntensity = sigma*powl(Twest,4)/M_PI ;
        double eastWallIntensity = sigma*powl(Teast,4)/M_PI ;
        auto firstCoordinate = centerlineCoordinates.begin() ;
        auto secondCoordinate = next(firstCoordinate) ;
        double cellSpacing = secondCoordinate->first-firstCoordinate->first ;
        double westWallLocation = firstCoordinate->first-0.5*cellSpacing ;

        for(const auto& [x,val]: centerlineCoordinates){

            double analyticalG = 0 ;
            double analyticalFlux = 0 ;
            double distanceFromWestWall = x-westWallLocation ;
            double distanceFromEastWall = plateSeparation-distanceFromWestWall ;

            for(int ordinate=0; ordinate<20; ++ordinate){

                double mu = 0.5*(gaussPoints[ordinate]+1.0) ;
                double weight = 0.5*gaussWeights[ordinate] ;

                double westIntensity = westWallIntensity*exp(-absorptionCoefficient*distanceFromWestWall/mu) ;
                double eastIntensity = eastWallIntensity*exp(-absorptionCoefficient*distanceFromEastWall/mu) ;

                analyticalG += 2.0*M_PI*weight*(westIntensity+eastIntensity) ;
                analyticalFlux += 2.0*M_PI*weight*mu*(westIntensity-eastIntensity) ;
            }

            analyticalData.centerlineG[x] = analyticalG ;
            analyticalData.centerlineFlux[x] = analyticalFlux ;
        }

        double westWallHeatFlux = 0 ;
        double eastWallHeatFlux = 0 ;

        for(int ordinate=0; ordinate<20; ++ordinate){

            double mu = 0.5*(gaussPoints[ordinate]+1.0) ;
            double weight = 0.5*gaussWeights[ordinate] ;

            double eastIntensityAtWestWall = eastWallIntensity*exp(-absorptionCoefficient*plateSeparation/mu) ;
            double westIntensityAtEastWall = westWallIntensity*exp(-absorptionCoefficient*plateSeparation/mu) ;

            westWallHeatFlux += 2.0*M_PI*weight*mu*(westWallIntensity-eastIntensityAtWestWall) ;
            eastWallHeatFlux += 2.0*M_PI*weight*mu*(eastWallIntensity-westIntensityAtEastWall) ;
        }

        analyticalData.wallHeatFlux["west"] = westWallHeatFlux ;
        analyticalData.wallHeatFlux["east"] = eastWallHeatFlux ;

        return analyticalData ;
    }

    double calculateCenterlineRelativeL2Error(const map<double,double>& numericalCenterline, const map<double,double>& analyticalCenterline){

        if(numericalCenterline.size() != analyticalCenterline.size()){
            throw runtime_error("Numerical and analytical centerline data sizes do not match.") ;
        }

        double errorSum = 0 ;
        double analyticalSum = 0 ;

        for(const auto& [coordinate,numericalValue]: numericalCenterline){

            auto analyticalValue = analyticalCenterline.find(coordinate) ;

            if(analyticalValue == analyticalCenterline.end()){
                throw runtime_error("Numerical and analytical centerline coordinates do not match.") ;
            }

            errorSum += powl(numericalValue-analyticalValue->second,2) ;
            analyticalSum += powl(analyticalValue->second,2) ;
        }

        if(analyticalSum == 0){
            throw runtime_error("Analytical centerline norm is zero.") ;
        }

        return sqrt(errorSum/analyticalSum) ;
    }


} ;

