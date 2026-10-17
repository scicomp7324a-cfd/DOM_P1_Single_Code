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
# include <algorithm>
# include <cerrno>
# include <chrono>
# include <cstdint>
# include <iomanip>
# include <sstream>
# include <sys/resource.h>
# include <sys/types.h>
# include <sys/wait.h>
# include <unistd.h>

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

# include "AnalyticalEvaluator.hpp"

using namespace std ;
namespace fs = std::filesystem;

class DOMP1ComparisionFiniteParallelPlateMediumOpen{

    public:

        DOMP1ComparisionFiniteParallelPlateMediumOpen(){}

        void test(){
            unifiedP1DOMMemoryAccuracyExperiment() ;
        }

        void p1ParallelPlateSimulationPureRadiationNoParticipatingMediumOpen(){

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

            string caseName = "ParallelPlateNoMediumVacuumOpen" ;
            //vector<string> subCaseNames = {"300"} ; // East Wall Temperature is taken here as sub case name as it is varying

            paths.createOutputDirectory(caseName) ;

            // Stores the mesh and field file lists used to initialise each simulation case.
            vector<string> meshFiles, fieldFiles ;
            
            // Lists the mesh sizes used in the study.
            // string dimensions = "20x5x5_0.01x1x1" ;
            // vector<string> subCaseNames = {"10x5x5", "20x5x5", "40x5x5"} ;
            vector<string> subCaseNames = {"40x5x5"} ;

            // Absorption coefficient value. No participating medium
            double a = 0.0 ;

            // Zero scattering
            double sigmaS = 0 ;

            // Zero conductivity
            double k = 0 ;

            // Emissivity values along the walls
            vector<double> emissivity = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0} ; 

            // East Wall temperatures 
            int Twest = 1000 ;
            int Teast = 300 ;

            // Diagnostics 
            ofstream heat_flux_transfer_info(paths.getOutputDirectory(caseName) / "heat_diagnostics.csv") ;
            heat_flux_transfer_info << "Mesh,"
                                    << "west_wall_temperature,"
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
                                    << "opposite_wall_heat_transfer_west," 
                                    << "heat_rate_analytical_A,"
                                    << "heat_rate_analytical_B,"
                                    << "heat_loss_surrounding_analytical_A,"
                                    << "heat_loss_surrounding_analytical_B,"
                                    << "heat_loss_flux_surrounding_analytical_A,"
                                    << "heat_loss_flux_surrounding_analytical_B,"
                                    << "gross_heat_transfer_AB_analytical,"
                                    << "gross_heat_transfer_BA_analytical,"
                                    << "net_heat_transfer_analytical,"
                                    << "error_heat_flux_plate_A,"
                                    << "error_heat_rate_plate_A," 
                                    << "energy_imbalance" << endl ;


            AnalyticalEvaluator analyticalValue(1, 0.01, 1, 1, 1.0, 1.0) ;
            
            for(string subCaseName: subCaseNames){

                analyticalValue.calculateHeatFluxAndTransfer(Twest, Teast) ;
                paths.createOutputDirectory(caseName, subCaseName) ;

                // Retrieves the mesh files for the current mesh configuration.
                meshFiles = paths.getMeshFilenames(caseName, subCaseName) ;
                cout << "reached here" << endl ;
                // Defines the field files for the initial G and T data.
                fieldFiles = {paths.getTFieldFile(caseName, to_string(Teast)), paths.getGFieldFile(caseName, to_string(Teast))} ;
                cout << "reached here too" << endl ;
                cout << fieldFiles[0] << endl ;
                cout << fieldFiles[1] << endl ;

                // Constructs the SIMPLE solver for the current case and executes the full simulation.
                P1 p1(meshFiles, fieldFiles, a, sigmaS, k, emissivity, subCaseName, caseName, subCaseName) ;
                auto start = std::chrono::high_resolution_clock::now();
                p1.runP1PureRadiationNoParticipatingMedium() ;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;

                p1.heatData()->calculateWestEastTransferForParallelPlates() ;
                //map<string,double> GAnalytical = p1.heatData()->parallelPlateHeatFluxBlackbodyAnalyticalValue(Twest,stod(subCaseName)) ;
                double error_heat_flux_plate_A = abs(p1.heatData()->avgqw()["west"] - analyticalValue.q1Total())/abs(analyticalValue.q1Total()) ;
                double error_heat_rate_plate_A = abs(p1.heatData()->totalQw()["west"] - analyticalValue.Q1Total())/abs(analyticalValue.Q1Total()) ;
                double energy_imbalance = abs(p1.heatData()->totalQw()["east"] + p1.heatData()->totalQw()["west"] + 
                                              p1.heatData()->totalQw()["north"] + p1.heatData()->totalQw()["south"] + 
                                              p1.heatData()->totalQw()["top"] + p1.heatData()->totalQw()["bottom"])/
                                              max(
                                                max(abs(p1.heatData()->totalQw()["east"]),abs(p1.heatData()->totalQw()["west"])), 
                                                abs(p1.heatData()->totalQw()["north"] + p1.heatData()->totalQw()["south"] + p1.heatData()->totalQw()["top"] + p1.heatData()->totalQw()["bottom"])
                                                ) ;

                heat_flux_transfer_info << subCaseName << ","
                                        << Twest << ","
                                        << Teast << ","
                                        << elapsed.count() << ","
                                        << p1.heatData()->avgqw()["east"] << ","
                                        << analyticalValue.q2Total() << ","
                                        << p1.heatData()->avgqw()["west"] << ","
                                        << analyticalValue.q1Total() << ","
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
                                        << analyticalValue.Q1Total() << ","
                                        << analyticalValue.Q2Total() << ","
                                        << analyticalValue.Q1Surrounding() << "," 
                                        << analyticalValue.Q2Surrounding()<< "," 
                                        << analyticalValue.q1Surrounding() << "," 
                                        << analyticalValue.q2Surrounding() << "," 
                                        << analyticalValue.Q12Gross() << "," 
                                        << analyticalValue.Q21Gross() << "," 
                                        << analyticalValue.Q12Net() << "," 
                                        << error_heat_flux_plate_A << "," 
                                        << error_heat_rate_plate_A << ","
                                        << energy_imbalance << endl ;


                ofstream centerline_G(paths.getOutputDirectory(caseName, subCaseName) / "centerlineG.csv") ;
                centerline_G << "west_wall_temperature," << Twest << endl ;
                centerline_G << "east_wall_temperature," << Teast << endl << endl ;
                centerline_G << "x_centerline_coordinate, G" << endl ;
                for(auto& [key,val]: p1.heatData()->centerlineG()){ centerline_G << key << "," << val << endl ; }
                centerline_G << endl ;

            }

        }

        void domParallelPlateSimulationPureRadiationNoParticipatingMediumOpen(){

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

            string caseName = "ParallelPlateNoMediumVacuumOpenDOM" ;
            //vector<string> subCaseNames = {"300","400","500","600","700","800","900","1000"} ; // East Wall Temperature is taken here as sub case name as it is varying
            //vector<string> subCaseNames = {"300"} ; // East Wall Temperature is taken here as sub case name as it is varying

            paths.createOutputDirectory(caseName) ;

            // Stores the mesh and field file lists used to initialise each simulation case.
            vector<string> meshFiles, fieldFiles ;
            
            // Lists the mesh sizes used in the study.
            //string dimensions = "10x10x10" ;
            vector<string> subCaseNames = {"10x5x5", "20x5x5", "40x5x5"} ;

            // Absorption coefficient value. No participating medium
            double a = 0.0 ;

            // Zero scattering
            double sigmaS = 0 ;

            // Zero conductivity
            double k = 0 ;

            // Emissivity values along the walls
            vector<double> emissivity = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0} ; 

            // East Wall temperatures 
            int Twest = 1000 ;
            int Teast = 300 ;

            int ordinate = 8 ;

            // Diagnostics 
            ofstream heat_flux_transfer_info(paths.getOutputDirectory(caseName) / "heat_diagnostics.csv") ;
            heat_flux_transfer_info << "Mesh,"
                                    << "Ordinate,"
                                    << "west_wall_temperature,"
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
                                    << "opposite_wall_heat_transfer_west,"
                                    << "heat_rate_analytical_A,"
                                    << "heat_rate_analytical_B,"
                                    << "heat_loss_surrounding_analytical_A,"
                                    << "heat_loss_surrounding_analytical_B,"
                                    << "heat_loss_flux_surrounding_analytical_A,"
                                    << "heat_loss_flux_surrounding_analytical_B,"
                                    << "gross_heat_transfer_AB_analytical,"
                                    << "gross_heat_transfer_BA_analytical,"
                                    << "net_heat_transfer_analytical,"
                                    << "error_heat_flux_plate_A,"
                                    << "error_heat_rate_plate_A," 
                                    << "energy_imbalance" << endl ;


            AnalyticalEvaluator analyticalValue(1, 0.01, 1, 1, 1.0, 1.0) ;

            for(string subCaseName: subCaseNames){

                analyticalValue.calculateHeatFluxAndTransfer(Twest, Teast) ;
                paths.createOutputDirectory(caseName, subCaseName) ;

                // Retrieves the mesh files for the current mesh configuration.
                meshFiles = paths.getMeshFilenames(caseName, subCaseName) ;

                // Defines the field files for the initial G and T data.
                fieldFiles = {paths.getTFieldFile(caseName, to_string(Teast))} ;

                // Constructs the SIMPLE solver for the current case and executes the full simulation.
                DOM dom(meshFiles, fieldFiles, a, k, sigmaS, emissivity, subCaseName, caseName, subCaseName, ordinate) ;
                auto start = std::chrono::high_resolution_clock::now();
                dom.runDOMPureRadiationNoParticipatingMedium() ;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;

                dom.heatData()->calculateWestEastTransferForParallelPlates() ;
                //map<string,double> GAnalytical = dom.heatData()->parallelPlateHeatFluxBlackbodyAnalyticalValue(Twest,Teast) ;

                double error_heat_flux_plate_A = abs(dom.heatData()->averageWallNormalHeatFlux()["west"] - analyticalValue.q1Total())/abs(analyticalValue.q1Total()) ;
                double error_heat_rate_plate_A = abs(dom.heatData()->totalWallHeatTransfer()["west"] - analyticalValue.Q1Total())/abs(analyticalValue.Q1Total()) ;
                double energy_imbalance = abs(dom.heatData()->totalWallHeatTransfer()["east"] + dom.heatData()->totalWallHeatTransfer()["west"] +
                                              dom.heatData()->totalWallHeatTransfer()["north"] + dom.heatData()->totalWallHeatTransfer()["south"] +
                                              dom.heatData()->totalWallHeatTransfer()["top"] + dom.heatData()->totalWallHeatTransfer()["bottom"])/
                                              max(abs(dom.heatData()->totalWallHeatTransfer()["north"] + dom.heatData()->totalWallHeatTransfer()["south"] +
                                              dom.heatData()->totalWallHeatTransfer()["top"] + dom.heatData()->totalWallHeatTransfer()["bottom"]),
                                              max(
                                                abs(dom.heatData()->totalWallHeatTransfer()["east"]),
                                                abs(dom.heatData()->totalWallHeatTransfer()["west"]))
                                                ) ;

                heat_flux_transfer_info << subCaseName << ","
                                        << ordinate << ","
                                        << Twest << ","
                                        << Teast << ","
                                        << elapsed.count() << ","
                                        << dom.heatData()->averageWallNormalHeatFlux()["east"] << ","
                                        << analyticalValue.q2Total() << ","
                                        << dom.heatData()->averageWallNormalHeatFlux()["west"] << ","
                                        << analyticalValue.q1Total() << ","
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
                                        << analyticalValue.Q1Total() << ","
                                        << analyticalValue.Q2Total() << ","
                                        << analyticalValue.Q1Surrounding() << "," 
                                        << analyticalValue.Q2Surrounding()<< "," 
                                        << analyticalValue.q1Surrounding() << "," 
                                        << analyticalValue.q2Surrounding() << "," 
                                        << analyticalValue.Q12Gross() << "," 
                                        << analyticalValue.Q21Gross() << "," 
                                        << analyticalValue.Q12Net() << "," 
                                        << error_heat_flux_plate_A << "," 
                                        << error_heat_rate_plate_A << ","
                                        << energy_imbalance << endl ;

                ofstream centerline_G(paths.getOutputDirectory(caseName, subCaseName) / "centerlineGInX.csv") ;
                centerline_G << "west_wall_temperature," << Twest << endl ;
                centerline_G << "east_wall_temperature," << Teast << endl << endl ;
                centerline_G << "x_centerline_coordinate, G" << endl ;
                for(auto& [key,val]: dom.heatData()->centerlineGInX()){ centerline_G << key << "," << val << endl ; }
                centerline_G << endl ;

            }

        }

        void unifiedP1DOMMemoryAccuracyExperiment(){

            /*
                Runs every P1 and DOM configuration in an independent child
                process. This is necessary because ru_maxrss is cumulative
                within one process and cannot be reset between simulations.

                The existing P1.hpp, DOM.hpp and all other codebase files are
                used without modification.
            */

            struct ProcessResult{
                int success = 0 ;
                long long numberOfCells = 0 ;
                int ordinateOrder = 0 ;
                long long numberOfOrdinates = 0 ;
                long long resolution = 0 ;
                double plateAHeatRate = 0.0 ;
                double analyticalPlateAHeatRate = 0.0 ;
                double relativeErrorPlateA = 0.0 ;
                double energyImbalance = 0.0 ;
                double peakRSSMiB = 0.0 ;
            } ;

            struct NamedResult{
                string processName ;
                string configuration ;
                string resolutionDefinition ;
                ProcessResult values ;
            } ;

            auto getPeakRSSMiB = []() -> double{

                struct rusage usage{} ;

                if(getrusage(RUSAGE_SELF, &usage) != 0){
                    throw runtime_error("getrusage failed while obtaining peak RSS.") ;
                }

                #ifdef __APPLE__
                    return static_cast<double>(usage.ru_maxrss)/(1024.0*1024.0) ;
                #else
                    return static_cast<double>(usage.ru_maxrss)/1024.0 ;
                #endif
            } ;

            auto getNumberOfCells = [](const string& meshName) -> long long{

                long long nx = 0 ;
                long long ny = 0 ;
                long long nz = 0 ;
                char separator1 = '\0' ;
                char separator2 = '\0' ;

                stringstream parser(meshName) ;
                parser >> nx >> separator1 >> ny >> separator2 >> nz ;

                if(!parser || separator1 != 'x' || separator2 != 'x' ||
                   nx <= 0 || ny <= 0 || nz <= 0){
                    throw runtime_error("Invalid mesh configuration name: " + meshName) ;
                }

                return nx*ny*nz ;
            } ;

            auto writeAll = [](int fileDescriptor, const void* data, size_t bytes) -> bool{

                const char* current = static_cast<const char*>(data) ;

                while(bytes > 0){

                    const ssize_t written = ::write(fileDescriptor, current, bytes) ;

                    if(written < 0){
                        if(errno == EINTR){
                            continue ;
                        }
                        return false ;
                    }

                    current += written ;
                    bytes -= static_cast<size_t>(written) ;
                }

                return true ;
            } ;

            auto readAll = [](int fileDescriptor, void* data, size_t bytes) -> bool{

                char* current = static_cast<char*>(data) ;

                while(bytes > 0){

                    const ssize_t received = ::read(fileDescriptor, current, bytes) ;

                    if(received == 0){
                        return false ;
                    }

                    if(received < 0){
                        if(errno == EINTR){
                            continue ;
                        }
                        return false ;
                    }

                    current += received ;
                    bytes -= static_cast<size_t>(received) ;
                }

                return true ;
            } ;

            PathConfig paths ;

            try{
                paths.validate() ;
                paths.ensureOutputDir() ;
            }
            catch(const exception& e){
                cerr << "Path error: " << e.what() << endl ;
                throw ;
            }

            const string p1CaseName = "ParallelPlateNoMediumVacuumOpen" ;
            const string domCaseName = "ParallelPlateNoMediumVacuumOpenDOM" ;
            const string comparisonOutputName = "DOMP1FiniteParallelPlateMemoryAccuracy" ;

            const vector<string> meshConfigurations = {
                "10x10x10",
                "20x20x20",
                "40x40x40"
            } ;

            const double a = 0.0 ;
            const double sigmaS = 0.0 ;
            const double k = 0.0 ;
            const vector<double> emissivity = {
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            } ;

            const int Twest = 1000 ;
            const int Teast = 300 ;
            const vector<int> ordinateOrders = {
                2, 4, 6, 8
            } ;

            paths.createOutputDirectory(comparisonOutputName) ;

            auto runOneCase = [&](const string& processName,
                                  const string& meshConfiguration,
                                  int requestedOrdinateOrder) -> ProcessResult{

                int pipeDescriptor[2] ;

                if(::pipe(pipeDescriptor) != 0){
                    throw runtime_error("Failed to create pipe for isolated memory measurement.") ;
                }

                cout.flush() ;
                cerr.flush() ;

                const pid_t processId = ::fork() ;

                if(processId < 0){
                    ::close(pipeDescriptor[0]) ;
                    ::close(pipeDescriptor[1]) ;
                    throw runtime_error("Failed to fork an isolated simulation process.") ;
                }

                if(processId == 0){

                    ::close(pipeDescriptor[0]) ;

                    ProcessResult result{} ;

                    try{
                        AnalyticalEvaluator analyticalValue(
                            1.0, 0.01, 1.0, 1.0, 1.0, 1.0
                        ) ;

                        analyticalValue.calculateHeatFluxAndTransfer(
                            Twest, Teast
                        ) ;

                        result.numberOfCells =
                            getNumberOfCells(meshConfiguration) ;

                        if(processName == "P1"){

                            paths.createOutputDirectory(
                                p1CaseName, meshConfiguration
                            ) ;

                            const vector<string> meshFiles =
                                paths.getMeshFilenames(
                                    p1CaseName, meshConfiguration
                                ) ;

                            const vector<string> fieldFiles = {
                                paths.getTFieldFile(
                                    p1CaseName, to_string(Teast)
                                ),
                                paths.getGFieldFile(
                                    p1CaseName, to_string(Teast)
                                )
                            } ;

                            P1 p1(
                                meshFiles,
                                fieldFiles,
                                a,
                                sigmaS,
                                k,
                                emissivity,
                                meshConfiguration,
                                p1CaseName,
                                meshConfiguration
                            ) ;

                            p1.runP1PureRadiationNoParticipatingMedium() ;
                            p1.heatData()->
                                calculateWestEastTransferForParallelPlates() ;

                            const auto& wallHeatRate =
                                p1.heatData()->totalQw() ;

                            const double east = wallHeatRate.at("east") ;
                            const double west = wallHeatRate.at("west") ;
                            const double north = wallHeatRate.at("north") ;
                            const double south = wallHeatRate.at("south") ;
                            const double top = wallHeatRate.at("top") ;
                            const double bottom = wallHeatRate.at("bottom") ;

                            const double surrounding =
                                north + south + top + bottom ;

                            const double total =
                                east + west + surrounding ;

                            const double normalisingHeatRate = max(
                                max(abs(east), abs(west)),
                                abs(surrounding)
                            ) ;

                            result.ordinateOrder = 0 ;
                            result.numberOfOrdinates = 0 ;
                            result.resolution = result.numberOfCells ;
                            result.plateAHeatRate = west ;
                            result.analyticalPlateAHeatRate =
                                analyticalValue.Q1Total() ;
                            result.relativeErrorPlateA =
                                abs(
                                    result.plateAHeatRate -
                                    result.analyticalPlateAHeatRate
                                )/
                                abs(result.analyticalPlateAHeatRate) ;
                            result.energyImbalance =
                                normalisingHeatRate > 0.0 ?
                                abs(total)/normalisingHeatRate : 0.0 ;
                        }
                        else if(processName == "DOM"){

                            const string domOutputConfiguration =
                                meshConfiguration + "_S" +
                                to_string(requestedOrdinateOrder) ;

                            paths.createOutputDirectory(
                                domCaseName, domOutputConfiguration
                            ) ;

                            const vector<string> meshFiles =
                                paths.getMeshFilenames(
                                    domCaseName, meshConfiguration
                                ) ;

                            const vector<string> fieldFiles = {
                                paths.getTFieldFile(
                                    domCaseName, to_string(Teast)
                                )
                            } ;

                            DOM dom(
                                meshFiles,
                                fieldFiles,
                                a,
                                k,
                                sigmaS,
                                emissivity,
                                meshConfiguration,
                                domCaseName,
                                domOutputConfiguration,
                                requestedOrdinateOrder
                            ) ;

                            dom.runDOMPureRadiationNoParticipatingMedium() ;
                            dom.heatData()->
                                calculateWestEastTransferForParallelPlates() ;

                            const auto& wallHeatRate =
                                dom.heatData()->totalWallHeatTransfer() ;

                            const double east = wallHeatRate.at("east") ;
                            const double west = wallHeatRate.at("west") ;
                            const double north = wallHeatRate.at("north") ;
                            const double south = wallHeatRate.at("south") ;
                            const double top = wallHeatRate.at("top") ;
                            const double bottom = wallHeatRate.at("bottom") ;

                            const double surrounding =
                                north + south + top + bottom ;

                            const double total =
                                east + west + surrounding ;

                            const double normalisingHeatRate = max(
                                max(abs(east), abs(west)),
                                abs(surrounding)
                            ) ;

                            result.ordinateOrder =
                                requestedOrdinateOrder ;
                            result.numberOfOrdinates =
                                static_cast<long long>(
                                    requestedOrdinateOrder
                                )*
                                static_cast<long long>(
                                    requestedOrdinateOrder + 2
                                ) ;
                            result.resolution =
                                result.numberOfCells*
                                result.numberOfOrdinates ;
                            result.plateAHeatRate = west ;
                            result.analyticalPlateAHeatRate =
                                analyticalValue.Q1Total() ;
                            result.relativeErrorPlateA =
                                abs(
                                    result.plateAHeatRate -
                                    result.analyticalPlateAHeatRate
                                )/
                                abs(result.analyticalPlateAHeatRate) ;
                            result.energyImbalance =
                                normalisingHeatRate > 0.0 ?
                                abs(total)/normalisingHeatRate : 0.0 ;
                        }
                        else{
                            throw runtime_error(
                                "Unsupported process name: " + processName
                            ) ;
                        }

                        result.peakRSSMiB = getPeakRSSMiB() ;
                        result.success = 1 ;
                    }
                    catch(const exception& e){
                        cerr << processName << " "
                             << meshConfiguration ;

                        if(processName == "DOM"){
                            cerr << "_S" << requestedOrdinateOrder ;
                        }

                        cerr << " failed: " << e.what() << endl ;
                        result.success = 0 ;
                    }

                    const bool sent = writeAll(
                        pipeDescriptor[1],
                        &result,
                        sizeof(ProcessResult)
                    ) ;

                    ::close(pipeDescriptor[1]) ;
                    cout.flush() ;
                    cerr.flush() ;
                    ::_exit(sent && result.success == 1 ? 0 : 1) ;
                }

                ::close(pipeDescriptor[1]) ;

                ProcessResult result{} ;
                const bool received = readAll(
                    pipeDescriptor[0],
                    &result,
                    sizeof(ProcessResult)
                ) ;

                ::close(pipeDescriptor[0]) ;

                int processStatus = 0 ;

                if(::waitpid(processId, &processStatus, 0) < 0){
                    throw runtime_error("waitpid failed for isolated simulation.") ;
                }

                if(!received ||
                   !WIFEXITED(processStatus) ||
                   WEXITSTATUS(processStatus) != 0 ||
                   result.success != 1){
                    throw runtime_error(
                        processName + " simulation failed for " +
                        meshConfiguration
                    ) ;
                }

                return result ;
            } ;

            vector<NamedResult> results ;

            for(const string& meshConfiguration: meshConfigurations){
                results.push_back({
                    "P1",
                    meshConfiguration,
                    "N_cells",
                    runOneCase("P1", meshConfiguration, 0)
                }) ;
            }

            for(const string& meshConfiguration: meshConfigurations){

                for(const int ordinateOrder: ordinateOrders){

                    results.push_back({
                        "DOM",
                        meshConfiguration + "_S" +
                            to_string(ordinateOrder),
                        "N_cells_x_N_ordinates",
                        runOneCase(
                            "DOM",
                            meshConfiguration,
                            ordinateOrder
                        )
                    }) ;
                }
            }

            const fs::path csvPath =
                paths.getOutputDirectory(comparisonOutputName)/
                "unified_memory_accuracy_domp1.csv" ;

            ofstream output(csvPath) ;

            if(!output){
                throw runtime_error(
                    "Failed to open unified CSV file: " +
                    csvPath.string()
                ) ;
            }

            output
                << "process_name,"
                << "configuration,"
                << "resolution_definition,"
                << "resolution,"
                << "number_of_cells,"
                << "ordinate_order,"
                << "number_of_ordinates,"
                << "peak_RSS_MiB,"
                << "plate_A_heat_rate_W,"
                << "plate_A_analytical_heat_rate_W,"
                << "relative_error_plate_A,"
                << "energy_imbalance\n" ;

            output << setprecision(15) ;

            for(const NamedResult& result: results){
                output
                    << result.processName << ","
                    << result.configuration << ","
                    << result.resolutionDefinition << ","
                    << result.values.resolution << ","
                    << result.values.numberOfCells << ","
                    << result.values.ordinateOrder << ","
                    << result.values.numberOfOrdinates << ","
                    << result.values.peakRSSMiB << ","
                    << result.values.plateAHeatRate << ","
                    << result.values.analyticalPlateAHeatRate << ","
                    << result.values.relativeErrorPlateA << ","
                    << result.values.energyImbalance
                    << "\n" ;
            }

            cout << "Unified P1/DOM comparison written to: "
                 << csvPath << endl ;
        }



} ;