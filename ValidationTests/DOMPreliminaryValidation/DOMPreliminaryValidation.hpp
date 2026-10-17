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

class DOMPreliminaryValidation{

    public:

    void test(){
        //domParallelPlateSimulationPureRadiationNoParticipatingMediumClosed() ;
        domParallelPlateSimulationPureRadiationParticipatingMedium() ;
        domSquarePlateSimulationRadiationConduction() ;
    }

    // Parallel plate simulation of pure radiation with no participating medium
    void domParallelPlateSimulationPureRadiationNoParticipatingMediumClosed(){

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

    string parentFolder = "DOMPreliminaryValidation" ;
    string caseName = "ParallelPlateNoMediumVacuumClosedDOM" ;
    vector<string> subCaseNames = {"300","400","500","600","700","800","900","1000"} ; // East Wall Temperature is taken here as sub case name as it is varying
    //vector<string> subCaseNames = {"300"} ; // East Wall Temperature is taken here as sub case name as it is varying

    paths.createOutputDirectory(parentFolder, caseName) ;

    // Stores the mesh and field file lists used to initialise each simulation case.
    vector<string> meshFiles, fieldFiles ;
    
    // Lists the mesh sizes used in the study.
    string dimensions = "300x1x1" ;

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
    vector<int> ordinates = {2,4,6,8} ;

    // Diagnostics 
    ofstream heat_flux_transfer_info(paths.getOutputDirectory(parentFolder,caseName) / "heat_diagnostics.csv") ;
    heat_flux_transfer_info << "west_wall_temperature,"
                            << "east_wall_temperature,"
                            << "ordinate,"
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

        for(int ordinate: ordinates){

            const string ordinateCase = "S" + to_string(ordinate);
            paths.createOutputDirectory(parentFolder, caseName, subCaseName + "/" + ordinateCase) ;

            // Retrieves the mesh files for the current mesh configuration.
            meshFiles = paths.getMeshFilenames(parentFolder, caseName, dimensions) ;

            // Defines the field files for the initial G and T data.
            fieldFiles = {paths.getTFieldFile(parentFolder, caseName, subCaseName)} ;

            // Constructs the SIMPLE solver for the current case and executes the full simulation.
            string domCaseName = parentFolder + "/" + caseName ;
            string domSubCaseName = subCaseName + "/" + ordinateCase ;

            DOM dom(meshFiles, fieldFiles, a, k, sigmaS, emissivity, dimensions, domCaseName, domSubCaseName, ordinate) ;
            auto start = std::chrono::high_resolution_clock::now();
            dom.runDOMPureRadiationNoParticipatingMedium() ;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            dom.heatData()->calculateWestEastTransferForParallelPlates() ;
            map<string,double> GAnalytical = dom.heatData()->parallelPlateHeatFluxBlackbodyAnalyticalValue(Twest,stod(subCaseName)) ;

            heat_flux_transfer_info << Twest << ","
                                    << subCaseName << ","
                                    << ordinate << ","
                                    << elapsed.count() << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["east"] << ","
                                    << GAnalytical["east"] << ","
                                    << dom.heatData()->averageWallNormalHeatFlux()["west"] << ","
                                    << GAnalytical["west"] << ","
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
                                    << dom.heatData()->oppositeWallHeatTransfer()["west"] << endl ;

            ofstream centerline_G(paths.getOutputDirectory(parentFolder, caseName, domSubCaseName) / "centerlineGInX.csv") ;
            centerline_G << "west_wall_temperature," << Twest << endl ;
            centerline_G << "east_wall_temperature," << subCaseName << endl ;
            centerline_G << "ordinateS," << ordinate << endl << endl ;
            centerline_G << "x_centerline_coordinate, G" << endl ;
            for(auto& [key,val]: dom.heatData()->centerlineGInX()){ centerline_G << key << "," << val << endl ; }
            centerline_G << endl ;
            
        }

    }

}

    void domParallelPlateSimulationPureRadiationParticipatingMedium(){
    /*
        Heat transfer between two parallel plates behaving as infinite plates.
        The west and east walls are maintained at 1000 K and 500 K,
        respectively, while the absorption coefficient is varied.
        The medium is non-scattering and non-conducting.
    */

    PathConfig paths;

    try {
        paths.validate();
        paths.ensureOutputDir();
    }
    catch (const exception& e) {
        cerr << "Path error: " << e.what() << endl;
        throw runtime_error(e.what());
    }

    const string parentFolder = "DOMPreliminaryValidation";
    const string caseName =
        "ParallelPlateParticipatingMediumClosedDOM";

    paths.createOutputDirectory(parentFolder, caseName);

    vector<string> meshFiles;
    vector<string> fieldFiles;

    const string dimensions = "300x1x1";

    /*
        Since sigmaS = 0 and L = 1 m, the numerical values of the
        absorption coefficient and optical thickness are identical.
    */
    const vector<string> subCaseNames = {
        "0.1", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0",
        "3.5", "4.0", "5.0", "7.0", "9.0", "10.0"
    };

    const double sigmaS = 0.0;
    const double k = 0.0;

    const vector<double> emissivity = {
        1.0, 1.0, 0.0, 0.0, 0.0, 0.0
    };

    const double Twest = 1000.0;
    const double Teast = 500.0;

    const vector<int> ordinates = {2, 4, 6, 8};

    const filesystem::path caseOutputDirectory =
        paths.getOutputDirectory(parentFolder, caseName);

    ofstream heat_flux_transfer_info(
        caseOutputDirectory / "heat_diagnostics.csv"
    );

    if (!heat_flux_transfer_info) {
        throw runtime_error(
            "Unable to create heat_diagnostics.csv in "
            + caseOutputDirectory.string()
        );
    }

    heat_flux_transfer_info
        << "west_wall_temperature,"
        << "east_wall_temperature,"
        << "ordinate,"
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
        << "dimensionless_heat_flux_bottom"
        << '\n';

    /*
        Approximate digitised analytical data from Fig. 4(b) of
        Siegel and Howell.
    */
    ofstream analytical_dimensionless_heat_flux(
        caseOutputDirectory
        / "analytical_dimensionless_heat_flux.csv"
    );

    if (!analytical_dimensionless_heat_flux) {
        throw runtime_error(
            "Unable to create analytical_dimensionless_heat_flux.csv"
        );
    }

    analytical_dimensionless_heat_flux
        << "absorption_coefficient,dimension_heat_flux_east\n";

    const vector<double> tauL_heatFlux = {
        0.25, 0.50, 0.75, 1.00, 1.25, 1.50,
        1.75, 2.00, 2.25, 2.50, 2.75, 3.00
    };

    const vector<double> psi_heatFlux = {
        0.81, 0.70, 0.61, 0.55, 0.50, 0.45,
        0.42, 0.39, 0.36, 0.34, 0.32, 0.30
    };

    for (size_t i = 0; i < tauL_heatFlux.size(); ++i) {
        analytical_dimensionless_heat_flux
            << tauL_heatFlux[i] << ","
            << psi_heatFlux[i] << '\n';
    }

    ofstream analytical_dimensionless_centerline_temperature(
        caseOutputDirectory
        / "analytical_dimensionless_centerline_temperature.csv"
    );

    if (!analytical_dimensionless_centerline_temperature) {
        throw runtime_error(
            "Unable to create "
            "analytical_dimensionless_centerline_temperature.csv"
        );
    }

    analytical_dimensionless_centerline_temperature
        << "x/l,0.1,1.0,2.0,10.0\n";

    const vector<double> xOverL = {
        0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0
    };

    const vector<double> phi_tau_0p1 = {
        0.58, 0.55, 0.54, 0.54, 0.49, 0.46, 0.43
    };

    const vector<double> phi_tau_1p0 = {
        0.80, 0.75, 0.69, 0.58, 0.47, 0.36, 0.24
    };

    const vector<double> phi_tau_2p0 = {
        0.88, 0.82, 0.75, 0.60, 0.42, 0.29, 0.15
    };

    const vector<double> phi_tau_10p0 = {
        0.99, 0.91, 0.81, 0.58, 0.33, 0.16, 0.06
    };

    for (size_t i = 0; i < xOverL.size(); ++i) {
        analytical_dimensionless_centerline_temperature
            << xOverL[i] << ","
            << phi_tau_0p1[i] << ","
            << phi_tau_1p0[i] << ","
            << phi_tau_2p0[i] << ","
            << phi_tau_10p0[i] << '\n';
    }

    for (const string& subCaseName : subCaseNames) {

        for (const int ordinate : ordinates) {

            const string ordinateCase =
                "S" + to_string(ordinate);

            const string domSubCaseName =
                subCaseName + "/" + ordinateCase;

            paths.createOutputDirectory(
                parentFolder,
                caseName,
                domSubCaseName
            );

            meshFiles = paths.getMeshFilenames(
                parentFolder,
                caseName,
                dimensions
            );

            /*
                The common initial fields remain in the case-level input
                folder and are shared by all absorption coefficients and
                ordinate orders.
            */
            fieldFiles = {
                paths.getTFieldFile(parentFolder, caseName),
                paths.getGFieldFile(parentFolder, caseName)
            };

            const string domCaseName =
                parentFolder + "/" + caseName;

            const double absorptionCoefficient =
                stod(subCaseName);

            DOM dom(
                meshFiles,
                fieldFiles,
                absorptionCoefficient,
                k,
                sigmaS,
                emissivity,
                dimensions,
                domCaseName,
                domSubCaseName,
                ordinate
            );

            const auto start =
                chrono::high_resolution_clock::now();

            dom.runDOMPureRadiationParticipatingMedium();

            const auto end =
                chrono::high_resolution_clock::now();

            const chrono::duration<double> elapsed =
                end - start;

            dom.heatData()
                ->calculateWestEastTransferForParallelPlates();

            const map<string, double> GAnalytical =
                dom.heatData()
                    ->parallelPlateHeatFluxBlackbodyAnalyticalValue(
                        Twest,
                        Teast
                    );

            const auto& wallFlux =
                dom.heatData()->averageWallNormalHeatFlux();

            const auto& wallTransfer =
                dom.heatData()->totalWallHeatTransfer();

            heat_flux_transfer_info
                << Twest << ","
                << Teast << ","
                << ordinate << ","
                << elapsed.count() << ","
                << absorptionCoefficient << ","
                << wallFlux.at("east") << ","
                << GAnalytical.at("east") << ","
                << wallFlux.at("west") << ","
                << GAnalytical.at("west") << ","
                << wallFlux.at("north") << ","
                << wallFlux.at("south") << ","
                << wallFlux.at("top") << ","
                << wallFlux.at("bottom") << ","
                << (
                    wallFlux.at("east")
                    + wallFlux.at("west")
                    + wallFlux.at("north")
                    + wallFlux.at("south")
                    + wallFlux.at("top")
                    + wallFlux.at("bottom")
                ) << ","
                << (
                    wallFlux.at("east")
                    + wallFlux.at("west")
                ) << ","
                << dom.heatData()
                       ->oppositeWallHeatFlux()
                       .at("east") << ","
                << dom.heatData()
                       ->oppositeWallHeatFlux()
                       .at("west") << ","
                << wallTransfer.at("east") << ","
                << wallTransfer.at("west") << ","
                << wallTransfer.at("north") << ","
                << wallTransfer.at("south") << ","
                << wallTransfer.at("top") << ","
                << wallTransfer.at("bottom") << ","
                << (
                    wallTransfer.at("east")
                    + wallTransfer.at("west")
                    + wallTransfer.at("north")
                    + wallTransfer.at("south")
                    + wallTransfer.at("top")
                    + wallTransfer.at("bottom")
                ) << ","
                << (
                    wallTransfer.at("east")
                    + wallTransfer.at("west")
                ) << ","
                << dom.heatData()
                       ->oppositeWallHeatTransfer()
                       .at("east") << ","
                << dom.heatData()
                       ->oppositeWallHeatTransfer()
                       .at("west") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("east") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("west") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("north") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("south") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("top") << ","
                << dom.heatData()
                       ->dimensionlessQw()
                       .at("bottom")
                << '\n';

            const filesystem::path runOutputDirectory =
                paths.getOutputDirectory(
                    parentFolder,
                    caseName,
                    domSubCaseName
                );

            ofstream centerline_G(
                runOutputDirectory / "centerlineGInX.csv"
            );

            ofstream centerline_T(
                runOutputDirectory / "centerlineT.csv"
            );

            ofstream dimensionless_centerline_T(
                runOutputDirectory
                / "dimensionlessCenterlineT.csv"
            );

            if (
                !centerline_G
                || !centerline_T
                || !dimensionless_centerline_T
            ) {
                throw runtime_error(
                    "Unable to create centreline output files in "
                    + runOutputDirectory.string()
                );
            }

            centerline_G
                << "west_wall_temperature," << Twest << '\n'
                << "east_wall_temperature," << Teast << '\n'
                << "absorption_coefficient,"
                << absorptionCoefficient << '\n'
                << "ordinateS," << ordinate << "\n\n"
                << "x_over_L,G\n";

            centerline_T
                << "west_wall_temperature," << Twest << '\n'
                << "east_wall_temperature," << Teast << '\n'
                << "absorption_coefficient,"
                << absorptionCoefficient << '\n'
                << "ordinateS," << ordinate << "\n\n"
                << "x_over_L,T\n";

            dimensionless_centerline_T
                << "west_wall_temperature," << Twest << '\n'
                << "east_wall_temperature," << Teast << '\n'
                << "absorption_coefficient,"
                << absorptionCoefficient << '\n'
                << "ordinateS," << ordinate << "\n\n"
                << "x_over_L,dimensionless_T\n";

            for (
                const auto& [coordinate, value]
                : dom.heatData()->centerlineGInX()
            ) {
                centerline_G
                    << coordinate << "," << value << '\n';
            }

            for (
                const auto& [coordinate, value]
                : dom.heatData()->centerlineTemperatureInX()
            ) {
                centerline_T
                    << coordinate << "," << value << '\n';
            }

            for (
                const auto& [coordinate, value]
                : dom.heatData()
                      ->dimensionlessCenterlineTemperatureInX()
            ) {
                dimensionless_centerline_T
                    << coordinate << "," << value << '\n';
            }
        }
    }
}
    
    void domSquarePlateSimulationRadiationConduction()
    {
        /*
            Coupled radiation--conduction in a square enclosure.
            The bottom wall is maintained at 1000 K, while the top,
            west and east walls are maintained at 500 K. The
            conductivity is varied through the prescribed Stark numbers.
        */

        PathConfig paths;

        try {
            paths.validate();
            paths.ensureOutputDir();
        }
        catch (const exception& e) {
            cerr << "Path error: " << e.what() << endl;
            throw runtime_error(e.what());
        }

        const string parentFolder = "DOMPreliminaryValidation";
        const string caseName =
            "SquarePlateRadiationConductionClosedDOM";

        paths.createOutputDirectory(parentFolder, caseName);

        vector<string> meshFiles;
        vector<string> fieldFiles;

        const string dimensions = "40x5x40";

        const double a = 1.0;
        const double sigmaS = 0.0;
        const double sigma = 5.670374419e-8;

        const vector<double> emissivity = {
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        };

        const double Tbottom = 1000.0;
        const double Ttop = Tbottom / 2.0;
        const double Twest = Tbottom / 2.0;
        const double Teast = Tbottom / 2.0;

        /*
            Keep the Stark-number ordering explicit so that each
            conductivity remains permanently associated with the
            correct benchmark curve.
        */
        const vector<double> starkNumbers = {
            1.0, 0.1, 0.05
        };

        const double sc =
            4.0
            * sigma
            * (
                powl(Tbottom, 3)
                - powl(Ttop, 3)
            )
            / a;

        vector<string> subCaseNames;
        subCaseNames.reserve(starkNumbers.size());

        for (const double starkNumber : starkNumbers) {
            subCaseNames.push_back(
                to_string(sc * starkNumber)
            );
        }

        const vector<int> ordinates = {2, 4, 6, 8};

        const filesystem::path caseOutputDirectory =
            paths.getOutputDirectory(parentFolder, caseName);

        ofstream heat_flux_transfer_info(
            caseOutputDirectory / "heat_diagnostics.csv"
        );

        if (!heat_flux_transfer_info) {
            throw runtime_error(
                "Unable to create heat_diagnostics.csv in "
                + caseOutputDirectory.string()
            );
        }

        heat_flux_transfer_info
            << "bottom_wall_temperature,"
            << "top_wall_temperature,"
            << "west_wall_temperature,"
            << "east_wall_temperature,"
            << "ordinate,"
            << "runtime,"
            << "absorption_coefficient,"
            << "stark_number,"
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
            << "dimensionless_heat_flux_bottom"
            << '\n';

        ofstream analytical_dimensionless_centerline_temperature(
            caseOutputDirectory
            / "analytical_dimensionless_centerline_temperature.csv"
        );

        if (!analytical_dimensionless_centerline_temperature) {
            throw runtime_error(
                "Unable to create "
                "analytical_dimensionless_centerline_temperature.csv"
            );
        }

        analytical_dimensionless_centerline_temperature
            << "y/l,N0.05,N0.1,N1\n";

        const vector<double> yOverLy_fig7 = {
            0.00, 0.05, 0.10, 0.20, 0.30, 0.40,
            0.50, 0.60, 0.70, 0.80, 0.90, 1.00
        };

        const vector<double> TbyTw_N_1_fig7 = {
            1.000, 0.951, 0.901, 0.815, 0.737, 0.677,
            0.629, 0.592, 0.562, 0.539, 0.519, 0.500
        };

        const vector<double> TbyTw_N_0p1_fig7 = {
            1.000, 0.952, 0.905, 0.822, 0.757, 0.709,
            0.664, 0.624, 0.593, 0.561, 0.534, 0.500
        };

        const vector<double> TbyTw_N_0p05_fig7 = {
            1.000, 0.953, 0.907, 0.828, 0.768, 0.726,
            0.685, 0.649, 0.613, 0.578, 0.545, 0.500
        };

        for (size_t i = 0; i < yOverLy_fig7.size(); ++i) {
            analytical_dimensionless_centerline_temperature
                << yOverLy_fig7[i] << ","
                << TbyTw_N_0p05_fig7[i] << ","
                << TbyTw_N_0p1_fig7[i] << ","
                << TbyTw_N_1_fig7[i] << '\n';
        }

        for (
            size_t caseIndex = 0;
            caseIndex < subCaseNames.size();
            ++caseIndex
        ) {
            const string& subCaseName =
                subCaseNames[caseIndex];

            const double conductivity =
                stod(subCaseName);

            const double starkNumber =
                starkNumbers[caseIndex];

            for (const int ordinate : ordinates) {

                const string ordinateCase =
                    "S" + to_string(ordinate);

                const string domSubCaseName =
                    subCaseName + "/" + ordinateCase;

                paths.createOutputDirectory(
                    parentFolder,
                    caseName,
                    domSubCaseName
                );

                meshFiles = paths.getMeshFilenames(
                    parentFolder,
                    caseName,
                    dimensions
                );

                fieldFiles = {
                    paths.getTFieldFile(parentFolder, caseName),
                    paths.getGFieldFile(parentFolder, caseName)
                };

                const string domCaseName =
                    parentFolder + "/" + caseName;

                DOM dom(
                    meshFiles,
                    fieldFiles,
                    a,
                    conductivity,
                    sigmaS,
                    emissivity,
                    dimensions,
                    domCaseName,
                    domSubCaseName,
                    ordinate
                );

                const auto start =
                    chrono::high_resolution_clock::now();

                dom.runDOMRadiationConduction();

                const auto end =
                    chrono::high_resolution_clock::now();

                const chrono::duration<double> elapsed =
                    end - start;

                dom.heatData()
                    ->calculateWestEastTransferForParallelPlates();

                const auto& wallFlux =
                    dom.heatData()->averageWallNormalHeatFlux();

                const auto& wallTransfer =
                    dom.heatData()->totalWallHeatTransfer();

                heat_flux_transfer_info
                    << Tbottom << ","
                    << Ttop << ","
                    << Twest << ","
                    << Teast << ","
                    << ordinate << ","
                    << elapsed.count() << ","
                    << a << ","
                    << starkNumber << ","
                    << conductivity << ","
                    << wallFlux.at("east") << ","
                    << wallFlux.at("west") << ","
                    << wallFlux.at("north") << ","
                    << wallFlux.at("south") << ","
                    << wallFlux.at("top") << ","
                    << wallFlux.at("bottom") << ","
                    << (
                        wallFlux.at("east")
                        + wallFlux.at("west")
                        + wallFlux.at("north")
                        + wallFlux.at("south")
                        + wallFlux.at("top")
                        + wallFlux.at("bottom")
                    ) << ","
                    << (
                        wallFlux.at("east")
                        + wallFlux.at("west")
                    ) << ","
                    << dom.heatData()
                        ->oppositeWallHeatFlux()
                        .at("east") << ","
                    << dom.heatData()
                        ->oppositeWallHeatFlux()
                        .at("west") << ","
                    << wallTransfer.at("east") << ","
                    << wallTransfer.at("west") << ","
                    << wallTransfer.at("north") << ","
                    << wallTransfer.at("south") << ","
                    << wallTransfer.at("top") << ","
                    << wallTransfer.at("bottom") << ","
                    << (
                        wallTransfer.at("east")
                        + wallTransfer.at("west")
                        + wallTransfer.at("north")
                        + wallTransfer.at("south")
                        + wallTransfer.at("top")
                        + wallTransfer.at("bottom")
                    ) << ","
                    << (
                        wallTransfer.at("east")
                        + wallTransfer.at("west")
                    ) << ","
                    << dom.heatData()
                        ->oppositeWallHeatTransfer()
                        .at("east") << ","
                    << dom.heatData()
                        ->oppositeWallHeatTransfer()
                        .at("west") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("east") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("west") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("north") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("south") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("top") << ","
                    << dom.heatData()
                        ->dimensionlessQw()
                        .at("bottom")
                    << '\n';

                const filesystem::path runOutputDirectory =
                    paths.getOutputDirectory(
                        parentFolder,
                        caseName,
                        domSubCaseName
                    );

                ofstream centerline_G(
                    runOutputDirectory / "centerlineGInX.csv"
                );

                ofstream centerline_T(
                    runOutputDirectory / "centerlineTInX.csv"
                );

                ofstream dimensionless_centerline_T(
                    runOutputDirectory
                    / "dimensionlessCenterlineTInX.csv"
                );

                ofstream centerline_G_In_Z(
                    runOutputDirectory / "centerlineGInZ.csv"
                );

                ofstream centerline_T_In_Z(
                    runOutputDirectory / "centerlineTInZ.csv"
                );

                ofstream dimensionless_centerline_T_In_Z(
                    runOutputDirectory
                    / "dimensionlessCenterlineTInZ.csv"
                );

                if (
                    !centerline_G
                    || !centerline_T
                    || !dimensionless_centerline_T
                    || !centerline_G_In_Z
                    || !centerline_T_In_Z
                    || !dimensionless_centerline_T_In_Z
                ) {
                    throw runtime_error(
                        "Unable to create centreline output files in "
                        + runOutputDirectory.string()
                    );
                }

                const auto writeMetadata =
                    [&](ofstream& output) {
                        output
                            << "west_wall_temperature,"
                            << Twest << '\n'
                            << "east_wall_temperature,"
                            << Teast << '\n'
                            << "top_wall_temperature,"
                            << Ttop << '\n'
                            << "bottom_wall_temperature,"
                            << Tbottom << '\n'
                            << "absorption_coefficient,"
                            << a << '\n'
                            << "stark_number,"
                            << starkNumber << '\n'
                            << "conductivity,"
                            << conductivity << '\n'
                            << "ordinateS,"
                            << ordinate << "\n\n";
                    };

                writeMetadata(centerline_G);
                writeMetadata(centerline_T);
                writeMetadata(dimensionless_centerline_T);
                writeMetadata(centerline_G_In_Z);
                writeMetadata(centerline_T_In_Z);
                writeMetadata(dimensionless_centerline_T_In_Z);

                centerline_G
                    << "x_over_L,G\n";

                centerline_T
                    << "x_over_L,T\n";

                dimensionless_centerline_T
                    << "x_over_L,dimensionless_T\n";

                centerline_G_In_Z
                    << "z_over_L,G\n";

                centerline_T_In_Z
                    << "z_over_L,T\n";

                dimensionless_centerline_T_In_Z
                    << "z_over_L,dimensionless_T\n";

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()->centerlineGInX()
                ) {
                    centerline_G
                        << coordinate << "," << value << '\n';
                }

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()->centerlineTemperatureInX()
                ) {
                    centerline_T
                        << coordinate << "," << value << '\n';
                }

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()
                        ->dimensionlessCenterlineTemperatureInX()
                ) {
                    dimensionless_centerline_T
                        << coordinate << "," << value << '\n';
                }

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()->centerlineGInZ()
                ) {
                    centerline_G_In_Z
                        << coordinate << "," << value << '\n';
                }

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()->centerlineTemperatureInZ()
                ) {
                    centerline_T_In_Z
                        << coordinate << "," << value << '\n';
                }

                for (
                    const auto& [coordinate, value]
                    : dom.heatData()
                        ->dimensionlessCenterlineTemperatureInZ()
                ) {
                    dimensionless_centerline_T_In_Z
                        << coordinate << "," << value << '\n';
                }
            }
        }
    }

} ;