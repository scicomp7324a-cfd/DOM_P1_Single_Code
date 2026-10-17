# DiscreteOrdinateModel_P1Approximation_Radiation_Model
This repository contains the source code of Discrete Ordinates Model and P1 Approximation submitted as part of MPhil dissertation at the University of Cambridge.

P1 and Discrete Ordinates Validation Code
==========================================

This folder contains the finite-volume P1 and Discrete Ordinates Model implementations used for the validation and model-comparison part of the thesis. The same mesh, sparse-matrix and finite-volume framework is shared by both radiation models.

Requirements
------------

The C++ code needs:

- CMake 3.16 or newer
- a C++17 compiler
- a Unix-like terminal environment

The plotting scripts use Python 3. The main Python packages are:

    numpy
    pandas
    matplotlib

They can be installed, if required, with:

    pip install numpy pandas matplotlib

Code layout
-----------

The main files and folders are:

    main.cpp
        Selects and runs the validation tests.

    P1.hpp
        P1 radiation model.

    DOM.hpp
        Discrete Ordinates Model.

    PathConfig.hpp
        Input and output path handling.

    Mesh/
        Mesh reading and face addressing.

    FVCore/
        Fields, boundary conditions, ordinate data and other finite-volume utilities.

    Discretisation/
        Diffusion, convection and source-term discretisation.

    SparseOperations/
        Sparse matrix, vector and addressing classes.

    Solver/
        Iterative linear solvers.

    ValidationTests/
        Individual validation cases and their plotting scripts.

    Files/
        Mesh files and initial G and T field files used by the cases.

    Output/
        Simulation results.

Building and running
--------------------

From this folder run:

    chmod +x run.sh
    ./run.sh

The script performs a clean Release build and runs the executable p1_app from inside the build directory.

The same steps can be done manually with:

    rm -rf build
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    cd build
    ./p1_app

CMake copies Files/ into build/Files/ during the build. The executable should therefore normally be run from the build directory. Results are written to ../Output/, which is the Output folder in the project root.

Tests
-----

The supplied main.cpp runs the validation groups one after another. These include:

- cold non-scattering parallel plates
- cold isotropically scattering parallel plates
- finite open parallel-plate model comparison
- hot-wall/cold-medium rectangular enclosure
- annular-source rectangular enclosure
- isotropic and anisotropic scattering
- preliminary P1 validation
- preliminary DOM validation

To run fewer cases, comment out the unwanted .test() calls in main.cpp and rebuild.

Plotting
--------

Each main validation group has its plotting script inside ValidationTests/. Run the C++ cases first, then run the required script from the project root.

P1 preliminary validation:

    python3 ValidationTests/P1PreliminaryValidation/plot_P1_preliminary_validation.py --output-root Output

DOM preliminary validation:

    python3 ValidationTests/DOMPreliminaryValidation/plot_DOM_all_ordinates.py --output-root Output

Cold isotropic scattering:

    python3 ValidationTests/DOMP1ComparisionColdScatteringParallelPlatePhysicalSeperation/plot_cold_isotropic_scattering.py Output

Isotropic and anisotropic scattering:

    python3 ValidationTests/DOMP1ComparisonAnisotropicAndIsotropicScattering/plot_anisotropic_scattering.py Output

Hot-wall/cold-medium enclosure:

    python3 ValidationTests/DOMP1ComparisionRectangularEnclosureHotWallColdMedium/plot_case3_hot_wall_cold_medium.py --output-root Output

Annular-source enclosure:

    python3 ValidationTests/DOMP1ComparisionRectangularEnclosureAnnularSource/plot_case4_annular_source.py --output-root Output

The finite open-plate memory/accuracy plotting script compares the P1/DOM CSV with the DTM CSV. Place both unified CSV files in the same directory and run, for example:

    python3 ValidationTests/DOMP1ComparisionFiniteParallelPlateMediumOpen/plot_memory_accuracy_comparison.py --input-directory Output/DOMP1FiniteParallelPlateMemoryAccuracy

The scripts save their figures inside the relevant case output folder unless another figure directory is supplied.

Files and paths
---------------

The Files/ directory contains the mesh and field data required by the validation cases and should not be moved relative to the project root. PathConfig.hpp assumes the executable is run from build/, where CMake has copied this folder.

VTK files written by some cases can be opened separately in ParaView.

