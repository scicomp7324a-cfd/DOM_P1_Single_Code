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

# include "DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation.hpp"
# include "DOMP1ComparisionColdScatteringParallelPlatePhysicalSeperation.hpp"
# include "DOMP1ComparisionRectangularEnclosureAnnularSource.hpp"
# include "DOMP1ComparisionRectangularEnclosureHotWallColdMedium.hpp"
# include "DOMP1ComparisonAnisotropicAndIsotropicScattering.hpp"
# include "DOMP1ComparisionFiniteParallelPlateMediumOpen.hpp"
# include "DOMPreliminaryValidation.hpp"
# include "P1PreliminaryValidation.hpp"

namespace fs = std::filesystem;
using namespace std ;

int main(){

    DOMP1ComparisionColdNonScatteringParallelPlatePhysicalSeperation DOMP1CCNSPPPS ;
    DOMP1CCNSPPPS.test() ;

    DOMP1ComparisionColdScatteringParallelPlatePhysicalSeperation DOMP1CCSPPPS ;
    DOMP1CCSPPPS.test() ;

    DOMP1ComparisionFiniteParallelPlateMediumOpen DOMP1CFPPMO ;
    DOMP1CFPPMO.test() ;

    DOMP1ComparisionRectangularEnclosureHotWallColdMedium DOMP1CREHWCM ;
    DOMP1CREHWCM.test() ;

    DOMP1ComparisionRectangularEnclosureAnnularSource DOMP1CREAS ;
    DOMP1CREAS.test() ;

    DOMP1ComparisonAnisotropicAndIsotropicScattering DOMP1CAAIS ;
    DOMP1CAAIS.test() ;

    P1PreliminaryValidation P1PV ;
    P1PV.test() ;

    DOMPreliminaryValidation DOMPV ;
    DOMPV.test() ;

}







