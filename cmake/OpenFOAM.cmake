# Check valid OpenFOAM
if(DEFINED ENV{WM_PROJECT_DIR})
    message(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR "The OpenFOAM bashrc is not sourced")
endif(DEFINED ENV{WM_PROJECT_DIR})

set(OpenFOAM_VERSION $ENV{WM_PROJECT_VERSION})
set(OpenFOAM_DIR $ENV{WM_PROJECT_DIR})
set(OpenFOAM_LIB_DIR $ENV{FOAM_LIBBIN})
set(OpenFOAM_SRC_DIR $ENV{FOAM_SRC})

message(STATUS "OpenFOAM lib path: ${OpenFOAM_LIB_DIR}")
message(STATUS "OpenFOAM src path: ${OpenFOAM_SRC_DIR}")

# Read build relevant OpenFOAM environment variables and cache them
set(of_wm_project_dir "$ENV{WM_PROJECT_DIR}" CACHE PATH "Path to OpenFOAM project folder.")
set(of_wm_arch "$ENV{WM_ARCH}" CACHE STRING "Architecture. Usually linux64.")
set(of_wm_arch_option "$ENV{WM_ARCH_OPTION}" CACHE STRING "Information if 32 or 64 bit operating system.")
set(of_wm_precision_option "$ENV{WM_PRECISION_OPTION}" CACHE STRING "Flag if to use single precision (SP) or double precision (DP).")
set(of_wm_label_size "$ENV{WM_LABEL_SIZE}" CACHE STRING "Size in bit to use as label type. Can be either 32 or 64.")
set(of_wm_compile_option "$ENV{WM_COMPILE_OPTION}" CACHE STRING "OpenFOAM build type: Opt, Debug, Prof.")
set(of_wm_compiler "$ENV{WM_COMPILER}" CACHE STRING "Compiler used for OpenFOAM build.")
set(of_wm_label_option "$ENV{WM_LABEL_OPTION}" CACHE STRING "Concrete Type used for label. Either Int32 or Int64.")
# WM_ARCH + WM_COMPILER + WM_PRECISION_OPTION + WM_LABEL_OPTION + WM_COMPILE_OPTION
set(of_wm_options "${of_wm_arch}${of_wm_compiler}${of_wm_precision_option}${of_wm_label_option}" CACHE STRING "Name of subfolder which contains compiled OpenFOAM libraries" FORCE)

# Compile definitions required for OpenFOAM
add_compile_definitions(
    WM_LABEL_SIZE=${of_wm_label_size}
    WM_${of_wm_precision_option}
    WM_ARCH_OPTION=${of_wm_arch_option}
    ${of_wm_arch}
    OPENFOAM="$ENV{WM_PROJECT_VERSION}" # Figures out OF version. TM. 
    NoRepository
)

# Required to make linking to OpenFOAM libraries work
set(CMAKE_EXE_LINKER_FLAGS "-Xlinker --add-needed -Xlinker --no-as-needed")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -DDEBUG -D_USE_MATH_DEFINES -DBOOST_MATH_INSTRUMENT")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -Wno-unused-variable -Wno-deprecated -D_USE_MATH_DEFINES")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -ggdb3 -O3 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -DDEBUG -D_USE_MATH_DEFINES")

# =====================================================================================================
# Of course you can change some parameters here, e.g., some pre-settings in ~/.OpenFOAM/prefs.h
set(PATH_LIB_OPENMPI "sys-openmpi")
# set(DEFINITIONS_COMPILE "-std=c++11 -m64 -pthread -ftrapping-math -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -Wno-unknown-pragmas -O3 -DNoRepository -ftemplate-depth-100 -fPIC")

# Compiling configure
# add_definitions("${DEFINITIONS_COMPILE}")
# # ======== OS specific setting =============
# if(APPLE)
#     add_definitions(" -Ddarwin64 ")
# else()
#     add_definitions("-Dlinux64")
# endif(APPLE)

# ==========================================

# # Required OpenFOAM libraries
find_library(OF_OpenFOAM OpenFOAM PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_finiteVolume finiteVolume PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_meshTools meshTools PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_sampling sampling PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_turbulenceModels turbulenceModels  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_compressibleTurbulenceModels compressibleTurbulenceModels  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_lagrangian lagrangian  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_lagrangianFunctionObjects lagrangianFunctionObjects  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_specie specie PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_solidThermo solidThermo PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_solidChemistryModel solidChemistryModel  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_compressibleTransportModels compressibleTransportModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_thermophysicalProperties thermophysicalProperties PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_reactionThermophysicalModels reactionThermophysicalModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_SLGThermo SLGThermo PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_Pstream Pstream PATHS ${OpenFOAM_LIB_DIR}/${PATH_LIB_OPENMPI} REQUIRED)
find_library(OF_distributed distributed PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_fvOptions fvOptions  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_surfMesh surfMesh PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_regionModels regionModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_fluidThermophysicalModels fluidThermophysicalModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_chemistryModel chemistryModel PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_dynamicMesh dynamicMesh  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_dynamicFvMesh dynamicFvMesh  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_combustionModels combustionModels  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_ODE ODE  PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_fieldFunctionObjects fieldFunctionObjects PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_topoChangerFvMesh topoChangerFvMesh PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_decompositionMethods decompositionMethods PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_pyrolysisModels pyrolysisModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)
find_library(OF_radiationModels radiationModels PATHS ${OpenFOAM_LIB_DIR} REQUIRED)

find_library(OF_kahipDecomp kahipDecomp PATHS ${OpenFOAM_LIB_DIR}/dummy REQUIRED)
find_library(OF_metisDecomp metisDecomp PATHS ${OpenFOAM_LIB_DIR}/dummy REQUIRED)
find_library(OF_ptscotchDecomp ptscotchDecomp PATHS ${OpenFOAM_LIB_DIR}/dummy REQUIRED)
find_library(OF_scotchDecomp scotchDecomp PATHS ${OpenFOAM_LIB_DIR}/dummy REQUIRED)


include_directories(
    ${OpenFOAM_SRC_DIR}/OpenFOAM/lnInclude
    ${OpenFOAM_SRC_DIR}/OSspecific/POSIX/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/sampling/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/distributionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solid/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidChemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/combustionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalFunctions/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/radiation/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/ODE/lnInclude
    ${OpenFOAM_SRC_DIR}/surfMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/parallel/distributed/lnInclude
    ${OpenFOAM_SRC_DIR}/dynamicMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/dynamicFvMesh/lnInclude
    )

