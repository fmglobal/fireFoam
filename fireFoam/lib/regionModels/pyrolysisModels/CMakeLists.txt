set(lib additionalPyrolysisModels)

set(SOURCES
    reactingOneDimFF/reactingOneDimFF.C
    )

include_directories(lnInclude)

add_custom_command(
    OUTPUT lnInclude/uptodate
    DEPENDS CMakeLists.txt
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND rm -rf lnInclude
    COMMAND wmakeLnInclude .
    COMMAND touch lnInclude/uptodate
    )

add_custom_target(${lib}_lnInclude DEPENDS lnInclude/uptodate)
add_library(${lib} SHARED ${SOURCES})
add_dependencies(${lib} ${lib}_lnInclude)
add_dependencies(${lib} surfaceFilmModels)

target_include_directories(${lib} PUBLIC
    ${FF_PKG_DIR}/regionModels/surfaceFilmModels/lnInclude

    ${FF_LIB_DIR}/regionCoupling/lnInclude
    ${FF_LIB_DIR}/fvPatchFieldsPyrolysis/lnInclude
    ${FF_LIB_DIR}/fvPatchFields/lnInclude

    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermoTools/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidChemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/radiation/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/OpenFOAM/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_finiteVolume}
    ${OF_OpenFOAM}
    ${OF_meshTools}
    ${OF_chemistryModel}
    ${OF_specie}
    ${OF_fluidThermophysicalModels}
    ${OF_solidChemistryModel}
    ${OF_solidThermo}
    ${OF_turbulenceModels}
    ${OF_compressibleTurbulenceModels}
    ${OF_regionModels}
    ${OF_reactionThermophysicalModels}
    ${OF_fvPatchFieldsPyrolysis}
    ${OF_radiation}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libsurfaceFilmModels.so
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
