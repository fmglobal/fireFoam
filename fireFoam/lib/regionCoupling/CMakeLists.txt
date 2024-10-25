set(lib regionCouplingFF)

set(SOURCES
    derivedFvPatchFields/filmPyrolysisTemperatureCoupledFF/filmPyrolysisTemperatureCoupledFFFvPatchScalarField.C
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
    ${FF_PKG_DIR}/thermophysicalModels/radiation/lnInclude
    ${FF_PKG_DIR}/regionModels/surfaceFilmModels/lnInclude

    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermoTools/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/thermoTools/lnInclude
    ${OpenFOAM_SRC_DIR}/sampling/lnInclude
    ${OpenFOAM_SRC_DIR}/AMIInterpolation/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solid/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidChemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/properties/liquidProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/properties/liquidProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/properties/liquidMixtureProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/properties/solidProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/properties/solidMixtureProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/turbulenceModel/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/RAS/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/LES/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/LES/LESdeltas/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_regionModels}
    ${OF_solidChemistryModel}
    ${OF_finiteVolume}
    ${OF_meshTools}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libsurfaceFilmModels.so
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
