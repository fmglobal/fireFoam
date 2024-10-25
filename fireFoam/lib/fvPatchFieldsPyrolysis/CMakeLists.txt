set(lib fvPatchFieldsPyrolysis)

set(SOURCES
    pyrolysisTemperatureFixedRadiation/pyrolysisTemperatureFixedRadiationFvPatchScalarField.C
    pyrolysisVelocityCoupled/pyrolysisVelocityCoupledFvPatchVectorField.C
    pyrolysisTemperatureCoupled/pyrolysisTemperatureCoupledFvPatchScalarField.C
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
add_dependencies(${lib} RadiationBCs_FMG surfaceFilmModels)

target_include_directories(${lib} PUBLIC

    ${FF_PKG_DIR}/regionModels/surfaceFilmModels/lnInclude
    ${FF_PKG_DIR}/thermophysicalModels/radiation/lnInclude

    ${FF_LIB_DIR}/thermophysicalModels/radiation/derivedFvPatchFields/lnInclude

    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermoTools/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidSpecie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solid/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basicSolidThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidChemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/turbulenceModel
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/turbulenceModel/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/RAS/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/LES/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/LES/LESdeltas/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels//LES/LESfilters/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/sampling/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_specie}
    ${OF_finiteVolume}
    ${OF_regionModels}
    ${OF_pyrolysisModels}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRadiationBCs_FMG.so
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
