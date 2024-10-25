set(lib patchFieldsFireFOAM)

set(SOURCES
    constHTemperatureFvPatchScalarField/constHTemperatureFvPatchScalarField.C
    inertPlateOneSideHeating/inertPlateOneSideHeatingFvPatchScalarField.C
    alphatPyrolysisMeshCorrected/alphatPyrolysisMeshCorrectedFvPatchScalarField.C
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

target_include_directories(${lib} PUBLIC
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/turbulenceModel
    ${OpenFOAM_SRC_DIR}/turbulenceModels/compressible/RAS/lnInclude
    ${OpenFOAM_SRC_DIR}/triSurface/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/LES/LESdeltas/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels/LES/LESfilters/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidThermo/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_finiteVolume}
    ${OF_meshTools}
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
