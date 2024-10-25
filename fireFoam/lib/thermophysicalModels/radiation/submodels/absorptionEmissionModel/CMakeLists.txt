set(lib greyMeanAbsorptionEmission_FMG)

set(SOURCES
    constRadFractionEmission/constRadFractionEmission.C
    RVFGasSootEmission/RVFGasSootEmission.C
    greyMeanSolidAbsorptionEmissionFF/greyMeanSolidAbsorptionEmissionFF.C
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
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude

    ${FF_LIB_DIR}/thermophysicalModels/combustionModels/lnInclude

    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/surfMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/parallel/distributed/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/combustionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/radiation/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/LES/LESdeltas/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/LES/LESfilters/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/turbulenceModel
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/LES/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_finiteVolume}
    ${OF_combustionModels}
    ${OF_distributed}
    ${OF_reactionThermophysicalModels}
    ${OF_radiation}
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
