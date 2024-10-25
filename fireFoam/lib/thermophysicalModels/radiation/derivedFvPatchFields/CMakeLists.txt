set(lib RadiationBCs_FMG)

set(SOURCES
    radiationCoupledBaseFF/radiationCoupledBaseFF.C
    greyDiffusiveRadiationFF/greyDiffusiveRadiationFFMixedFvPatchScalarField.C
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

    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/radiation/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/surfMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/parallel/distributed/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_finiteVolume}
    ${OF_meshTools}
    ${OF_regionModels}
    ${OF_pyrolysisModels}
    ${OF_radiation}
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
)
