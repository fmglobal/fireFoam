add_executable(${PROJECT_NAME} fireFoam.C)
add_dependencies(${PROJECT_NAME} patchFieldsFireFOAM)
add_dependencies(${PROJECT_NAME} combustionModelsEDC)
add_dependencies(${PROJECT_NAME} greyMeanAbsorptionEmission_FMG)
add_dependencies(${PROJECT_NAME} fvPatchFieldsPyrolysis)
add_dependencies(${PROJECT_NAME} surfaceFilmModels)
add_dependencies(${PROJECT_NAME} additionalPyrolysisModels)
add_dependencies(${PROJECT_NAME} regionCouplingFF)
#add_dependencies(${PROJECT_NAME} lagrangianIntermediate)
add_dependencies(${PROJECT_NAME} RadiationBCs_FMG)

target_include_directories(${PROJECT_NAME} PUBLIC

    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include

    ${FF_PKG_DIR}/regionModels/surfaceFilmModels/lnInclude
    #${FF_PKG_DIR}/lagrangian/intermediateNew/lnInclude


    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/regionFaModels/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteArea/lnInclude
    ${OpenFOAM_SRC_DIR}/faOptions/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
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
    ${OpenFOAM_SRC_DIR}/combustionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/intermediate/lnInclude
    ${OpenFOAM_SRC_DIR}/ODE/lnInclude
    ${OpenFOAM_SRC_DIR}/surfMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/parallel/distributed/lnInclude
    ${OpenFOAM_SRC_DIR}/dynamicFvMesh/lnInclude
    )

target_link_libraries(${PROJECT_NAME} PUBLIC

    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libpatchFieldsFireFOAM.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libcombustionModelsEDC.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libgreyMeanAbsorptionEmission_FMG.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libfvPatchFieldsPyrolysis.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libsurfaceFilmModels.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libadditionalPyrolysisModels.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libregionCouplingFF.so
    #${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/liblagrangianIntermediate.so
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRadiationBCs_FMG.so

    ${CMAKE_DL_LIBS}

    ${OF_OpenFOAM}
    ${OF_Pstream}
    ${OF_finiteVolume}
    ${OF_fvOptions}
    ${OF_meshTools}
    ${OF_sampling}
    ${OF_specie}
    ${OF_dynamicFvMesh}
    ${OF_fluidThermophysicalModels}
    ${OF_reactionThermophysicalModels}
    ${OF_SLGThermo}
    ${OF_chemistryModel}
    ${OF_solidChemistryModel}
    ${OF_combustionModels}
    ${OF_turbulenceModels}
    ${OF_compressibleTurbulenceModels}
    ${OF_thermophysicalProperties}
    ${OF_lagrangian}
    ${OF_ODE}
    ${OF_regionModels}
    ${OF_fieldFunctionObjects}
    ${OF_lagrangianFunctionObjects}
    ${OF_pyrolysisModels}
    ${OF_radiation}
    )

install(TARGETS             ${PROJECT_NAME}
    RUNTIME DESTINATION $ENV{FOAM_USER_APPBIN}
)
