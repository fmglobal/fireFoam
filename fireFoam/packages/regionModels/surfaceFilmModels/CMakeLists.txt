set(lib surfaceFilmModels)

list(APPEND SOURCES
    surfaceFilmModel/surfaceFilmModel.C
    surfaceFilmModel/surfaceFilmModelNew.C
    surfaceFilmRegionModel/surfaceFilmRegionModel.C
    noFilm/noFilm.C
    kinematicSingleLayer/kinematicSingleLayer.C
    thermoSingleLayer/thermoSingleLayer.C
    )

list(APPEND SOURCES
    submodels/filmSubModelBase.C
    )

list(APPEND SOURCES
    submodels/kinematic/force/force/force.C
    submodels/kinematic/force/force/forceNew.C
    submodels/kinematic/force/forceList/forceList.C
    submodels/kinematic/force/contactAngleForces/contactAngleForce/contactAngleForce.C
    submodels/kinematic/force/contactAngleForces/distribution/distributionContactAngleForce.C
    submodels/kinematic/force/contactAngleForces/temperatureDependent/temperatureDependentContactAngleForce.C
    submodels/kinematic/force/contactAngleForces/perturbedTemperatureDependent/perturbedTemperatureDependentContactAngleForce.C
    submodels/kinematic/force/thermocapillaryForce/thermocapillaryForce.C
    )

list(APPEND SOURCES
    submodels/kinematic/injectionModel/injectionModel/injectionModel.C
    submodels/kinematic/injectionModel/injectionModel/injectionModelNew.C
    submodels/kinematic/injectionModel/injectionModelList/injectionModelList.C
    submodels/kinematic/injectionModel/drippingInjection/drippingInjection.C
    submodels/kinematic/injectionModel/BrunDrippingInjection/BrunDrippingInjection.C
    submodels/kinematic/injectionModel/curvatureSeparation/curvatureSeparation.C
    )

list(APPEND SOURCES
    submodels/kinematic/transferModels/transferModel/transferModel.C
    submodels/kinematic/transferModels/transferModel/transferModelNew.C
    submodels/kinematic/transferModels/transferModelList/transferModelList.C
    )

list(APPEND SOURCES
    submodels/kinematic/filmThermoModel/filmThermoModel/filmThermoModel.C
    submodels/kinematic/filmThermoModel/filmThermoModel/filmThermoModelNew.C
    submodels/kinematic/filmThermoModel/constantFilmThermo/constantFilmThermo.C
    submodels/kinematic/filmThermoModel/liquidFilmThermo/liquidFilmThermo.C
    )

list(APPEND SOURCES
    submodels/kinematic/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModel.C
    submodels/kinematic/filmTurbulenceModel/filmTurbulenceModel/filmTurbulenceModelNew.C
    submodels/kinematic/filmTurbulenceModel/laminar/laminar.C
    )

list(APPEND SOURCES
    submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModel.C
    submodels/thermo/phaseChangeModel/phaseChangeModel/phaseChangeModelNew.C
    submodels/thermo/phaseChangeModel/noPhaseChange/noPhaseChange.C
    submodels/thermo/phaseChangeModel/standardPhaseChange/standardPhaseChange.C
    submodels/thermo/phaseChangeModel/solidification/solidification.C
    submodels/thermo/phaseChangeModel/waxSolventEvaporation/waxSolventEvaporation.C
    )

list(APPEND SOURCES
    submodels/thermo/massAbsorptionModel/massAbsorptionModel/massAbsorptionModel.C
    submodels/thermo/massAbsorptionModel/massAbsorptionModel/massAbsorptionModelNew.C
    submodels/thermo/massAbsorptionModel/noMassAbsorption/noMassAbsorption.C
    submodels/thermo/massAbsorptionModel/standardMassAbsorption/standardMassAbsorption.C
    )

list(APPEND SOURCES
    submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModel.C
    submodels/thermo/heatTransferModel/heatTransferModel/heatTransferModelNew.C
    submodels/thermo/heatTransferModel/constantHeatTransfer/constantHeatTransfer.C
    submodels/thermo/heatTransferModel/mappedConvectiveHeatTransfer/mappedConvectiveHeatTransfer.C
    )

list(APPEND SOURCES
    submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModel.C
    submodels/thermo/filmRadiationModel/filmRadiationModel/filmRadiationModelNew.C
    submodels/thermo/filmRadiationModel/noRadiation/noRadiation.C
    submodels/thermo/filmRadiationModel/constantRadiation/constantRadiation.C
    submodels/thermo/filmRadiationModel/standardRadiation/standardRadiation.C
    submodels/thermo/filmRadiationModel/rampingRadiation/rampingRadiation.C
    submodels/thermo/filmRadiationModel/standardIgnitorRadiation/standardIgnitorRadiation.C
    )

list(APPEND SOURCES
    submodels/thermo/filmViscosityModel/filmViscosityModel/filmViscosityModel.C
    submodels/thermo/filmViscosityModel/filmViscosityModel/filmViscosityModelNew.C
    submodels/thermo/filmViscosityModel/constantViscosity/constantViscosity.C
    submodels/thermo/filmViscosityModel/liquidViscosity/liquidViscosity.C
    submodels/thermo/filmViscosityModel/thixotropicViscosity/thixotropicViscosity.C
    submodels/thermo/filmViscosityModel/ArrheniusViscosity/ArrheniusViscosity.C
    submodels/thermo/filmViscosityModel/function1Viscosity/function1Viscosity.C
    submodels/thermo/filmViscosityModel/waxSolventViscosity/waxSolventViscosity.C
    )

list(APPEND SOURCES
    derivedFvPatchFields/inclinedFilmNusseltHeight/inclinedFilmNusseltHeightFvPatchScalarField.C
    derivedFvPatchFields/inclinedFilmNusseltInletVelocity/inclinedFilmNusseltInletVelocityFvPatchVectorField.C
    )

list(APPEND SOURCES
    derivedFvPatchFields/wallFunctions/alphatFilmWallFunction/alphatFilmWallFunctionFvPatchScalarField.C
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
#add_dependencies(${lib} pyrolysisModels)

target_include_directories(${lib} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/pyrolysisModels/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/solidChemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/chemistryModel/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basicSolidThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/distributionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/turbulenceModels
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/turbulenceModels/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/turbulenceModel/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/LES/lnInclude
    ${OpenFOAM_SRC_DIR}/TurbulenceModels/compressible/RAS/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_fluidThermophysicalModels}
    ${OF_specie}
    ${OF_thermophysicalProperties}
    ${OF_reactionThermophysicalModels}
    ${OF_SLGThermo}
    ${OF_distributionModels}
    ${OF_regionModels}
    ${OF_finiteVolume}
    ${OF_meshTools}
    ${OF_compressibleTransportModels}
    ${OF_compressibleTurbulenceModels}
    ${OF_pyrolysisModels}
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
    )
