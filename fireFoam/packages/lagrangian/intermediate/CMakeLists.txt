set(lib lagrangianIntermediate)

list(APPEND SOURCES
    clouds/baseClasses/kinematicCloud/kinematicCloud.C
    clouds/baseClasses/thermoCloud/thermoCloud.C
    clouds/baseClasses/reactingCloud/reactingCloud.C
    clouds/baseClasses/reactingMultiphaseCloud/reactingMultiphaseCloud.C
    clouds/baseClasses/reactingHeterogeneousCloud/reactingHeterogeneousCloud.C
    )

list(APPEND SOURCES
    parcels/derived/basicKinematicParcel/defineBasicKinematicParcel.C
    parcels/derived/basicKinematicParcel/makeBasicKinematicParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicKinematicCollidingParcel/defineBasicKinematicCollidingParcel.C
    parcels/derived/basicKinematicCollidingParcel/makeBasicKinematicCollidingParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicThermoParcel/defineBasicThermoParcel.C
    parcels/derived/basicThermoParcel/makeBasicThermoParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicReactingParcel/defineBasicReactingParcel.C
    parcels/derived/basicReactingParcel/makeBasicReactingParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicReactingMultiphaseParcel/defineBasicReactingMultiphaseParcel.C
    parcels/derived/basicReactingMultiphaseParcel/makeBasicReactingMultiphaseParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicHeterogeneousReactingParcel/defineBasicHeterogeneousReactingParcel.C
    parcels/derived/basicHeterogeneousReactingParcel/makeBasicHeterogeneousReactingParcelSubmodels.C
    )

list(APPEND SOURCES
    parcels/derived/basicKinematicMPPICParcel/defineBasicKinematicMPPICParcel.C
    parcels/derived/basicKinematicMPPICParcel/makeBasicKinematicMPPICParcelSubmodels.C
    )

list(APPEND SOURCES
    submodels/addOns/radiation/absorptionEmission/cloudAbsorptionEmission/cloudAbsorptionEmission.C
    submodels/addOns/radiation/scatter/cloudScatter/cloudScatter.C
    )

list(APPEND SOURCES
    submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.C
    submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionDataList.C
    )

list(APPEND SOURCES
    submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionData.C
    submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionDataIO.C
    submodels/Kinematic/InjectionModel/KinematicLookupTableInjection/kinematicParcelInjectionDataIOList.C
    submodels/Kinematic/InjectionModel/PatchInjection/patchInjectionBase.C
    )

list(APPEND SOURCES
    submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionData.C
    submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionDataIO.C
    submodels/Thermodynamic/InjectionModel/ThermoLookupTableInjection/thermoParcelInjectionDataIOList.C
    )

list(APPEND SOURCES
    submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionData.C
    submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionDataIO.C
    submodels/Reacting/InjectionModel/ReactingLookupTableInjection/reactingParcelInjectionDataIOList.C
    )

list(APPEND SOURCES
    submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionData.C
    submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIO.C
    submodels/ReactingMultiphase/InjectionModel/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIOList.C
    )

list(APPEND SOURCES
    submodels/MPPIC/ParticleStressModels/ParticleStressModel/ParticleStressModel.C
    submodels/MPPIC/ParticleStressModels/HarrisCrighton/HarrisCrighton.C
    submodels/MPPIC/ParticleStressModels/Lun/Lun.C
    submodels/MPPIC/ParticleStressModels/exponential/exponential.C
    )

list(APPEND SOURCES
    submodels/MPPIC/CorrectionLimitingMethods/CorrectionLimitingMethod/CorrectionLimitingMethod.C
    submodels/MPPIC/CorrectionLimitingMethods/noCorrectionLimiting/noCorrectionLimiting.C
    submodels/MPPIC/CorrectionLimitingMethods/absolute/absolute.C
    submodels/MPPIC/CorrectionLimitingMethods/relative/relative.C
    )

list(APPEND SOURCES
    submodels/MPPIC/TimeScaleModels/TimeScaleModel/TimeScaleModel.C
    submodels/MPPIC/TimeScaleModels/equilibrium/equilibrium.C
    submodels/MPPIC/TimeScaleModels/nonEquilibrium/nonEquilibrium.C
    submodels/MPPIC/TimeScaleModels/isotropic/isotropic.C
    )

list(APPEND SOURCES
    integrationScheme/integrationScheme/integrationScheme.C
    integrationScheme/integrationScheme/integrationSchemeNew.C
    integrationScheme/Euler/Euler.C
    integrationScheme/analytical/analytical.C
    )

list(APPEND SOURCES
    phaseProperties/phaseProperties/phaseProperties.C
    phaseProperties/phaseProperties/phasePropertiesIO.C
    phaseProperties/phasePropertiesList/phasePropertiesList.C
    )

list(APPEND SOURCES
    clouds/Templates/KinematicCloud/cloudSolution/cloudSolution.C
    submodels/CloudFunctionObjects/CloudFunctionObject/cloudFunctionObjectTools.C
    )

list(APPEND SOURCES
    submodels/MPPIC/AveragingMethods/makeAveragingMethods.C
    )

list(APPEND SOURCES
    conversion/ensight/ensightOutputCloud.C
    conversion/vtk/foamVtkLagrangianWriter.C
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

    ${CMAKE_CURRENT_SOURCE_DIR}/lnInclude

    ${OpenFOAM_SRC_DIR}/thermophysicalModels/radiation/lnInclude
    ${FF_PKG_DIR}/regionModels/surfaceFilmModels/lnInclude

    ${OpenFOAM_SRC_DIR}/lagrangian/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/lagrangian/distributionModels/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels/compressible/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/specie/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/basic/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/thermophysicalProperties/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/reactionThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/thermophysicalModels/SLGThermo/lnInclude
    ${OpenFOAM_SRC_DIR}/transportModels
    ${OpenFOAM_SRC_DIR}/regionModels/regionModel/lnInclude
    ${OpenFOAM_SRC_DIR}/dynamicMesh/lnInclude
    #${OpenFOAM_SRC_DIR}/dynamicFvMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/surfMesh/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteVolume/lnInclude
    ${OpenFOAM_SRC_DIR}/fileFormats/lnInclude
    ${OpenFOAM_SRC_DIR}/meshTools/lnInclude
    ${OpenFOAM_SRC_DIR}/parallel/distributed/lnInclude
    ${OpenFOAM_SRC_DIR}/regionFaModels/lnInclude
    ${OpenFOAM_SRC_DIR}/finiteArea/lnInclude
    ${OpenFOAM_SRC_DIR}/faOptions/lnInclude
    )

target_link_libraries(${lib} PUBLIC
    ${OF_lagrangian}
    ${OF_distributionModels}
    ${OF_specie}
    ${OF_compressibleTransportModels}
    ${OF_fluidThermophysicalModels}
    ${OF_thermophysicalProperties}
    ${OF_reactionThermophysicalModels}
    ${OF_SLGThermo}
    ${OF_incompressibleTransportModels}
    ${OF_regionModels}
    ${OF_dynamicMesh}
#    ${OF_dynamicFvMesh}
    ${OF_regionFaModels}
    ${OF_finiteArea}
    ${OF_finiteVolume}
    ${OF_fileFormats}
    ${OF_surfMesh}
    ${OF_meshTools}
    ${OF_radiation}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libsurfaceFilmModels.so
    )

install(TARGETS             ${lib}
    LIBRARY DESTINATION $ENV{FOAM_USER_LIBBIN}
    )

# Define source and destination directories  
set(SOURCE_DIR $ENV{FOAM_SRC}/lagrangian/intermediate)  
set(DEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})  
  
# Glob all files in the source directory  
file(GLOB_RECURSE FILES_TO_COPY RELATIVE ${SOURCE_DIR} ${SOURCE_DIR}/*)  
  
# Define a custom command to copy files only if different  
foreach(FILE ${FILES_TO_COPY})  
    # Create necessary directories  
    get_filename_component(DIR ${FILE} PATH)  
    add_custom_command(  
        OUTPUT ${DEST_DIR}/${FILE}  
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DEST_DIR}/${DIR}  
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SOURCE_DIR}/${FILE} ${DEST_DIR}/${FILE}  
        DEPENDS ${SOURCE_DIR}/${FILE}  
        COMMENT "Copying ${FILE} from ${SOURCE_DIR} to ${DEST_DIR}"  
    )  
endforeach()  
  
# Create a custom target that depends on the copied files  
add_custom_target(copy_files ALL DEPENDS ${FILES_TO_COPY}) 

# Ensure the copy_files target is built before your main target
add_dependencies(${lib} copy_files)
add_dependencies(${lib}_lnInclude copy_files)
