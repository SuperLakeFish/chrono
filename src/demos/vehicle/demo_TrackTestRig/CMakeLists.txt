#=============================================================================
# CMake configuration file for the SUSPENSION_TEST_RIG demo.
# This example program works only with Irrlicht run-time visualization
#=============================================================================

if(NOT ENABLE_MODULE_IRRLICHT)
    return()
endif()

#--------------------------------------------------------------
# List all model files for this demo

set(DEMOS
    demo_VEH_TrackTestRig
)

#--------------------------------------------------------------
# Compiler and linker flags

include_directories(${CH_IRRLICHTINC})

set(COMPILER_FLAGS "${CH_CXX_FLAGS}  ${CH_IRRLICHT_CXX_FLAGS}")
set(LINKER_FLAGS "${CH_LINKERFLAG_EXE}")
list(APPEND LIBS "")

if(ENABLE_MODULE_FEA)
    list(APPEND LIBS "ChronoEngine_fea")
endif()

if(ENABLE_MODULE_MKL)
    include_directories(${CH_MKL_INCLUDES})
    set(COMPILER_FLAGS "${COMPILER_FLAGS} ${CH_MKL_CXX_FLAGS}")
    set(LINKER_FLAGS "${LINKER_FLAGS} ${CH_MKL_LINK_FLAGS}")
    list(APPEND LIBS "ChronoEngine_mkl")
endif()

#--------------------------------------------------------------
# Add executables

foreach(DEMO ${DEMOS})

    message(STATUS "...add ${DEMO}")
    
    add_executable(${DEMO} ${DEMO}.cpp)
    source_group("" FILES ${DEMO}.cpp)
    set_target_properties(${DEMO} PROPERTIES
                          FOLDER demos 
                          COMPILE_FLAGS "${COMPILER_FLAGS}"
                          LINK_FLAGS "${LINKER_FLAGS}")
    target_link_libraries(${DEMO} 
                          ChronoEngine
                          ChronoEngine_irrlicht
                          ChronoEngine_vehicle
                          ChronoModels_vehicle
                          ${LIBS})
    install(TARGETS ${DEMO} DESTINATION ${CH_INSTALL_DEMO})

endforeach()
