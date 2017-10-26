#=============================================================================
# CMake configuration file for the SUSPENSION_TEST_RIG demo.
# This example program works only with Irrlicht run-time visualization
#=============================================================================

IF(NOT ENABLE_MODULE_IRRLICHT)
    RETURN()
ENDIF()

IF(NOT ENABLE_MODULE_MKL)
    RETURN()
ENDIF()

#--------------------------------------------------------------
# List all model files for this demo

SET(DEMO demo_VEH_TrackTestRig)
SOURCE_GROUP("" FILES ${DEMO}.cpp)

#--------------------------------------------------------------
# Additional include directories

INCLUDE_DIRECTORIES(${CH_MKL_INCLUDES})
INCLUDE_DIRECTORIES(${CH_IRRLICHTINC})
INCLUDE_DIRECTORIES(CMAKE_CURRENT_SOURCE_DIR)

#--------------------------------------------------------------
# Add executable

MESSAGE(STATUS "...add ${DEMO}")

ADD_EXECUTABLE(${DEMO} ${DEMO}.cpp ${MODEL_FILES})
SET_TARGET_PROPERTIES(${DEMO} PROPERTIES
                      FOLDER demos 
                      COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_MKL_CXX_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_EXE} ${CH_MKL_LINK_FLAGS}")
TARGET_LINK_LIBRARIES(${DEMO} 
                      ChronoEngine
                      ChronoEngine_fea
                      ChronoEngine_mkl
                      ChronoEngine_irrlicht
                      ChronoEngine_vehicle
                      ChronoModels_vehicle)
ADD_DEPENDENCIES(${DEMO}
                ChronoEngine
                ChronoEngine_fea
                ChronoEngine_mkl
                ChronoEngine_irrlicht)
INSTALL(TARGETS ${DEMO} DESTINATION ${CH_INSTALL_DEMO})

