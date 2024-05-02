########################################################
# cmake file for building Marlin ACTSTracking Package
# @author Karol Krizka
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( MarlinACTSTracking )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version

FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDCore )

FIND_PACKAGE( ROOT REQUIRED)

FIND_PACKAGE( Acts REQUIRED COMPONENTS Core PluginJson PluginTGeo )

INCLUDE(GNUInstallDirs)

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################
#ADD_SUBDIRECTORY(data)
#ADD_SUBDIRECTORY(example)

### LIBRARY #################################################################

# sources
CONFIGURE_FILE(src/config.h.in config.h)
ADD_LIBRARY(MarlinACTSTracking SHARED
  ./src/ACTSProcBase.cxx
  ./src/ACTSSeededCKFTrackingProc.cxx
  ./src/ACTSKFTrackingProc.cxx
  ./src/GeometryIdSelector.cxx
  ./src/GeometryIdMappingTool.cxx
  ./src/Helpers.cxx
  ./src/LUXEDuplicateRemoval.cxx
  ./src/TrackTruthProc.cxx
  )
TARGET_INCLUDE_DIRECTORIES(MarlinACTSTracking
  PRIVATE
    ./MarlinACTSTracking
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Marlin_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(MarlinACTSTracking
  ${ROOT_LIBRARIES} DD4hep::DDCore
  ActsCore ActsPluginTGeo ActsPluginJson)

INSTALL( TARGETS MarlinACTSTracking DESTINATION lib )

