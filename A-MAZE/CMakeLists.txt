PROJECT(AMAZEReader)
cmake_minimum_required(VERSION 2.8.12)

SET(PLUGIN_NAME AMAZE)
SET(PLUGIN_VERSION "1.0")

if (NOT ParaView_SOURCE_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

vtk_module_load(vtkhdf5)
vtk_module_load(vtkIOAMR)

INCLUDE_DIRECTORIES(
  ${vtkhdf5_INCLUDE_DIRS}
  ${VTK_INCLUDE_DIRS}
)

ADD_DEFINITIONS(-DH5_USE_16_API)

SET_SOURCE_FILES_PROPERTIES(
  vtkAMRAmazeReaderInternal
  WRAP_EXCLUDE
)

ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION} 
                    SERVER_MANAGER_XML AMAZEReader.xml
                    SERVER_SOURCES vtkAMRAmazeReaderInternal.cxx
                    SERVER_MANAGER_SOURCES vtkAMRAmazeReader.cxx
#vtkSortAMRPoints.cxx
#vtkIonizationRayReader2.cxx

#vtkMetaAMRReader.cxx
#vtkExtractAMRStars.cxx
#                    GUI_RESOURCE_FILES CSCSAMRReader.xml
REQUIRED_ON_SERVER
)


#SUBDIRS(Testing)


include_directories(${VTK_INCLUDE_DIRS})
set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})

target_link_libraries(${PLUGIN_NAME} LINK_PRIVATE  ${vtkhdf5_LIBRARIES})

#ADD_EXECUTABLE(TestJet3D TestJet3D.cxx)
#TARGET_LINK_LIBRARIES(TestJet3D  vtkIOAMR vtkIOMPIParallel ${PLUGIN_NAME})
#set_target_properties(TestJet3D PROPERTIES COMPILE_FLAGS "-DAMAZE_EXPORTS")

ADD_EXECUTABLE(TestAMAZEReader TestAMAZEReader.cxx)
TARGET_LINK_LIBRARIES(TestAMAZEReader ${PLUGIN_NAME} ${VTK_LIBRARIES})
set_target_properties(TestAMAZEReader PROPERTIES COMPILE_FLAGS "-DAMAZE_EXPORTS")
