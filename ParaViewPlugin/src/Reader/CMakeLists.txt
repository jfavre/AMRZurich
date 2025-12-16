set(classes
  vtkCSCSAMRReader)

set(sources
  vtkAMAZEReader.cxx)
  
set(private_headers
  vtkCSCSAMRReader.h)

ADD_DEFINITIONS(-DH5_USE_16_API)

set_source_files_properties(
  vtkAMAZEReader
  WRAP_EXCLUDE
)

vtk_module_add_module(AMAZEReader
  CLASSES ${classes}
  SOURCES ${sources}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  CSCSAMRReader_Server.xml)
