set(classes
  vtkCSCSAMRReader
  vtkAMAZEReader
  vtkAMAZEFileSeriesReader)

set_source_files_properties(
  vtkAMAZEReader
  WRAP_EXCLUDE
)

vtk_module_add_module(AMAZEReader
  CLASSES ${classes}
  )
paraview_add_server_manager_xmls(
  XMLS  CSCSAMRReader_Server.xml
  )
