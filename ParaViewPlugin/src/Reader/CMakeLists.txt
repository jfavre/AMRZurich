set(classes
  vtkAMRAmazeReader
  vtkAMRAmazeReaderInternal
  vtkAMAZEFileSeriesReader)

set_source_files_properties(
  vtkAMRAmazeReaderInternal
  WRAP_EXCLUDE
)

vtk_module_add_module(AMAZEReader
  CLASSES ${classes}
  )
paraview_add_server_manager_xmls(
  XMLS AmazeReader_Server.xml
  )
