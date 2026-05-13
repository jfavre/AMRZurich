set(classes
  vtkAMRAmazeReader
  )
set(private_classes
  vtkAMRAmazeReaderInternal
  )

vtk_module_add_module(AMAZEReader2
  CLASSES ${classes}
  PRIVATE_CLASSES ${private_classes}
  )
  
paraview_add_server_manager_xmls(
  XMLS AmazeReader_Server.xml
  )
