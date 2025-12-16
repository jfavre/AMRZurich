add_executable(TestJet3D TestJet3D.cxx)
target_link_libraries(TestJet3D
	PUBLIC AMAZEReader
        PRIVATE VTK::RenderingOpenGL2)

