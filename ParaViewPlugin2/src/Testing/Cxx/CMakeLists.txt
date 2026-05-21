add_executable(TestJet3D TestJet3D.cxx)
target_link_libraries(TestJet3D
	PUBLIC AMAZEReader2
        PRIVATE VTK::RenderingOpenGL2)

add_executable(TestAMRValidity TestAMRValidity.cxx)
target_link_libraries(TestAMRValidity
	PUBLIC AMAZEReader2
        )
