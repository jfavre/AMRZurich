add_executable(TestJet3D TestJet3D.cxx)
target_link_libraries(TestJet3D
	PUBLIC AMAZEReader
        PRIVATE VTK::RenderingOpenGL2
        )
        
add_executable(TestAMRValidity TestAMRValidity.cxx)
target_link_libraries(TestAMRValidity
	PUBLIC AMAZEReader
        )

