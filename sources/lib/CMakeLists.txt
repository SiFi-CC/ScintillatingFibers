file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

ROOT_GENERATE_DICTIONARY(G__ScintillatingFibers ${headers} LINKDEF LinkDef.h)

add_library(ScintillatingFibers SHARED ${sources} G__ScintillatingFibers.cxx)
target_link_libraries(ScintillatingFibers DesktopDigitizer6 sqlite3 ${ROOT_LIBRARIES} CmdLineArgs ${FITTERFACTORY_LIBRARIES})

set_target_properties(ScintillatingFibers PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${VERSION_MAJOR}
	PUBLIC_HEADER "${headers}"
	INTERFACE_${CMAKE_PROJECT_NAME}_MAJOR_VERSION ${VERSION_MAJOR}
	)
