file(GLOB sources src/*.cc)
file(GLOB headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/include include/*.hh)
file(GLOB _headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} include/*.hh)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include ${DESKTOPDIGITIZER6_INCLUDE_DIR})

link_directories(${CMDLINEARGS_LIBRARY_DIR})

ROOT_GENERATE_DICTIONARY(G__ScintillatingFibers ${headers} LINKDEF LinkDef.h)

add_library(ScintillatingFibers SHARED ${sources} G__ScintillatingFibers.cxx)
target_link_libraries(ScintillatingFibers DesktopDigitizer6 SiFi Fibers sqlite3 ROOT::Core ROOT::Gpad SiFi::CmdLineArgs RT::FitterFactory)

set_target_properties(ScintillatingFibers PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${VERSION_MAJOR}
	PUBLIC_HEADER "${_headers}"
	INTERFACE_${CMAKE_PROJECT_NAME}_MAJOR_VERSION ${VERSION_MAJOR}
)


install(TARGETS ScintillatingFibers
        EXPORT ${CMAKE_PROJECT_NAME}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

