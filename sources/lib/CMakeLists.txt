#file(GLOB sources src/*.cc)
#file(GLOB headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/include include/*.hh)
#file(GLOB _headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} include/*.hh)

set(SRCS
    src/SFInfo.cc
    src/SFDrawCommands.cc
    src/SFResults.cc
    src/SFData.cc
    src/SFAttenuation.cc
    src/SFAttenuationModel.cc
    src/SFEnergyRes.cc
)

set(HDRS
    include/SFInfo.hh
    include/SFDrawCommands.hh
    include/SFResults.hh
    include/SFData.hh
    include/SFAttenuation.hh
    include/SFAttenuationModel.hh
    include/SFEnergyRes.hh
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

ROOT_GENERATE_DICTIONARY(G__ScintillatingFibers ${headers} LINKDEF LinkDef.h)

add_library(ScintillatingFibers SHARED ${SRCS} ${HDRS} G__ScintillatingFibers.cxx)
target_link_libraries(ScintillatingFibers DesktopDigitizer6 SiFi Fibers sqlite3 ${ROOT_LIBRARIES} CmdLineArgs ${FITTERFACTORY_LIBRARIES})

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

