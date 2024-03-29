cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ScintillatingFibers)

set (VERSION_MAJOR 2)
set (VERSION_MINOR 0)
set (VERSION_PATCH 0)
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

include(GNUInstallDirs)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/Modules $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS Core Spectrum)
#set(ROOTSYS $ENV{ROOTSYS})

find_package(CmdLineArgs 2.0.0 REQUIRED)
find_package(FitterFactory)
find_package(SiFi)
find_package(DesktopDigitizer6 REQUIRED)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${SIFI_INCLUDE_DIR})
include_directories(BEFORE sources/lib/include)

link_directories(${SIFI_LIBRARY_DIR})
link_directories(${FITTERFACTORY_LIBRARY_DIR})
link_directories(${DESKTOPDIGITIZER6_LIBRARY_DIR})

set(DEF_CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME})
set(CMAKE_INSTALL_CMAKEDIR ${DEF_CMAKE_INSTALL_CMAKEDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(MultiDimensionalFactory REQUIRED)

add_subdirectory(sources)
	
#install(TARGETS ScintillatingFibers
#	EXPORT ${CMAKE_PROJECT_NAME}Targets
#	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#)

#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/libScintillatingFibers_rdict.pcm
#	DESTINATION ${CMAKE_INSTALL_LIBDIR}
#)
	
install(TARGETS data attenuation energyres lightout peakfin posres stability tconst temp timeres model energyreco posreco  
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

export(PACKAGE ${CMAKE_PROJECT_NAME})

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/Config.cmake.in ${CMAKE_PROJECT_NAME}Config.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
	PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_BINDIR)

write_basic_package_version_file(
	${CMAKE_PROJECT_NAME}ConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

export(EXPORT ${CMAKE_PROJECT_NAME}Targets
	FILE ${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}Targets.cmake
)

install(EXPORT ${CMAKE_PROJECT_NAME}Targets
	FILE ${CMAKE_PROJECT_NAME}Targets.cmake
	DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake
	${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake
	DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
)

#find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/sources/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("DOXYGEN BUILD STARTED")
    
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "GENERATING API DOCUMENATATION WITH DOXYGEN"
        VERBATIM)
else(DOXYGEN_FOUND)
    message("DOXYGEN MUST BE INSTALLED TO GENERATE THE DOCUMENTATION")
endif(DOXYGEN_FOUND)
