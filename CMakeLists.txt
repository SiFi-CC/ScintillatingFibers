cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ScintillatingFibers)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/Modules $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Spectrum)
set(ROOTSYS $ENV{ROOTSYS})

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/sources/include)
include_directories($ENV{DD6PATH}sources/include)
include_directories($ENV{DD6PATH}sources/include)

add_subdirectory(sources)

find_package(Doxygen)

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
