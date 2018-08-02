#--- CMakeLists.txt for SF package.
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ScintillatingFibers)


#--- Finding ROOT installation
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locating the ROOT package and defining a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Minuit)


#---Setting global variables e.g. where sources and includes are
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#---Defining useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories($ENV{DD6PATH}sources/include)
ROOT_GENERATE_DICTIONARY(G__ScintillatingFibers ${headers} LINKDEF LinkDef.h)


#---Creating a shared library with generated dictionary
add_library(ScintillatingFibers SHARED scifi.cc ${sources} G__ScintillatingFibers.cxx)
target_link_libraries(ScintillatingFibers ${ROOT_LIBRARIES} $ENV{DD6PATH}build/libDesktopDigitizer6.so libsqlite3.so)


#---Creating a main program using the library
add_executable(scifi scifi.cc)
target_link_libraries(scifi ScintillatingFibers)

add_executable(data data.cc)
target_link_libraries(data ScintillatingFibers)

add_executable(attenuation attenuation.cc)
target_link_libraries(attenuation ScintillatingFibers)

add_executable(timeres timeres.cc)
target_link_libraries(timeres ScintillatingFibers)

add_executable(tconst tconst.cc)
target_link_libraries(tconst ScintillatingFibers)

add_executable(lightout lightout.cc)
target_link_libraries(lightout ScintillatingFibers)

#---Locating Doxygen and generating documentation 
find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../sources/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
else(DOXYGEN_FOUND)
    message("Deoxygen must be installed to generate the documentation")
endif(DOXYGEN_FOUND)
