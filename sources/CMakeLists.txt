file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

ROOT_GENERATE_DICTIONARY(G__ScintillatingFibers ${headers} LINKDEF LinkDef.h)

add_library(ScintillatingFibers SHARED ${sources} G__ScintillatingFibers.cxx)
target_link_libraries(ScintillatingFibers ${ROOT_LIBRARIES} $ENV{DD6PATH}build/libDesktopDigitizer6.so libsqlite3.so)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/bin)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/bin)
