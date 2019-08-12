add_executable(data data.cc)
target_link_libraries(data ScintillatingFibers)

add_executable(attenuation attenuation.cc)
target_link_libraries(attenuation ScintillatingFibers)

add_executable(energyres energyres.cc)
target_link_libraries(energyres ScintillatingFibers)
