add_executable(data data.cc)
target_link_libraries(data ScintillatingFibers)

add_executable(attenuation attenuation.cc)
target_link_libraries(attenuation ScintillatingFibers)

add_executable(energyres energyres.cc)
target_link_libraries(energyres ScintillatingFibers)

add_executable(lightout lightout.cc)
target_link_libraries(lightout ScintillatingFibers)

add_executable(timeres timeres.cc)
target_link_libraries(timeres ScintillatingFibers)

add_executable(tconst tconst.cc)
target_link_libraries(tconst ScintillatingFibers)

add_executable(posres posres.cc)
target_link_libraries(posres ScintillatingFibers)

add_executable(temp temp.cc)
target_link_libraries(temp ScintillatingFibers)

add_executable(peakfin peakfin.cc)
target_link_libraries(peakfin ScintillatingFibers)

add_executable(stability stability.cc)
target_link_libraries(stability ScintillatingFibers)
