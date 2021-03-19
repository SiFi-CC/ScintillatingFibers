add_executable(data data.cc)
target_link_libraries(data ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(attenuation attenuation.cc)
target_link_libraries(attenuation ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(energyres energyres.cc)
target_link_libraries(energyres ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(lightout lightout.cc)
target_link_libraries(lightout ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(timeres timeres.cc)
target_link_libraries(timeres ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(tconst tconst.cc)
target_link_libraries(tconst ${FITTERFACTORY_LIBRARIES} ScintillatingFibers SiFi Fibers)

add_executable(posres posres.cc)
target_link_libraries(posres ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(temp temp.cc)
target_link_libraries(temp ${FITTERFACTORY_LIBRARIES} ScintillatingFibers SiFi Fibers)

add_executable(peakfin peakfin.cc)
target_link_libraries(peakfin ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(model model.cc)
target_link_libraries(model ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(energyreco energyreco.cc)
target_link_libraries(energyreco ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(posreco posreco.cc)
target_link_libraries(posreco ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)

add_executable(stability stability.cc)
target_link_libraries(stability ${FITTERFACTORY_LIBRARIES} ScintillatingFibers MultiDimensionalFactory SmartFactory jsoncpp RootTools SiFi Fibers)
