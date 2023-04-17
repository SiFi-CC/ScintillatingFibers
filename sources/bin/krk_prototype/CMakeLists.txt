add_executable(attenuation_krk attenuation_krk.cc)
target_link_libraries(attenuation_krk ScintillatingFibers SiFi Fibers FitterFactory ROOT::Core ROOT::Gpad ROOT::Spectrum)

add_executable(posres_krk posres_krk.cc)
target_link_libraries(posres_krk ScintillatingFibers SiFi Fibers FitterFactory ROOT::Core ROOT::Gpad ROOT::Spectrum)

add_executable(energyres_krk energyres_krk.cc)
target_link_libraries(energyres_krk ScintillatingFibers SiFi Fibers FitterFactory ROOT::Core ROOT::Gpad ROOT::Spectrum)
