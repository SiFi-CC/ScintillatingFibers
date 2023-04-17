add_executable(peakfin_pmi peakfin_pmi.cc)
target_link_libraries(peakfin_pmi ScintillatingFibers SmartFactory jsoncpp RootTools SiFi Fibers FitterFactory ROOT::Core ROOT::Gpad ROOT::Spectrum)

add_executable(attenuation_pmi attenuation_pmi.cc)
target_link_libraries(attenuation_pmi ScintillatingFibers SmartFactory jsoncpp RootTools SiFi Fibers FitterFactory ROOT::Core ROOT::Gpad ROOT::Spectrum)
