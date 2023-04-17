// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             SFCountMap.hh             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************
#ifndef __SFCountMap_H_
#define __SFCounrMap_H_ 1

#include "SFResults.hh"
#include "SFTools.hh"

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>

#include <iostream>

/// This class is meant for analysis of measurements recorded with
/// the entire detector module (test bench: PMI, description: Module 
/// series). It allows to plot number of counts registered in the 
/// fiber as a function of the fiber number. Two fiber sides are plotted 
/// separately. In order to obtain number of counts the charge spectra 
/// are integrated with the ROOT's method IntegralAndError(). 

namespace SFCountMap
{
    SFResults* DrawCounts(TString side, std::vector<TH1D*> spectra,
                          std::vector<int> fiberNo);
};

#endif /* __SFCountMap_H_ */
