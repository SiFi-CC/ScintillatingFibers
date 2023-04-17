// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFEnergyRes.hh              *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************
#ifndef __SFEnergyRes_H_
#define __SFEnergyRes_H_ 1

#include "SFPeakFinder.hh"
#include "SFResults.hh"

#include <TF1.h>
#include <TGraphErrors.h>

#include <iostream>
#include <cmath>
#include <memory>

namespace SFEnergyRes
{
    SFResults* CalculateEnergyResSingle(TH1D* h, TString path);
    SFResults* CalculateEnergyResSeries(TString side, double pos_uncert,
                                        std::vector<double> positions,
                                        std::vector<TH1D*> spectra,
                                        std::vector<TString> path);
};

#endif
