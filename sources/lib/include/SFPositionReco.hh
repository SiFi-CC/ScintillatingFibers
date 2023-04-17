// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionReco.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFPositionReco_H_
#define __SFPositionReco_H_ 1

#include "SDDSamples.h"
#include "SLoop.h"
#include "SCategoryManager.h"

#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFDrawCommands.hh"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TTree.h>
#include <TF1.h>
#include <TFile.h>

#include <iostream>
#include <tuple>
#include <optional>

namespace SFPositionReco
{
    auto ReconstructPosition(SDDSamples * samples, TF1* fun_mlr,
                             float BL_sigma_cut, float qmin, float qmax) 
                             -> std::optional<double>;
    
    SFResults* ReconstructPositionDist(SFChAddr addr, TH1D* spectrum_av,
                                       TF1* fun_mlr, TString path, 
                                       double BL_sigma_cut,
                                       TString collimator);
    
    auto ReconstructPositionDistAll(SFChAddr addr, 
                                    std::vector<TH1D*> spectrum_av, 
                                    std::vector<double> positions,
                                    std::vector<TString> path,
                                    TF1* fun_mlr, 
                                    double BL_sigma_cut,
                                    double pos_uncert, 
                                    TString collimator)
                                    -> std::tuple<SFResults*,std::vector<TH1D*>>;
                                    
    SFResults* ReconstructPositionDistSum(SFChAddr addr, 
                                          std::vector<TH1D*> spectrum_av,
                                          std::vector<double> positions,
                                          std::vector<TString> path,
                                          TF1* fun_mlr, double BL_sigma_cut, 
                                          TString collimator);
};

#endif
