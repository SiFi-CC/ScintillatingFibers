// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionRes.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFPositionRes_H_
#define __SFPositionRes_H_ 1

#include "SDDSamples.h"
#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFTools.hh"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TObject.h>
#include <TTree.h>

#include <iostream>
#include <tuple>
#include <optional>

class TH1;
class TF1;

namespace SFTools {

auto calculatePosition(SDDSamples * samples, TF1* mlr, float BL_sigma_cut, float qmin, float qmax) -> std::optional<float>;
    
auto calculatePositionResolutions(TH1* spectrum) -> std::tuple<float,std::vector<float>>;

auto calculateAveragePositionResolutions(std::vector<TH1*>) -> void;

};

namespace SFPositionRes
{
/*
    TGraphErrors* fPosVsMLRGraph;

    std::vector<TH1D*> fQRatios;
    std::vector<TH1D*> fPosRecoPol3Dist;
    std::vector<TH1D*> fPosRecoPol1Dist;
    std::vector<TH1D*> fSpecAv;

    bool LoadRatios(void);


    bool AnalyzePositionRes(void);

    std::vector<TH1D*> GetPositionRecoDist(TString type);
    std::vector<TH1D*> GetSpectra(void);*/
    
    double ReconstructPosition(SDDSamples * samples, TF1* fun_mlr,
                             float BL_sigma_cut, float qmin, float qmax);
    
    SFResults* ReconstructPositionDist(SFChAddr addr, SLoop *loop, TH1D* spectrum_av,
                                       TF1* fun_mlr, double BL_sigma_cut);
    
    auto ReconstructPositionDistAll(SFChAddr addr, SLoop* loop, TH1D* spectrum_av, 
                                    TF1* fun_mlr, std::vector<double> positions,
                                    double BL_sigma_cut, double pos_uncert, TString suffix)
                                    -> std::tuple<SFResults*,std::vector<TH1D*>>;

};

#endif
