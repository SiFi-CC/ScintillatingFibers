// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFLightOutput.hh           *
// *             Jonas Kasper              *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFLightOutput_H_
#define __SFLightOutput_H_ 1

#include "SFPeakFinder.hh"
#include "SFResults.hh"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TFile.h>

#include <iostream>

namespace SFLightOutput
{
    SFResults* CalculateLightOut(TGraphErrors *gLightOutL, TGraphErrors *gLightOutR);
    
    SFResults* CalculateLightOut(TString side, std::vector<double> positions,
                                 std::vector<TH1D*> spectra, std::vector<TString> path,
                                 double pos_uncert, double fiberLen, double overvol, 
                                 TString fiber, TString sipm, SFResults* attResults);

    SFResults* CalculateLightCol(TGraphErrors *gLightColL, TGraphErrors *gLightColR);
    
    SFResults* CalculateLightCol(TString side, std::vector<double> positions,
                                 std::vector<TH1D*> spectra, std::vector<TString> path, double pos_uncert);

    double GetCrossTalk(TString SiPM, double overvol);
    
    auto GetPDE(TString fiber, TString SiPM, double overvol)
                -> std::tuple<double, TCanvas*>;
};

#endif
