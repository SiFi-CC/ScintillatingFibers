// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFPositionRecoModel.hh        *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#ifndef __SFPositionRecoModel_H_
#define __SFPositionRecoModel_H_ 1

#include "SFResults.hh"
#include "SFPeakFinder.hh"
#include "SFDrawCommands.hh"
#include "SFAttenuationModel.hh"

#include "SDDSamples.h"
#include "SLoop.h"
#include "SCategoryManager.h"

#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TH1D.h>

#include <tuple>
#include <optional>

namespace SFPositionRecoModel
{

    auto CalculateMLR(TGraphErrors *gAttCorrL,
                      TGraphErrors *gAttCorrR)
                      -> std::tuple<TGraphErrors*, TGraphErrors*>;
                      
    auto CalculateRecoCoefficients(TGraphErrors *gMLRCorr,
                                   double fiberLen)
                                   -> std::tuple<TGraphErrors*, double, double>;
    
    auto ReconstructPosition(SDDSamples *samples,
                             SFResults *results,
                             TMatrixD covMatrix,
                             TF1 *fPrReco,
                             TF1 *fPlReco,
                             std::vector<double> &parsForErrors,
                             float BL_sigma_cut,
                             float xmin, float xmax)
                             -> std::optional<std::tuple<double, double>>;
                             
    SFResults* ReconstructPositionDist(SFChAddr addr,
                                       SLoop *loop,
                                       TH1D* spectrum_av,
                                       TH1D* spectrum_l,
                                       TH1D* spectrum_r,
                                       SFResults *results,
                                       TMatrixD covMatrix,
                                       TF1 *fPrReco,
                                       TF1 *fPlReco,
                                       std::vector<double> &parsForErrors,
                                       TString path,
                                       double BL_sigma_cut,
                                       TString collimator);
    
    auto ReconstructPositionDistAll(SFChAddr addr, 
                                    std::vector<SLoop*> loop,
                                    std::vector<TH1D*> spectra_av,
                                    std::vector<TH1D*> spectra_l,
                                    std::vector<TH1D*> spectra_r,
                                    std::vector<double> positions,
                                    std::vector<TString> path,
                                    TMatrixD covMatrix,
                                    TF1 *fPrReco,
                                    TF1 *fPlReco,
                                    TGraphErrors *gAttCorrL,
                                    TGraphErrors *gAttCorrR,
                                    double BL_sigma_cut,
                                    double pos_uncert,
                                    double fiberLen,
                                    TString collimator,
                                    TString suffix)
                                    -> std::tuple<SFResults*, 
                                                  std::vector<TH1D*>,
                                                  std::vector<TH1D*>>;
                                                  
    SFResults* ReconstructPositionDistSum(SFChAddr addr, 
                                          std::vector<SLoop*> loop,
                                          std::vector<TH1D*> spectra_av,
                                          std::vector<TH1D*> spectra_l,
                                          std::vector<TH1D*> spectra_r,
                                          std::vector<double> positions,
                                          std::vector<TString> path,
                                          TMatrixD covMatrix,
                                          TF1 *fPrReco,
                                          TF1 *fPlReco,
                                          TGraphErrors *gAttCorrL,
                                          TGraphErrors *gAttCorrR,
                                          double BL_sigma_cut,
                                          double pos_uncert,
                                          TString collimator,
                                          TString suffix);

};

#endif
