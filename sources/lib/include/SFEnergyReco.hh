// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFEnergyReco.hh            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#ifndef __SFEnergyReco_H_
#define __SFEnergyReco_H_ 1

#include "SFAttenuationModel.hh"
#include "SFData.hh"
#include "SFResults.hh"

#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TObject.h>

class SFEnergyReco : public TObject
{

  private:
    int     fSeriesNo;
    SFData* fData;

    TGraphErrors* fMAttCh0Graph;
    TGraphErrors* fMAttCh1Graph;
    TGraphErrors* fMAttCh0CorrGraph;
    TGraphErrors* fMAttCh1CorrGraph;
    
    TF2* fPlRecoFun;
    TF2* fPrRecoFun;

    TGraphErrors* fEnergyRecoGraph;
    TGraphErrors* fEnergyRecoCorrGraph;
    
    TGraphErrors* fEnergyAlphaGraph;
    TGraphErrors* fEnergyAlphaCorrGraph;
    
    TGraphErrors* fEnergyRecoSpecGraph;
    TGraphErrors* fEnergyRecoSpecCorrGraph;
    
    TGraphErrors* fEnergyResGraph;
    TGraphErrors* fEnergyResCorrGraph;
    
    std::vector<TH1D*> fEnergySpectra;
    std::vector<TH1D*> fEnergySpectraCorr;
    
    SFResults* fResultsExp;
    SFResults* fResultsCorr;
    
    double fEref;
    
  public:
    SFEnergyReco(int seriesNo);
    ~SFEnergyReco();

    bool CalculateAlpha(void);
    bool EnergyReco(void);
    bool EnergyRecoByEvent(void);

    std::vector<SFResults*> GetResults(void);
    std::vector<TH1D*>      GetEnergySpectra(TString type);

    void Print(void);

    ClassDef(SFEnergyReco, 1)
};

#endif
