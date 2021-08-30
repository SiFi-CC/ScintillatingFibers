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
    SFAttenuationModel* fModel;

    TGraphErrors* fMAttCh0Graph;
    TGraphErrors* fMAttCh1Graph;
    TGraphErrors* fMAttCh0CorrGraph;
    TGraphErrors* fMAttCh1CorrGraph;
    
    TF2* fPlRecoFun;
    TF2* fPrRecoFun;
    
    std::vector<TH1D*> fEnergySpectra;
    std::vector<TH1D*> fEnergySpectraCorr;
    std::vector<TH1D*> fEnergyUncertDistCorr;
    
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
    std::vector<TH1D*>      GetErrorDistributions(void) { return fEnergyUncertDistCorr; };

    void Print(void);

    ClassDef(SFEnergyReco, 1)
};

#endif
