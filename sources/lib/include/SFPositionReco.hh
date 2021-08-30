// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFPositionReco.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#ifndef __SFPositionReco_H_
#define __SFPositionReco_H_ 1

#include "SFAttenuationModel.hh"
#include "SFAttenuation.hh"
#include "SFPositionRes.hh"
#include "SFData.hh"
#include "SFResults.hh"

#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TObject.h>

class SFPositionReco : public TObject
{

  private:
    int     fSeriesNo;
    SFData* fData;
    SFAttenuationModel* fModel;

    TGraphErrors* fMAttCh0CorrGraph;
    TGraphErrors* fMAttCh1CorrGraph;

    TGraphErrors* fMLRGraph;
    TGraphErrors* fMLRCorrGraph;

    TGraphErrors* fAGraph;
    
    TGraphErrors* fPosRecoGraph;
    TGraphErrors* fPosRecoCorrGraph;
    
    TGraphErrors* fPosResiduals;
    TGraphErrors* fPosResidualsCorr;
    
    TGraphErrors *fPosRecoDiff;
    TGraphErrors *fPosRecoDiffCorr;
    
    TGraphErrors* fPosResGraph;
    TGraphErrors* fPosResCorrGraph;
    
    TF2* fPlRecoFun;
    TF2* fPrRecoFun;
    
    TH1D* fPosRecoAll;
    
    std::vector<TH1D*> fRecoPositionsHist;
    std::vector<TH1D*> fRecoPositionsCorrHist;
    std::vector<TH1D*> fRecoPositionsUncertCorrHist;

    SFResults* fResultsExp;
    SFResults* fResultsCorr;

  public:
    SFPositionReco(int seriesNo);
    ~SFPositionReco();

    bool CalculateMLR(void);
    bool CalculatePosRecoCoefficients(void);
    bool PositionReco(void);

    std::vector<SFResults*> GetResults(void);
    std::vector<TH1D*>      GetPositionDistributions(TString type);
    std::vector<TH1D*>      GetErrorDistributions(void) { return fRecoPositionsUncertCorrHist; };

    void Print(void);

    ClassDef(SFPositionReco, 1)
};

#endif
