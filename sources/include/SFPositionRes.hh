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
#include "TObject.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "SFData.hh"
#include "SFTools.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

class SFPositionRes : public TObject{
    
private:
    int                  fSeriesNo;
    double               fMLRPosResSeries;
    double               fMLRPosResSeriesErr;
    double               fYPosResSeries;
    double               fYPosResSeriesErr;
    
    SFData              *fData;
    SFAttenuation       *fAtt;
    
    TGraphErrors        *fMLRMeanVsPos;
    TGraphErrors        *fMLRPosResVsPos;
    TGraphErrors        *fYMeanVsPos;
    TGraphErrors        *fYPosResVsPos;
    
    std::vector <TH1D*>  fQRatios;
    std::vector <TH1D*>  fQRatiosCorr;
    std::vector <TH1D*>  fSpecCh0;
    std::vector <TH1D*>  fSpecCh1;
    
    std::vector <double> fMLRPosRes;
    std::vector <double> fMLRPosResErr;
    std::vector <double> fYPosRes;
    std::vector <double> fYPosResErr;
    
    bool LoadRatios(void);
    
public:
    SFPositionRes(int seriesNo);
    ~SFPositionRes();
    
    bool          AnalyzePositionResMLR(void); 
    bool          AnalyzePositionResY(void);
    
    TGraphErrors *GetMeanGraph(TString opt);
    TGraphErrors *GetPositionResGraph(TString opt);
    
    std::vector <TH1D*>  GetRatios(TString opt);
    std::vector <TH1D*>  GetSpectra(int ch);
    std::vector <double> GetPositionResSeries(TString opt);
    std::vector <double> GetPositionRes(TString opt);
    std::vector <double> GetPositionResError(TString opt);
    
    void Print();
    
    
  ClassDef(SFPositionRes,1)
};

#endif
