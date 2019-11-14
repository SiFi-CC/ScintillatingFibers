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
#include "TTree.h"
#include "DDSignal.hh"
#include "SFData.hh"
#include "SFTools.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>
#include "TCanvas.h"

class SFPositionRes : public TObject{
    
private:
    int                  fSeriesNo;
    double               fPosResSeries;
    double               fPosResSeriesErr;
    
    SFData              *fData;
    SFAttenuation       *fAtt;
    
    TGraphErrors        *fPosRecoVsPos;
    TGraphErrors        *fPosResVsPos;
    TGraphErrors        *fMLRvsPos;
    
    std::vector <TH1D*>  fQRatios;
    std::vector <TH1D*>  fPosRecoDist;
    std::vector <TH1D*>  fSpecAv;
    
    std::vector <double> fPosRes;
    std::vector <double> fPosResErr;
    std::vector <double> fPosReco;
    std::vector <double> fPosRecoErr;
    
    bool LoadRatios(void);
    
public:
    SFPositionRes(int seriesNo);
    ~SFPositionRes();
    
    bool          AnalyzePositionRes(void); 
    
    TGraphErrors *GetPositionRecoGraph(void);
    TGraphErrors *GetPositionResGraph(void);
    TGraphErrors *GetAttenuationCurve(void);
    
    std::vector <TH1D*>  GetRatios(void);
    std::vector <TH1D*>  GetPositionRecoDist(void);
    std::vector <TH1D*>  GetSpectra(void);
    std::vector <double> GetPositionResSeries(void);
    std::vector <double> GetPositionRes(void);
    std::vector <double> GetPositionResError(void);
    std::vector <double> GetPositionReco(void);
    std::vector <double> GetPositionRecoError(void);
    
    void Print();
    
    
  ClassDef(SFPositionRes,1)
};

#endif
