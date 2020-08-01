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
#include "SDDSamples.h"
#include "SFData.hh"
#include "SFTools.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>
#include "TCanvas.h"

struct PositionResResults{
    
    double fPosRes    = -1;
    double fPosResErr = -1;
    
    std::vector <double> fPosResAll;
    std::vector <double> fPosResAllErr;
    
    std::vector <double> fPosReco;
    std::vector <double> fPosRecoErr;    
};

class SFPositionRes : public TObject{
    
private:
    int                  fSeriesNo;
    SFData              *fData;
    SFAttenuation       *fAtt;
    
    TGraphErrors        *fPosRecoVsPos;
    TGraphErrors        *fPosResVsPos;
    TGraphErrors        *fMLRvsPos;
    TGraphErrors        *fResiduals;
    
    std::vector <TH1D*>  fQRatios;
    std::vector <TH1D*>  fPosRecoDist;
    std::vector <TH1D*>  fSpecAv;
    
    PositionResResults   fResults;
    
    bool LoadRatios(void);
    
public:
    SFPositionRes(int seriesNo);
    ~SFPositionRes();
    
    bool          AnalyzePositionRes(void); 
    
    TGraphErrors *GetPositionRecoGraph(void);
    TGraphErrors *GetPositionResGraph(void);
    TGraphErrors *GetAttenuationCurve(void);
    TGraphErrors *GetResiduals(void);
    
    std::vector <TH1D*>  GetRatios(void);
    std::vector <TH1D*>  GetPositionRecoDist(void);
    std::vector <TH1D*>  GetSpectra(void);
    
    PositionResResults   GetResults(void) { return fResults; };
    
    void Print();
    
    ClassDef(SFPositionRes,1)
};

#endif
