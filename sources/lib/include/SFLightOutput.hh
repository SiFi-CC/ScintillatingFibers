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
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFTools.hh"
#include "TCanvas.h"
#include "TLatex.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

struct LightResults{
    
    double fRes    = -1;
    double fResErr = -1;
    
    double fResCh0    = -1;
    double fResCh0Err = -1;
    
    double fResCh1    = -1;
    double fResCh1Err = -1;
};

class SFLightOutput : public TObject{
 
private:
    int     fSeriesNo;
    double  fPDE;
    double  fCrossTalk;
    
    TGraphErrors *fLightOutGraph;
    TGraphErrors *fLightOutCh0Graph;
    TGraphErrors *fLightOutCh1Graph;
    
    TGraphErrors *fLightColGraph;
    TGraphErrors *fLightColCh0Graph;
    TGraphErrors *fLightColCh1Graph;
    
    TCanvas       *fInputData;
    
    std::vector <TH1D*> fSpectraCh0;
    std::vector <TH1D*> fSpectraCh1;
    std::vector <SFPeakFinder*> fPFCh0;
    std::vector <SFPeakFinder*> fPFCh1;
    
    SFData        *fData;
    SFAttenuation *fAtt;
    
    LightResults fLightOutResults;
    LightResults fLightColResults;
    
public:
  SFLightOutput(int seriesNo);
  ~SFLightOutput();
  
  bool CalculateLightOut(void);
  bool CalculateLightOut(int ch);
  
  bool CalculateLightCol(void);
  bool CalculateLightCol(int ch);
  
  double GetCrossTalk(void);
  double GetPDE(void);
  
  LightResults GetLOResults(void) { return fLightOutResults; };
  LightResults GetLCResults(void) { return fLightColResults; };
  TCanvas*     GetInputData(void) { return fInputData; };
  
  TGraphErrors*        GetLightOutputGraph(void);
  TGraphErrors*        GetLightOutputGraph(int ch);
  TGraphErrors*        GetLightColGraph(void);
  TGraphErrors*        GetLightColGraph(int ch);
  
  std::vector <TH1D*>  GetSpectra(int ch);
  
  void Print(void);
  
  ClassDef(SFLightOutput,1)
};

#endif
