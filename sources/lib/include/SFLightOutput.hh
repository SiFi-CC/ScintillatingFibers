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

/// Structure containing numerical results of the 
/// light output and light collection analysis.

struct SFLightResults{
    
    double fRes    = -1;     ///< Result (LO/LC) for summed channels
    double fResErr = -1;     ///< Uncertainty of the result for summed channel
    
    double fResCh0    = -1;  ///< Result (LO/LC) for channel 0 
    double fResCh0Err = -1;  ///< Uncertainty of the result for channel 0
     
    double fResCh1    = -1;  ///< Result (LO/LC) for channel 1
    double fResCh1Err = -1;  ///< Uncertainty of the result for channel 1
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
    
    SFLightResults fLightOutResults;
    SFLightResults fLightColResults;
    
public:
  SFLightOutput(int seriesNo);
  ~SFLightOutput();
  
  bool CalculateLightOut(void);
  bool CalculateLightOut(int ch);
  
  bool CalculateLightCol(void);
  bool CalculateLightCol(int ch);
  
  double GetCrossTalk(void);
  double GetPDE(void);
  
  SFLightResults GetLOResults(void) { return fLightOutResults; };
  SFLightResults GetLCResults(void) { return fLightColResults; };
  TCanvas*       GetInputData(void) { return fInputData; };
  
  TGraphErrors*        GetLightOutputGraph(void);
  TGraphErrors*        GetLightOutputGraph(int ch);
  TGraphErrors*        GetLightColGraph(void);
  TGraphErrors*        GetLightColGraph(int ch);
  
  std::vector <TH1D*>  GetSpectra(int ch);
  
  void Print(void);
  
  ClassDef(SFLightOutput,1)
};

#endif
