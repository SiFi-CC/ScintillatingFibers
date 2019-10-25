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
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

class SFLightOutput : public TObject{
 
private:
    int     fSeriesNo;
    double  fPDE;
    double  fCrossTalk;
    
    double  fLightOut;
    double  fLightOutErr;
    double  fLightOutCh0;
    double  fLightOutCh0Err;
    double  fLightOutCh1;
    double  fLightOutCh1Err;
    
    double  fLightCol;
    double  fLightColErr;
    double  fLightColCh0;
    double  fLightColCh0Err;
    double  fLightColCh1;
    double  fLightColCh1Err;
    
    TGraphErrors *fLightOutGraph;
    TGraphErrors *fLightOutCh0Graph;
    TGraphErrors *fLightOutCh1Graph;
    
    TGraphErrors *fLightColGraph;
    TGraphErrors *fLightColCh0Graph;
    TGraphErrors *fLightColCh1Graph;
    
    std::vector <TH1D*> fSpectraCh0;
    std::vector <TH1D*> fSpectraCh1;
    std::vector <SFPeakFinder*> fPFCh0;
    std::vector <SFPeakFinder*> fPFCh1;
    
    SFData        *fData;
    SFAttenuation *fAtt;
    
public:
  SFLightOutput(int seriesNo);
  ~SFLightOutput();
  
  bool CalculateLightOut(void);
  bool CalculateLightOut(int ch);
  
  bool CalculateLightCol(void);
  bool CalculateLightCol(int ch);
  
  std::vector <double> GetLightOutput(void);
  std::vector <double> GetLightOutput(int ch);
  std::vector <double> GetLightCol(void);
  std::vector <double> GetLightCol(int ch);
  
  TGraphErrors*        GetLightOutputGraph(void);
  TGraphErrors*        GetLightOutputGraph(int ch);
  TGraphErrors*        GetLightColGraph(void);
  TGraphErrors*        GetLightColGraph(int ch);
  
  std::vector <TH1D*>  GetSpectra(int ch);
  
  void Print(void);
  
  ClassDef(SFLightOutput,1)
};

#endif
