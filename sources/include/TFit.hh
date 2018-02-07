// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              TFit.hh                  *
// *            Jonas Kasper               *
// *    kasper@physik.rwth-aachen.de       *        
// *          Created in 2018              *
// *                                       *
// *****************************************


#ifndef __TFit_H_
#define __TFit_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TH1D.h"
#include "TROOT.h"
#include <iostream>
#include <stdlib.h>

class TFit : public TObject{
  
private:
  
  TString name;
  Double_t position;
  Double_t resolution;
  
  // Spectrum that is analysed
  TH1D* spectrum;
  
  // Needed Templates
  TH1D* background;
  TH1D* pp511;
  TH1D* compton511;
  TH1D* pp1200;
  TH1D* compton1200;
  
  // Fitting function called in Constructer
  void Fitting();
  
public:
  TFit();
  TFit(TH1D* Spectrum,TH1D* Background);
  TFit(TH1D* Spectrum,TH1D* Background, TString Name);
  ~TFit();
  
  // Getter for fitted spectrum
  TCanvas* DrawFittedSpectrum();
  // Getter of fit data
  Double_t* GetFitData();
  
  ClassDef(TFit,1)
  
};

#endif