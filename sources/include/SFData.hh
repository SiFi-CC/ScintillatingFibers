// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFData.hh                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFData_H_
#define __SFData_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TRandom.h"
#include <iostream>
#include <stdlib.h>

class SFData : public TObject{
  
private:
  int      fSeriesNo;
  int      fNpoints;
  TString  fFiber;
  TString  fDesc;
  TString  *fNames;
  double   *fPositions;
  TH1D     *fSpectra[9];
  TH1D     *fSpectrum;
  
  static TString fNames_1[9];
  static TString fNames_2[9];
  static TString fNames_3[9];
  static TString fNames_4[9];
  static TString fNames_5[9];
  static TString fNames_6[5];
  static TString fNames_7[1];
  static TString fNames_8[1];
  static TString fNames_9[2];
  
  static double fPositions_1[9];
  static double fPositions_2[9];
  static double fPositions_3[9];
  static double fPositions_4[9];
  static double fPositions_5[9];
  static double fPositions_6[5];
  
  TString GetSelection(int ch, TString type);
  
public:
  SFData();
  SFData(int seriesNo);
  ~SFData();
  
  bool     SetDetails(int seriesNo);
  TH1D*    GetSpectrum(int ch, TString type, TString cut, double position);
  TH1D**   GetSpectra(int ch, TString type, TString cut);
  void     Reset(void);
  void     Print(void);
  
  int      GetNpoints(void){ return fNpoints; };
  TString  GetFiber(void){ return fFiber; };
  TString  GetDescription(void){ return fDesc; };
  TString* GetNames(void) { return fNames; };
  double*  GetPositions(void) { return fPositions; };
  
  ClassDef(SFData,1)
  
};

#endif