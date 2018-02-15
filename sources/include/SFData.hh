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
#include "TProfile.h"
#include "DDSignal.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

/// Class to access experiemntal data. It is possible to get single
/// histograms and vectors of histograms of requested type. Single
/// signals and averaged signals are also possible to access. Cutting
/// functionality for all has been implemented.

class SFData : public TObject{
  
private:
  int            fSeriesNo;		///< Series number
  int            fNpoints;		///< Number of measurements in the series
  TString        fFiber;		///< Scintillating fiber type e.g. LuAG (1)
  TString        fDesc;			///< Description of the measurement series
  TString        *fNames;		///< Array with names of measurements
  double         *fPositions;		///< Array with positions of radioactive source in mm
  TH1D           *fSpectrum;		///< Single requested spectrum
  TProfile       *fSignalProfile;	///< Average of n requested signals
  TH1D           *fSignal;		///< Histogram of single chosen signal
  vector <TH1D*> fSpectra;		///< Vector with all spectra from this series (requested type e.g. fPE)
  vector <TH1D*> fRatios;		///< Vector of ratio histograms for all measurements (requested type)
  
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
  static double fPositions_7[1];
  static double fPositions_8[1];
  static double fPositions_9[2];
    
  TString GetSelection(int ch, TString type);
  int     GetIndex(double position);
  bool    InterpretCut(DDSignal *sig, TString cut);
  
public:
  SFData();
  SFData(int seriesNo);
  ~SFData();
  
  bool           SetDetails(int seriesNo);
  TH1D*          GetSpectrum(int ch, TString type, TString cut, double position);
  vector <TH1D*> GetSpectra(int ch, TString type, TString cut);
  vector <TH1D*> GetRatios(TString selection, TString cut);
  TProfile*      GetSignalAverage(int ch, double position, TString cut, int number, bool bl);
  TH1D*          GetSignal(int ch, double position, TString cut, int number, bool bl);  
  void           Reset(void);
  void           Print(void);
  
  /// Returns number of measurements in the series
  int      GetNpoints(void){ return fNpoints; };
  ///Returns fiber type
  TString  GetFiber(void){ return fFiber; };
  ///Returns description of the series
  TString  GetDescription(void){ return fDesc; };
  ///Returns an array containing names of all measurements in this series
  TString* GetNames(void) { return fNames; };
  ///Returns an array containing source positions for all measurements in this series
  double*  GetPositions(void) { return fPositions; };
  
  ClassDef(SFData,1)
  
};

#endif