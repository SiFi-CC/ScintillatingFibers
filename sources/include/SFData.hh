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
#include "TH2D.h"
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
  TString	 fThreshold;		///< Flag to identify which tree should be opened.
					///< Default value - "ft"
  TString        *fNames;		///< Array with names of measurements
  double         *fPositions;		///< Array with positions of radioactive source in mm
  TH1D           *fSpectrum;		///< Single requested spectrum
  TH1D           *fHist;		///< Custom histogram of requested type
  TH2D           *fHist2D;		///< 2D correlation histogram of requested type 
  TProfile       *fSignalProfile;	///< Average of n requested signals
  TH1D           *fSignal;		///< Histogram of single chosen signal
  vector <TH1D*> fSpectra;		///< Vector with all spectra from this series (requested type e.g. fPE)
  vector <TH1D*> fHists;		///< Vector of custom histograms for all measurements (requested type)
  vector <TH2D*> fHists2D;		///< Vector of 2D correlation histograms for all measurements (requested type)
  
  static TString fNames_1[9];		///< Measurements names, series 1
  static TString fNames_2[9];		///< Measurements names, series 2
  static TString fNames_3[9];		///< Measurements names, series 3
  static TString fNames_4[9];		///< Measurements names, series 4
  static TString fNames_5[9];		///< Measurements names, series 5
  static TString fNames_6[5];		///< Measurements names, series 6
  static TString fNames_7[1];		///< Measurements names, series 7
  static TString fNames_8[1];		///< Measurements names, series 8
  static TString fNames_9[2];		///< Measurements names, series 9
  
  static double fPositions_1[9];	///< Source positions for series 1 [mm]
  static double fPositions_2[9];	///< Source positions for series 2 [mm]
  static double fPositions_3[9];	///< Source positions for series 3 [mm]
  static double fPositions_4[9];	///< SOurce positions for series 4 [mm]
  static double fPositions_5[9];	///< Source positions for series 5 [mm]
  static double fPositions_6[5];	///< SOurce positions for series 6 [mm]
  static double fPositions_7[1];	///< Source positions for series 7 [mm]
  static double fPositions_8[1];	///< Source positions for series 8 [mm]
  static double fPositions_9[2];	///< Source positions for series 9 [mm]
    
  TString GetSelection(int ch, TString type);
  TString GetSelectionCustom(TString selection);
  int     GetIndex(double position);
  bool    InterpretCut(DDSignal *sig, TString cut);
  
public:
  SFData();
  SFData(int seriesNo);
  SFData(int seriesNo, TString threshold);
  ~SFData();
  
  bool           SetDetails(int seriesNo);
  bool           SetThreshold(TString threshold);
  TH1D*          GetSpectrum(int ch, TString type, TString cut, double position);
  TH1D*          GetCustomHistogram(TString selection, TString cut, double position);
  TH2D*          GetCorrHistogram(TString selection, TString cut, double position);
  vector <TH1D*> GetSpectra(int ch, TString type, TString cut);
  vector <TH1D*> GetCustomHistograms(TString type, TString cut);
  vector <TH2D*> GetCorrHistograms(TString type, TString cut);
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