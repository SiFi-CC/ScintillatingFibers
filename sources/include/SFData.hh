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
#include <sqlite3.h>

using namespace std;

/// Class to access experiemntal data. Data is acessed via SQLite3 
/// data base. It is possible to get single histograms and vectors 
/// of histograms of requested type. Single signals and averaged 
/// signals are also possible to access. Cutting functionality for 
/// all has been implemented.

class SFData : public TObject{
  
private:
  int            fSeriesNo;		///< Series number
  int            fNpoints;		///< Number of measurements in the series
  TString        fFiber;		///< Scintillating fiber type e.g. LuAG (1)
  TString        fSource;		///< Type of the radioactive source
  TString        fDesc;			///< Description of the measurement series
  TString	 fThreshold;		///< Flag to identify which tree should be opened.
					///< Default value - "ft"
  vector <TString> fNames;		///< Vector with names of measurements
  vector <double>  fPositions;		///< Vector with positions of radioactive source in mm
  vector <double>  fTimes;		///< Vector with times of measurement in s
  TH1D             *fSpectrum;		///< Single requested spectrum
  TH1D             *fHist;		///< Custom histogram of requested type
  TH2D             *fHist2D;		///< 2D correlation histogram of requested type 
  TProfile         *fSignalProfile;	///< Average of n requested signals
  TH1D             *fSignal;		///< Histogram of single chosen signal
  vector <TH1D*>   fSpectra;		///< Vector with all spectra from this series (requested type e.g. fPE)
  vector <TH1D*>   fHists;		///< Vector of custom histograms for all measurements (requested type)
  vector <TH2D*>   fHists2D;		///< Vector of 2D correlation histograms for all measurements (requested type)
  sqlite3          *fDB;		///< SQLite3 data base
    
  TString    GetSelection(int ch, TString type);
  TString    GetSelectionCustom(TString selection);
  int        GetIndex(double position);
  bool       InterpretCut(DDSignal *sig, TString cut);
  
public:
  SFData();
  SFData(int seriesNo);
  SFData(int seriesNo, TString threshold);
  ~SFData();
  
  bool           OpenDataBase(TString name);
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
  ///Returns source type
  TString  GetSource(void){ return fSource; };
  ///Returns a vector containing names of all measurements in this series
  vector <TString> GetNames(void){ return fNames; };
  ///Returns a vector containing source positions for all measurements in this series
  vector <double>  GetPositions(void){ return fPositions; };
  ///Returns a vector containing times of measurements.
  vector <double> GetTimes(void){ return fTimes; };
  
  ClassDef(SFData,1)
  
};

#endif