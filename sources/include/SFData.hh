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
#include "TProfile.h"
#include "TVectorT.h"
#include "DDSignal.hh"
#include "SFDrawCommands.hh"
#include "SFTools.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sqlite3.h>

/// Class to access experiemntal data. Data is acessed via SQLite3 
/// data base. It is possible to get single histograms and vectors 
/// of histograms of requested type. Single signals and averaged 
/// signals are also possible to access. Cutting functionality for 
/// all has been implemented.

class SFData : public TObject{
    
private:
  int              fSeriesNo;        ///< Series number
  int              fNpoints;         ///< Number of measurements in the series
  TString          fFiber;           ///< Scintillating fiber type e.g. LuAG (1)
  TString          fSource;          ///< Type of the radioactive source
  TString          fCollimator;      ///< Type of the used collimator: Lead or Electronic
  TString          fDesc;            ///< Description of the measurement series
  TString          fTestBench;       ///< Type of test bench: Krakow/Aachen/Simulation
  TH1D             *fSpectrum;       ///< Single requested spectrum
  TH1D             *fHist;           ///< Custom histogram of requested type
  TH2D             *fHist2D;         ///< 2D correlation histogram of requested type 
  TProfile         *fSignalProfile;  ///< Average of n requested signals
  TH1D             *fSignal;         ///< Histogram of single chosen signal
  std::vector <TString> fNames;      ///< Vector with names of measurements
  std::vector <double>  fPositions;  ///< Vector with positions of radioactive source in mm
  std::vector <int>  fTimes;         ///< Vector with times of measurement in s
  std::vector <int>  fStart;         ///< Vector containing starting times of measurements (in UNIX time)
  std::vector <int>  fStop;          ///< Vector containing stopping times of measurements (in UNIX time)
  std::vector <TH1D*>   fSpectra;    ///< Vector with all spectra from this series (requested type e.g. fPE)
  std::vector <TH1D*>   fHists;      ///< Vector of custom histograms for all measurements (requested type)
  std::vector <TH2D*>   fHists2D;    ///< Vector of 2D correlation histograms for all measurements (requested type)
  sqlite3               *fDB;        ///< SQLite3 data base
  
  int  gUnique = 0.;                 ///< Unique flag to identify histograms
  
  bool      InterpretCut(DDSignal *sig, TString cut);
  TProfile* GetSignalAverageKrakow(int ch, double position, TString cut, int number, bool bl);
  TProfile* GetSignalAverageAachen(int ch, double position, TString cut, int number);
  TH1D*     GetSignalKrakow(int ch, double position, TString cut, int number, bool bl);
  TH1D*     GetSignalAachen(int ch, double position, TString cut, int number);
  
public:
  SFData();
  SFData(int seriesNo);
  ~SFData();
  
  bool                OpenDataBase(TString name);
  bool                SetDetails(int seriesNo);
  TH1D*               GetSpectrum(int ch, SFSelectionType sel_type, TString cut, double position);
  TH1D*               GetCustomHistogram(SFSelectionType sel_type, TString cut, double position);
  TH2D*               GetCorrHistogram(SFSelectionType sel_type, TString cut, double position);
  std::vector <TH1D*> GetSpectra(int ch, SFSelectionType sel_type, TString cut);
  std::vector <TH1D*> GetCustomHistograms(SFSelectionType sel_type, TString cut);
  std::vector <TH2D*> GetCorrHistograms(SFSelectionType sel_type, TString cut);
  TProfile*           GetSignalAverage(int ch, double position, TString cut, int number, bool bl);
  TH1D*               GetSignal(int ch, double position, TString cut, int number, bool bl);
  void                Print(void);
  
  /// Returns number of measurements in the series
  int      GetNpoints(void){ return fNpoints; };
  /// Returns fiber type
  TString  GetFiber(void){ return fFiber; };
  /// Returns description of the series
  TString  GetDescription(void){ return fDesc; };
  /// Returns collimator type
  TString  GetCollimator(void){ return fCollimator; };
  /// Returns source type
  TString  GetSource(void){ return fSource; };
  /// Returns test bench type
  TString GetTestBench(void){ return fTestBench; };
  /// Returns a vector containing names of all measurements in this series
  std::vector <TString> GetNames(void){ return fNames; };
  /// Returns a vector containing source positions for all measurements in this series
  std::vector <double>  GetPositions(void){ return fPositions; };
  /// Returns a vector containing times of measurements.
  std::vector <int> GetTimes(void){ return fTimes; };
  /// Returns a vector containing starting times of measurements.
  std::vector <int> GetStartTimes(void){ return fStart; };
  /// Returns a vector containing stopping times of measurements.
  std::vector <int> GetStopTimes(void){ return fStop; };
  
  ClassDef(SFData,1)
  
};

#endif
