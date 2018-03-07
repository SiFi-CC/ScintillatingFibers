// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               SFMC.hh                 *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************


#ifndef __SFMC_H_
#define __SFMC_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

/// Class to access montecarlo data. It is possible to get single
/// histograms and vectors of histograms of requested simulated data. 

class SFMC : public TObject{
  
private:
  int            MCSeriesNo;		///< Series number of MC data
  int            MCNpoints;		///< Number of simulated measurements in the series
  TString        MCFiber;		///< Simulated Scintillating fiber 
  TString        MCDesc;			///< Description of the simulated data
  TString        *MCNames;		///< Array with names of simulated points
  double         *MCPositions;		///< Array with simulated positions of the radioactive source in mm
  TH1D           *MCSpectrum;		///< Single requested spectrum
  vector <TH1D*> MCSpectra;		///< Vector with all spectra from this simulated series (requested type e.g. Compton, ElectronBackground)
  
  static TString MCNames_1[9];
  
  static double MCPositions_1[9];
  
  TString GetSelection(TString type,TString cut);
  int     GetIndex(double position);
  
public:
  SFMC();
  SFMC(int seriesNo);
  ~SFMC();
  
  bool           SetDetails(int seriesNo);
  TH1D*          GetSpectrum(TString type, TString cut, double position);
  vector <TH1D*> GetSpectra(TString type, TString cut);
  void           Reset(void);
  void           Print(void);
  
  /// Returns number of measurements in the series
  int      GetNpoints(void){ return MCNpoints; };
  ///Returns fiber type
  TString  GetFiber(void){ return MCFiber; };
  ///Returns description of the series
  TString  GetDescription(void){ return MCDesc; };
  ///Returns an array containing names of all measurements in this series
  TString* GetNames(void) { return MCNames; };
  ///Returns an array containing source positions for all measurements in this series
  double*  GetPositions(void) { return MCPositions; };
  
  ClassDef(SFMC,1)
  
};

#endif
