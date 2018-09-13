// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPeakFinder.hh            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFPeakFinder_H_
#define __SFPeakFinder_H_ 1
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "BGFit.hh"
#include "SFData.hh"
#include <vector>
#include <iostream>

using namespace std;

/// Class to perform background subtraction from charge spectra.
/// It is possible to subtract background from 511 keV ans 1270 keV peak,
/// depending on peak ID given in the constructor. Backround is described 
/// as exponential function. Class finds peak position and sigma along with 
/// their errors and creates histogram with background-subtracted peak.

class SFPeakFinder : public TObject{
 
private:
  TH1D    *fSpectrum;	///< Analyzed experimental spectrum
  TH1D    *fPeak;	///< Histogram with the chosen peak after background subtraction
  TString fPeakID;	///< Flag to identify which peak should be analyzed
  double  fPosition;	///< Position of the peak, determined as mean of Gaussian fit
  double  fPosErr;	///< Error of the peak position
  double  fSigma;	///< Energy resolution, determined as sigma of the Gaussian fit
  double  fSigErr;	///< Error of the peak sigma
  double  fChi2ndf;	///< Chi2 /ndf of the fit at the peak itself 
  bool    fVerbose;	///< Print-outs level
  
  bool    Fit(void);
  TString GetFiberMaterial(void);
  TString GetMeasureType(void);
  
public:
  SFPeakFinder();
  SFPeakFinder(TH1D *spectrum, TString peakID, bool verbose);
  SFPeakFinder(TH1D *spectrum, TString peakID);
  ~SFPeakFinder();
  
  bool SetSpectrum(TH1D *spectrum, TString peakID);
  bool FindPeakRange(double &min, double &max);
  vector <double> GetParameter(void);
  void Print(void);
  void Clear(void);
  
  /// Sets level of print-outs: false - quiet, true - verbose.
  void SetVerbLevel(bool verbose) { fVerbose=verbose; };
  /// Returns position of the peak.
  double GetPeakPosition(void) { return fPosition; };
  /// Returns error on peak position.
  double GetPeakPosError(void) { return fPosErr; };
  /// Returns sigma of the peak.
  double GetPeakSigma(void)    { return fSigma; };
  /// Returns error on peak's sigma
  double GetPeakSigError(void) { return fSigErr; };
  /// Returns the chi2 od the peak fit 
  double GetChi2ndf(void) { return fChi2ndf; };
  /// Returns background-subtracted histogram. 
  TH1D*  GetPeak(void)         { return fPeak; };
  
  ClassDef(SFPeakFinder,1)
};

#endif

