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
#include "TSpectrum.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "SFData.hh"
#include "SFTools.hh"
#include <vector>
#include <iostream>


/// Class designed to locate 511 keV in the recorded spectrum.
/// With the provided set of functions it is possible to find
/// the peak using ROOT's TSpectrum (SFPeakFinder::FindPeakSpectrum()),
/// by fitting sum of functions describing shape of spectrum 
/// (SFPeakFinder::FindPeakFit()) and by backgournd subtraction
/// (SFPeakFinder::FindPeakNoBackground()). All of the methods
/// provide information about peak position and sigma along with
/// their uncertainties. The latter method creates histogram with 
/// background subtracted 511 keV peak. Additionally, it is possible
/// to get peak range. 
///
/// Testing mode implemented for some methods in order to check fitting 
/// and background subtraction quality.

class SFPeakFinder : public TObject{
 
private:
  TH1D    *fSpectrum;  ///< Analyzed experimental spectrum
  TH1D    *fPeak;      ///< Histogram with the chosen peak after background subtraction
  double  fPosition;   ///< Position of the peak, determined as mean of Gaussian fit
  double  fPosErr;     ///< Error of the peak position
  double  fSigma;      ///< Energy resolution, determined as sigma of the Gaussian fit
  double  fSigErr;     ///< Error of the peak sigma
  double  fChi2NDF;    ///< Chi2/NDF of the fit Gassian fit of the peak 
  bool    fVerbose;    ///< Print-outs level
  bool    fTests;      ///< Flag for testing mode
  
public:
  SFPeakFinder();
  SFPeakFinder(TH1D *spectrum, bool verbose);
  SFPeakFinder(TH1D *spectrum);
  SFPeakFinder(TH1D *spectrum, bool verbose, bool tests);
  ~SFPeakFinder();
  
  bool                 FindPeakSpectrum(void);
  bool                 FindPeakRange(double &min, double &max);
  bool                 FindPeakFit(void);
  bool                 FindPeakNoBackground(void);
  std::vector <double> GetParameters(void); 
  void                 Print(void);
  
  void SetSpectrum(TH1D *spectrum);
  /// Sets level of print-outs: false - quiet, true - verbose.
  void SetVerbLevel(bool verbose) { fVerbose = verbose; };
  /// Sets flag for testing mode.
  void SetTests(bool tests) { fTests = tests; };
  /// Returns position of the peak.
  double GetPeakPosition(void) { return fPosition; };
  /// Returns error on peak position.
  double GetPeakPosError(void) { return fPosErr; };
  /// Returns sigma of the peak.
  double GetPeakSigma(void)    { return fSigma; };
  /// Returns error on peak's sigma
  double GetPeakSigError(void) { return fSigErr; };
  /// Returns the chi2 od the peak fit 
  double GetChi2NDF(void)      { return fChi2NDF; };
  /// Returns background-subtracted histogram. 
  TH1D*  GetPeak(void)         { return fPeak; };
  
  ClassDef(SFPeakFinder,1)
};

#endif
