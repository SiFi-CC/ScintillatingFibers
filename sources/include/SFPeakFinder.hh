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
#include <vector>
#include <iostream>

using namespace std;

class SFPeakFinder : public TObject{
 
private:
  TH1D    *fSpectrum;	///< Analuzed experimental spectrum
  TH1D    *fPeak;	///< Histogram with the chosen peak after background subtraction
  TString fPeakID;	///< Flag to identify which peak should be analyzed
  double  fPosition;	///< Position of the peak, determined as mean of Gaussian fit
  double  fPosErr;	///< Error of the peak position
  double  fSigma;	///< Energy resolution, determined as sigma of the Gaussian fit
  double  fSigErr;	///< Error of the peak sigma
  bool    fVerbose;	///< Print-outs level
  
  bool Fit(void);
  bool FindPeakRange(double &min, double &max);
  
public:
  SFPeakFinder();
  SFPeakFinder(TH1D *spectrum, TString peakID, bool verbose);
  SFPeakFinder(TH1D *spectrum, TString peakID);
  ~SFPeakFinder();
  
  bool SetSpectrum(TH1D *spectrum, TString peakID);
  void SetVerbLevel(bool verbose) { fVerbose=verbose; };
  void Print(void);
  void Clear(void);
  
  double GetPeakPosition(void) { return fPosition; };
  double GetPeakPosError(void) { return fPosErr; };
  double GetPeakSigma(void)    { return fSigma; };
  double GetPeakSigError(void) { return fSigErr; };
  vector <double> GetParameter(void);
  TH1D*  GetPeak(void)         { return fPeak; };
  
  ClassDef(SFPeakFinder,1);
};

#endif

