// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFAttenuation.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFAttenuation_H_
#define __SFAttenuation_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include <iostream>

using namespace std;

/// Class to determine attenuation length. This class is suitable only for experimental 
/// series with different positions of source, i.e. 1, 2, 3, 4 and 5. Two methods of 
/// attenuation length determination are available: AttAveragedCh() - based on Pauwels 
/// et al., JINST 8 (2013) P09019, where averaged signal from both channels is analyzed, 
/// and AttSeparateCh() - where attenuation length is calulated for each channel separately.

class SFAttenuation : public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  double fAttnLen;		///< Attenuation length determined with averaged channels method [mm]
  double fAttnErr;		///< Error on attenuation length fAttnLen [mm]
  vector <TH1D*> fRatios;	///< Vector containing histograms of signal ratios from both channels, for whole series 
  TGraphErrors *fAttnGraph;	///< Attenuation graph i.e. ln(M_{FB}) vs. source position
  double fAttnLenCh0;		///< Attenuation length for channel 0
  double fAttnLenCh1;		///< Attenuation length for channel 1
  double fAttnErrCh0;		///< Error on attenuation length for channel 0
  double fAttnErrCh1;		///< Error on attenuation length for channel 1
  TGraphErrors *fAttnGraphCh0;	///< Attenuation graph for channel 0
  TGraphErrors *fAttnGraphCh1;	///< Attenuation graph for channel 1
  vector <TH1D*> fSpectraCh0;	///< Vector containing charge spectra from channel 0
  vector <TH1D*> fSpectraCh1;	///< Vector containing charche spectra from channel 1
  vector <TH1D*> fPeaksCh0;	///< Vector containing 511 keV peaks, channel 0
  vector <TH1D*> fPeaksCh1;	///< Vector containing 511 keV peaks, channel 1
  
public:
  SFAttenuation();
  SFAttenuation(int seriesNo);
  ~SFAttenuation();
  
  bool            AttAveragedCh(void);
  vector <TH1D*>  GetRatios(void);
  vector <double> GetAttenuation(void);
  TGraphErrors*   GetAttGraph(void);
  double          GetAttLength(void);
  double          GetAttError(void);
  
  bool            AttSeparateCh(int ch);
  vector <TH1D*>  GetSpectra(int ch);
  vector <TH1D*>  GetPeaks(int ch);
  vector <double> GetAttenuation(int ch);
  TGraphErrors*   GetAttGraph(int ch);
  double          GetAttLength(int ch);
  double          GetAttError(int ch);
  
  void Print(void);
  void Reset();
  
  ClassDef(SFAttenuation,1)
};

#endif