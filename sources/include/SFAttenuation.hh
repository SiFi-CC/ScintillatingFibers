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

/// Class to determine attenuation length. This class is suitable only for experimental 
/// series with different positions of source. Two methods of attenuation length
/// determination are available: AttAveragedCh() - based on Pauwels et al., JINST 8 (2013)
/// P09019, where combined signal from both channels is analyzed, and AttSeparateCh()
/// - where attenuation length is calulated for each channel separately.

class SFAttenuation : public TObject{
 
private:
  int                 fSeriesNo;       ///< Number of experimental series to be analyzed
  SFData*             fData;           ///< SFData object of the analyzed series
  double              fAttnLenPol1;    ///< Attenuation length determined with averaged channels method [mm]
  double              fAttnLenPol1Err; ///< Error on attenuation length fAttnLen [mm]
  double              fAttnLenPol3;
  double              fAttnLenPol3Err;
  std::vector <TH1D*> fRatios;          ///< Vector containing histograms of signal ratios from both channels,
                                    ///< for whole series 
  TGraphErrors* fAttnGraph;         ///< Attenuation graph i.e. ln(M_{FB}) vs. source position
  TGraphErrors* fSigmaGraph;
  double        fAttnLenCh0;        ///< Attenuation length for channel 0
  double        fAttnLenCh1;        ///< Attenuation length for channel 1
  double        fAttnErrCh0;        ///< Error on attenuation length for channel 0
  double        fAttnErrCh1;        ///< Error on attenuation length for channel 1
  TGraphErrors *fAttnGraphCh0;      ///< Attenuation graph for channel 0
  TGraphErrors *fAttnGraphCh1;      ///< Attenuation graph for channel 1
  std::vector <TH1D*> fSpectraCh0;  ///< Vector containing charge spectra from channel 0
  std::vector <TH1D*> fSpectraCh1;  ///< Vector containing charche spectra from channel 1
  std::vector <TH1D*> fPeaksCh0;    ///< Vector containing 511 keV peaks, channel 0 [not used at the moment]
  std::vector <TH1D*> fPeaksCh1;    ///< Vector containing 511 keV peaks, channel 1 [not used at the moment]
  
  public:
  SFAttenuation(int seriesNo);
  ~SFAttenuation();
  
  bool                 AttAveragedCh(void);
  bool                 Fit3rdOrder(void);
  std::vector <TH1D*>  GetRatios(void);
  std::vector <double> GetAttLenPol1(void);
  std::vector <double> GetAttLenPol3(void);
  TGraphErrors*        GetAttGraph(void);
  TGraphErrors*        GetSigmaGraph(void);
  
  bool                 AttSeparateCh(int ch);
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetPeaks(int ch);
  std::vector <double> GetAttenuation(int ch);
  TGraphErrors*        GetAttGraph(int ch);
  
  void Print(void);
  
  ClassDef(SFAttenuation,1)
};

#endif
