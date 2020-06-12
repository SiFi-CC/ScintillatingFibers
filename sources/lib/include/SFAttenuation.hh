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

/// Structure containing numerical results of attenuation
/// length analysis.

struct AttenuationResults{
    
    double fAttCombPol1    = -1;   ///< Attenuation length determined with combined channels method and 1st degree polynomial fit
    double fAttCombPol1Err = -1;   ///< Uncertainty of fAttCombPol1
    
    double fAttCombPol3    = -1;   ///< Attenuation length determined with combined channels method and 3rd degree polynomial fit
    double fAttCombPol3Err = -1;   ///< Uncertainty of fAttCombPol3
    
    double fAttCh0    = -1;   ///< Attenuation length of channel 0
    double fAttCh0Err = -1;   ///< Uncertainty of fAttCh0
    
    double fAttCh1    = -1;   ///< Attenuation length of channel 1
    double fAttCh1Err = -1;   ///< Uncertainty of fAttCh1
};

/// Class to determine attenuation length. This class is suitable only for experimental 
/// series with different positions of source. Two methods of attenuation length
/// determination are available: AttAveragedCh() - based on Pauwels et al., JINST 8 (2013)
/// P09019, where combined signal from both channels is analyzed, and AttSeparateCh(),
/// where attenuation length is calulated for each channel separately.

class SFAttenuation : public TObject{
 
private:
  int                 fSeriesNo;       ///< Number of experimental series to be analyzed
  SFData*             fData;           ///< SFData object of the analyzed series
  
  std::vector <TH1D*> fRatios;      ///< Vector containing histograms of ln(M_LR) distributions
  std::vector <TH1D*> fSpectraCh0;  ///< Vector containing charge spectra from channel 0
  std::vector <TH1D*> fSpectraCh1;  ///< Vector containing charche spectra from channel 1
  std::vector <TH1D*> fPeaksCh0;    ///< Vector containing 511 keV peaks, channel 0 [not used at the moment]
  std::vector <TH1D*> fPeaksCh1;    ///< Vector containing 511 keV peaks, channel 1 [not used at the moment]
  
  TGraphErrors *fAttnGraph;         ///< Attenuation graph i.e. ln(M_LR) vs. source position
  TGraphErrors *fSigmaGraph;        ///< Graph sigma of ln(M_LR) vs. source position
  TGraphErrors *fAttnGraphCh0;      ///< Attenuation graph for channel 0
  TGraphErrors *fAttnGraphCh1;      ///< Attenuation graph for channel 1
  
  AttenuationResults fResults;      ///< Results of attenuation analysis
  
  public:
  SFAttenuation(int seriesNo);
  ~SFAttenuation();
  
  bool                 AttAveragedCh(void);
  bool                 AttSeparateCh(int ch);
  bool                 Fit1stOrder(void);
  bool                 Fit3rdOrder(void);
  
  TGraphErrors*        GetAttGraph(void);
  TGraphErrors*        GetAttGraph(int ch);
  TGraphErrors*        GetSigmaGraph(void);
  
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetPeaks(int ch);
  std::vector <TH1D*>  GetRatios(void);
  
  AttenuationResults   GetResults(void) { return fResults; };
  
  void Print(void);
  
  ClassDef(SFAttenuation,1)
};

#endif
