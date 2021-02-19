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

#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"

#include <TF1.h>
#include <TGraphErrors.h>
#include <TObject.h>

#include <iostream>

/// Class to determine attenuation length. This class is suitable only for experimental
/// series with different positions of source. Two methods of attenuation length
/// determination are available: AttAveragedCh() - based on Pauwels et al., JINST 8 (2013)
/// P09019, where combined signal from both channels is analyzed, and AttSeparateCh(),
/// where attenuation length is calulated for each channel separately.

class SFAttenuation : public TObject
{

  private:
    int     fSeriesNo; ///< Number of experimental series to be analyzed
    SFData* fData;     ///< SFData object of the analyzed series

    std::vector<TH1D*> fRatios;     ///< Vector containing histograms of ln(M_LR) distributions
    std::vector<TH1D*> fSpectraCh0; ///< Vector containing charge spectra from channel 0
    std::vector<TH1D*> fSpectraCh1; ///< Vector containing charche spectra from channel 1
    std::vector<TH1D*> fPeaksCh0;   ///< Vector containing 511 keV peaks, channel 0 
                                    ///< [not used at the moment]
    std::vector<TH1D*> fPeaksCh1;   ///< Vector containing 511 keV peaks, channel 1 
                                    ///< [not used at the moment]

    TGraphErrors* fAttGraph;    ///< Attenuation graph i.e. ln(M_LR) vs. source position
    TGraphErrors* fSigmaGraph;  ///< Graph sigma of ln(M_LR) vs. source position
    TGraphErrors* fAttCh0Graph; ///< Attenuation graph for channel 0
    TGraphErrors* fAttCh1Graph; ///< Attenuation graph for channel 1

    SFResults* fResultsCh0;      ///< Results of attenuation analysis for channel 0
    SFResults* fResultsCh1;      ///< Results of attenuation analysis for channel 1
    SFResults* fResultsCombPol1; ///< Results of combined channels analysis (fitting pol1)
    SFResults* fResultsCombPol3; ///< Results of combined channels analysis (fitting pol3)

  public:
    SFAttenuation(int seriesNo);
    ~SFAttenuation();

    bool AttCombinedCh(void);
    bool AttSeparateCh(int ch);
    bool Fit1stOrder(void);
    bool Fit3rdOrder(void);

    std::vector<TH1D*> GetSpectra(int ch);
    std::vector<TH1D*> GetPeaks(int ch);
    std::vector<TH1D*> GetRatios(void);

    std::vector<SFResults*> GetResults(void);

    void Print(void);

    ClassDef(SFAttenuation, 1)
};

#endif
