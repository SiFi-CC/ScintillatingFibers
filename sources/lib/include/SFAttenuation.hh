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
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>

#include <iostream>

extern int iparSlAtt[];
extern int iparSrAtt[];

/// Structure for calculation of chi2. It is necessary for
/// the simultaneous fit of the simple exponanetial model. 
struct GlobalChi2Att
{
    /// Constructor.
    GlobalChi2Att(ROOT::Math::IMultiGenFunction& f1, ROOT::Math::IMultiGenFunction& f2)
        : fChi2_1(&f1), fChi2_2(&f2)
    {
    }

    /// Chi2 calculation.
    double operator()(const double* par) const
    {
        const int n = 4;
        double    p1[n];
        double    p2[n];

        for (int i = 0; i < n; ++i)
            p1[i] = par[iparSlAtt[i]];

        for (int i = 0; i < n; ++i)
            p2[i] = par[iparSrAtt[i]];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }

    const ROOT::Math::IMultiGenFunction* fChi2_1; ///< chi2_1
    const ROOT::Math::IMultiGenFunction* fChi2_2; ///< chi2_2
};

/// Class to determine attenuation length. This class is suitable only for experimental
/// series with different positions of source. Three methods of attenuation length
/// determination are available: AttCombinedCh() - based on Pauwels et al., JINST 8 (2013)
/// P09019, where combined signal from both channels is analyzed, AttSeparateCh(),
/// where attenuation length is calulated for each channel separately and FitSimultaneously(),
/// where exponential attenuation model is fitted simultaneously to both channels.

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

    TGraphErrors* fAttGraph; ///< Attenuation graph i.e. ln(M_LR) vs. source position
    TGraphErrors* fAttCh0;   ///< Attenuation graph for ch0 
    TGraphErrors* fAttCh1;   ///< Attenuation graph for ch1

    SFResults* fResultsCh0;      ///< Results of attenuation analysis for channel 0
    SFResults* fResultsCh1;      ///< Results of attenuation analysis for channel 1
    SFResults* fResultsCombPol1; ///< Results of combined channels analysis (fitting pol1)
    SFResults* fResultsCombPol3; ///< Results of combined channels analysis (fitting pol3)
    SFResults* fResultsExpSim;   ///< Results of simultaneous fitting of exponential model

  public:
    SFAttenuation(int seriesNo);
    ~SFAttenuation();

    bool AttCombinedCh(void);
    bool AttSeparateCh(int ch);
    bool Fit1stOrder(void);
    bool Fit3rdOrder(void);
    bool FitSimultaneously(void);

    std::vector<TH1D*> GetSpectra(int ch);
    std::vector<TH1D*> GetPeaks(int ch);
    std::vector<TH1D*> GetRatios(void);

    std::vector<SFResults*> GetResults(void);

    void Print(void);

    ClassDef(SFAttenuation, 1)
};

#endif
