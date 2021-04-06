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

#include "FitterFactory.h"
#include "SFData.hh"
#include "SFResults.hh"
#include "SFTools.hh"

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TObject.h>
#include <TSpectrum.h>

#include <fstream>
#include <iostream>
#include <vector>

/// Class searching for the 511 keV peak and performing peak fitting and background
/// subtraction on measured charge spectra.
/// Function fitted to the analyzed spectra has the following form:
/// \f[
/// f(Q) = \mathrm{signal} + \mathrm{background}
/// \f]
/// where signal is described with gaussian function:
/// \f[
/// f_{sig}(Q) = c \cdot \exp{\left(\frac{1}{2} \left(\frac{Q-\mu}{\sigma}\right)^2\right)}
/// \f]
/// where \f$c\f$, \f$\mu\f$ and \f$\sigma\f$ are parameters desribing 511 kev peak.
/// And background is defined as:
/// \f[
/// f_{bg}(Q) = p_0 + p_1 \cdot e^{(Q-p_2) \cdot p_3}
/// \f]

class SFPeakFinder : public TObject
{

  private:
    TH1D*      fSpectrum;  ///< Analyzed experimental spectrum
    TH1D*      fPeak;      ///< Histogram with the chosen peak after background subtraction
    TF1*       fFittedFun; ///< Function fitted to the analyzed spectrum: Gauss+background
    int        fID;
    bool       fVerbose;   ///< Print-outs level
    bool       fTests;     ///< Flag for testing mode
    SFResults* fResults; ///< Object containing parameters of 511 keV peak as determined by the fit

  public:
    SFPeakFinder();
    SFPeakFinder(TH1D* spectrum, int ID, bool verbose, bool tests);
    SFPeakFinder(TH1D* spectrum, bool verbose, bool tests);
    SFPeakFinder(TH1D* spectrum, bool verbose);
    SFPeakFinder(TH1D* spectrum);
    ~SFPeakFinder();
    
    TString InitSpectrumPerPosition(int seriesNo);
    TString InitSpectrumPerSeries(int seriesNo);
    TString Init(void);
    bool    FindPeakRange(double& min, double& max);
    bool    FindPeakFit(void);
    bool    SubtractBackground(void);
    void    SetSpectrum(TH1D* spectrum);
    void    Print(void);

    /// Returns structure containing parameters of the 511 keV peak.
    SFResults* GetResults(void) { return fResults; };
    /// Sets print-outs level.
    void SetVerbLevel(bool verbose) { fVerbose = verbose; };
    /// Sets testing mode.
    void SetTests(bool tests) { fTests = tests; };

    ClassDef(SFPeakFinder, 1)
};

#endif
