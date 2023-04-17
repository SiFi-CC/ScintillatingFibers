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
#include "SFResults.hh"

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
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

namespace SFPeakFinder
{

    SFResults* FindPeakFit(TH1D* spectrum, TString path,
                           bool verbose = 0, bool tests = 0);
    auto       FindPeakRange(TH1D* spectrum, TString path, TString colimator,
                             bool verbose = 0, bool tests = 0)
                             -> std::tuple<double, double>;
    SFResults* SubtractBackground(TH1D* spectrum, TString path, TString colimator, 
                                  bool verbose = 0, bool tests = 0);

};

#endif
