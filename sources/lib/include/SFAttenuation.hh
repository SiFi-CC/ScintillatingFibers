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

#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFTools.hh"
#include "SFInfo.hh"

#include <TF1.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>

#include <iostream>

class TGraphErrors;

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

namespace SFAttenuation
{
    SFResults* AttCombinedCh(TString fun, SFInfo* info, 
                             double pos_uncert,
                             std::vector<double> positions,
                             std::vector<TH1D*> spectra);
    
    SFResults* AttSeparateCh(char side, double pos_uncert, SFInfo *info,
                             std::vector<double> positions,
                             std::vector<TH1D*> spectra,
                             std::vector<TString> path);
    
    SFResults* AttSimultaneousFit(TGraphErrors *attL, TGraphErrors *attR);
};

#endif
