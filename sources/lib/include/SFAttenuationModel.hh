// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFAttenuationModel.hh         *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#ifndef __SFAttenuationModel_H_
#define __SFAttenuationModel_H_ 1

#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFResults.hh"

#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TMatrixT.h>
#include <TMatrixDfwd.h>

#include <cassert>

extern int iparSl[];
extern int iparSr[];

/// Structure for calculation of chi2. 
struct GlobalChi2Model
{
    /// Constructor.
    GlobalChi2Model(ROOT::Math::IMultiGenFunction& f1, ROOT::Math::IMultiGenFunction& f2)
        : fChi2_1(&f1), fChi2_2(&f2)
    {
    }

    /// Chi2 calculation.
    double operator()(const double* par) const
    {
        const int n = 6;
        double    p1[n];
        double    p2[n];

        for (int i = 0; i < n; ++i)
            p1[i] = par[iparSl[i]];

        for (int i = 0; i < n; ++i)
            p2[i] = par[iparSr[i]];

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }

    const ROOT::Math::IMultiGenFunction* fChi2_1; ///< chi2_1
    const ROOT::Math::IMultiGenFunction* fChi2_2; ///< chi2_2
};

/// This namespace fits exponential attenuation model with light reflection
/// to the experimental data. Based on the fit results and available
/// data primary light component is reconstructed. As one of the fitted 
/// parameters attenuation length is determined. Full uncertainties 
/// calculus is included.

namespace SFAttenuationModel
{
    double CalculateUncertainty(std::vector<double> params, 
                                TMatrixD covMatrix, TString side);
    SFResults* FitModel(TGraphErrors* left, TGraphErrors *right, double pos_uncert, double fiberLen);
};

#endif
