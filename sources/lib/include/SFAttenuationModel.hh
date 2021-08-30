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

/// This class fits exponential attenuation model with light reflection
/// to the experimental data. Based on the fit results and available
/// data primary light component is reconstructed. As one of the fitted 
/// parameters attenuation length is determined. Full uncertainties 
/// calculus is included.

class SFAttenuationModel : public TObject
{

  private:
    int     fSeriesNo; ///< Number of experimental series
    SFData* fData;     ///< SFData object of the experimental series

    TGraphErrors* fMAttCh0Graph;     ///< Experimental attenuation curve ch0
    TGraphErrors* fMAttCh1Graph;     ///< Experimental attenuation curve ch1
    TGraphErrors* fMAttCh0CorrGraph; ///< Corrected attenuation curve ch0
    TGraphErrors* fMAttCh1CorrGraph; ///< Corrected attenuation curve ch1

    TF2* fPlRecoFun;     ///< Reconstructed primary component (left/ch0)
    TF2* fPrRecoFun;     ///< Reconstructed primary component (right/ch1)
    
    TMatrixD fCovMatrix; ///< Covariation matrix
    
    SFResults*            fResults;       ///< Analysis results
    ROOT::Fit::FitResult* fFitterResults; ///< Fitting results

  public:
    SFAttenuationModel(int seriesNo);
    ~SFAttenuationModel();

    double CalculateUncertainty(std::vector<double> params, TString side);
    bool   FitModel(void);

    /// Returns results of the analysis.
    SFResults* GetResults(void) { return fResults; };

    void Print(void);

    ClassDef(SFAttenuationModel, 1)
};

#endif
