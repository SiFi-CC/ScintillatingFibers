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

struct GlobalChi2
{
    GlobalChi2(ROOT::Math::IMultiGenFunction& f1, ROOT::Math::IMultiGenFunction& f2)
        : fChi2_1(&f1), fChi2_2(&f2)
    {
    }

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

    const ROOT::Math::IMultiGenFunction* fChi2_1;
    const ROOT::Math::IMultiGenFunction* fChi2_2;
};

class SFAttenuationModel : public TObject
{

  private:
    int     fSeriesNo;
    SFData* fData;

    TGraphErrors* fMAttCh0Graph;
    TGraphErrors* fMAttCh1Graph;
    TGraphErrors* fMAttCh0CorrGraph;
    TGraphErrors* fMAttCh1CorrGraph;

    TF2* fPlRecoFun;
    TF2* fPrRecoFun;
    
    SFResults*            fResults;
    ROOT::Fit::FitResult* fFitterResults;

  public:
    SFAttenuationModel(int seriesNo);
    ~SFAttenuationModel();

    bool FitModel(void);

    SFResults* GetResults(void) { return fResults; };

    void Print(void);

    ClassDef(SFAttenuationModel, 1)
};

#endif
