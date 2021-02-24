// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionRes.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFPositionRes_H_
#define __SFPositionRes_H_ 1

#include "SDDSamples.h"
#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFTools.hh"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TObject.h>
#include <TTree.h>

#include <iostream>

class SFPositionRes : public TObject
{

  private:
    int            fSeriesNo;
    SFData*        fData;
    SFAttenuation* fAtt;

    TGraphErrors* fPosRecoVsPosGraph;
    TGraphErrors* fPosResVsPosGraph;
    TGraphErrors* fMLRvsPosGraph;
    TGraphErrors* fResidualGraph;

    std::vector<TH1D*> fQRatios;
    std::vector<TH1D*> fPosRecoDist;
    std::vector<TH1D*> fSpecAv;

    SFResults* fResults;

    bool LoadRatios(void);

  public:
    SFPositionRes(int seriesNo);
    ~SFPositionRes();

    bool AnalyzePositionRes(void);

    std::vector<TH1D*> GetRatios(void);
    std::vector<TH1D*> GetPositionRecoDist(void);
    std::vector<TH1D*> GetSpectra(void);

    SFResults* GetResults(void) { return fResults; };

    void Print();

    ClassDef(SFPositionRes, 1)
};

#endif
