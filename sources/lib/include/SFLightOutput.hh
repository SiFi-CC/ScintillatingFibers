// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFLightOutput.hh           *
// *             Jonas Kasper              *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFLightOutput_H_
#define __SFLightOutput_H_ 1

#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFTools.hh"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TObject.h>

#include <iostream>

class SFLightOutput : public TObject
{

  private:
    int    fSeriesNo;
    double fPDE;
    double fCrossTalk;

    TGraphErrors* fLightOutGraph;
    TGraphErrors* fLightOutCh0Graph;
    TGraphErrors* fLightOutCh1Graph;

    TGraphErrors* fLightColGraph;
    TGraphErrors* fLightColCh0Graph;
    TGraphErrors* fLightColCh1Graph;

    TCanvas* fInputData;

    std::vector<TH1D*>         fSpectraCh0;
    std::vector<TH1D*>         fSpectraCh1;
    std::vector<SFPeakFinder*> fPFCh0;
    std::vector<SFPeakFinder*> fPFCh1;

    SFData*        fData;
    SFAttenuation* fAtt;

    SFResults* fResultsLightOutCh0;
    SFResults* fResultsLightOutCh1;
    SFResults* fResultsLightOutSum;
    
    SFResults* fResultsLightColCh0;
    SFResults* fResultsLightColCh1;
    SFResults* fResultsLightColSum;

  public:
    SFLightOutput(int seriesNo);
    ~SFLightOutput();

    bool CalculateLightOut(void);
    bool CalculateLightOut(int ch);

    bool CalculateLightCol(void);
    bool CalculateLightCol(int ch);

    double GetCrossTalk(void);
    double GetPDE(void);

    std::vector<SFResults*> GetLOResults(void);
    std::vector<SFResults*> GetLCResults(void);
    
    TCanvas*   GetInputData(void) { return fInputData; };

    std::vector<TH1D*> GetSpectra(int ch);

    void Print(void);

    ClassDef(SFLightOutput, 1)
};

#endif
