// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFStabilityMon.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFStabilityMon_H_
#define __SFStabilityMon_H_ 1

#include "SFData.hh"
#include "SFDrawCommands.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"
#include "SFTools.hh"

#include <TF1.h>
#include <TGraphErrors.h>
#include <TObject.h>
#include <TString.h>

#include <vector>

class SFStabilityMon : public TObject
{

  private:
    int           fSeriesNo;
    SFData*       fData;
    TGraphErrors* fCh0PeakPosGraph;
    TGraphErrors* fCh1PeakPosGraph;
    TGraphErrors* fCh0ResidualGraph;
    TGraphErrors* fCh1ResidualGraph;

    SFResults* fResultsCh0;
    SFResults* fResultsCh1;

    std::vector<TH1D*> fSpecCh0;
    std::vector<TH1D*> fSpecCh1;

  public:
    SFStabilityMon(int seriesNo);
    ~SFStabilityMon();

    bool AnalyzeStability(int ch);

    std::vector<SFResults*> GetResults(void);
    std::vector<TH1D*>      GetSpectra(int ch);
    void                    Print();

    ClassDef(SFStabilityMon, 1)
};

#endif
