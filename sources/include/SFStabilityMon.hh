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
#include "TObject.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "SFData.hh"
#include "SFTools.hh"
#include "SFPeakFinder.hh"
#include "SFDrawCommands.hh"
#include <vector>

class SFStabilityMon : public TObject{
    
private:
    int           fSeriesNo;
    SFData       *fData;
    TGraphErrors *fCh0Graph;
    TGraphErrors *fCh1Graph;
    TGraphErrors *fCh0ResGraph;
    TGraphErrors *fCh1ResGraph;
    double        fCh0StdDev;
    double        fCh1StdDev;
    double        fCh0Mean;
    double        fCh1Mean;
    
    std::vector <TH1D*> fSpecCh0;
    std::vector <TH1D*> fSpecCh1; 
    
public:
    SFStabilityMon(int seriesNo);
    ~SFStabilityMon();
    
    bool AnalyzeStability(int ch);
    
    TGraphErrors*       GetPeakPosGraph(int ch);
    TGraphErrors*       GetResidualsGraph(int ch);
    double              GetStdDev(int ch);
    double              GetMean(int ch);
    std::vector <TH1D*> GetSpectra(int ch);
    void                Print();
    
    ClassDef(SFStabilityMon, 1)
    
};

#endif
