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
#include "TObject.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "SFData.hh"
#include "SFTools.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

class SFPositionRes : public TObject{
    
private:
    int                  fSeriesNo;
    double               fPosResSeries;
    double               fPosResSeriesErr;
    SFData              *fData;
    TGraphErrors        *fPosVsMean;
    TGraphErrors        *fPosVsSigmaAtt;
    TGraphErrors        *fRecoPosition;
    std::vector <TH1D*>  fQRatios;
    std::vector <TH1D*>  fSpecCh0;
    std::vector <TH1D*>  fSpecCh1;
    std::vector <double> fPosRes;
    std::vector <double> fPosResErr;
    SFAttenuation       *fAtt;
    
    bool LoadRatios();
    
public:
    SFPositionRes(int seriesNo);
    ~SFPositionRes();
    
    bool          AnalyzePositionRes(void); 
    bool          ReconstructPos(void);
    TGraphErrors *GetMeanGraph(void);
    TGraphErrors *GetPositionResGraph(void);
    TGraphErrors *GetRecoPosGraph(void);
    std::vector <TH1D*>  GetRatios(void);
    std::vector <TH1D*>  GetSpectra(int ch);
    std::vector <double> GetPositionResSeries(void);
    std::vector <double> GetPositionRes(void);
    std::vector <double> GetPositionResError(void);
    
    void Print();
    
    
  ClassDef(SFPositionRes,1)
};

#endif
