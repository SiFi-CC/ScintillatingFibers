// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimingRes.hh             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFTimingRes_H_
#define __SFTimingRes_H_ 1
#include "TObject.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFTools.hh"
#include <iostream>

/// Structure containig numerical results of the
/// timing resolution analysis.

struct SFTimingResResults{
    
    double fTimeRes    = -1;
    double fTimeResErr = -1;
    
    double fTimeResECut    = -1;
    double fTimeResECutErr = -1;
    
    std::vector <double> fTimeResAll;
    std::vector <double> fTimeResAllErr;
    
    std::vector <double> fTimeResECutAll;
    std::vector <double> fTimeResECutAllErr;
};

/// Class to determine timing resolution. Two methods available: "no cut" - cut imposed only
/// on scattered events, "with cut" - cut imposed on scattered events and additional energy 
/// cut on 511 keV peak. Results are returned as vectors containig numerical timing resolutions
/// and in a form of graph: mean of ch_0.T0-ch_1.T0 distribution vs. source position. 
/// IMPORTANT: timing resolution is determined as FWHM of lorentzian and gaussian functions
/// fitted to the ch_0.fT0-ch_1.fT0 distributions.

class SFTimingRes : public TObject{
 
private:
  int     fSeriesNo;   ///< Number of analyzed series.
  SFData *fData;       ///< Experimental series to be analyzed.

  std::vector <TH1D*>  fRatios;
  std::vector <TH1D*>  fSpecCh0;
  std::vector <TH1D*>  fSpecCh1;
  std::vector <TH1D*>  fT0Diff; 
  std::vector <TH1D*>  fT0DiffECut;
  
  TGraphErrors *fTResGraph;        
  TGraphErrors *fTResECutGraph;
  
  SFTimingResResults fResults;

  bool LoadRatios(void);
  
public:
  SFTimingRes(int seriesNo);
  ~SFTimingRes();
  
  bool AnalyzeWithECut(void);
  bool AnalyzeNoECut(void);
  void Print(void);
  
  std::vector <TH1D*>  GetRatios(void);
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetT0Diff(bool type);
  TGraphErrors*        GetTimingResGraph(bool type);
  SFTimingResResults   GetResults(void) { return fResults; };

  ClassDef(SFTimingRes,1)
};

#endif
