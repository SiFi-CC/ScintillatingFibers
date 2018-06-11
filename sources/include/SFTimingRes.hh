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
#include <iostream>

using namespace std;

/// Class to determine timing resolution. Two methods available: "no cut" - cut imposed only
/// on scattered events, "with cut" - cut imposed on scattered events and additional energy 
/// cut on 511 keV peak. Results are returned as vectors containig numerical timing resolutions
/// and in a form of graph: mean of ch_0.T0-ch_1.T0 distribution vs. source position. 

class SFTimingRes : public TObject{
 
private:
  int          fSeriesNo;	///< Number of analyzed series.
  TString      fThreshold;	///< Falg to identify which data should be analyzed.
  TString      fMethod;		///< Flag to determine type of analysis - with or without energy cut.
  SFData       *fData;		///< Experimental series to be analyzed.
  TGraphErrors *fT0Graph;	///< Graph ch_0.T0-ch_1.T0 vs. source position
  
  vector <TH1D*>  fRatios;	///< Ratio histograms necessary to impose cuts
  vector <TH1D*>  fPEch0;	///< PE spectra from ch0, necessary to impose 511 keV energy cut
  vector <TH1D*>  fPEch1;	///< PE spectra from ch1, necessary to impose 511 keV energy cut
  vector <TH1D*>  fT0Diff;	///< Histograms ch_0.T0-ch_1.T0 
  vector <double> fTimeRes;	///< Vector containing timing resolution values for whole series
  vector <double> fTimeResErr;	///< Vector containing errors on timing resolutions for whole series
  
  int GetIndex(double position);
  bool LoadRatios(void);
  bool AnalyzeWithECut(void);
  bool AnalyzeNoECut(void);
  
public:
  SFTimingRes();
  SFTimingRes(int seriesNo, TString threshold, TString method);
  ~SFTimingRes();

  vector <TH1D*>  GetT0Diff(void);
  TH1D*           GetT0Diff(double position);
  TGraphErrors*   GetT0Graph(void);
  vector <double> GetTimingResolution(double position);
  vector <double> GetTimingResolutions(void);
  vector <double> GetTimingResErrors(void);
  vector <TH1D*>  GetRatios(void);
  vector <TH1D*>  GetSpectra(int ch);
  void            Print(void);
  
  ClassDef(SFTimingRes,1)
};

#endif