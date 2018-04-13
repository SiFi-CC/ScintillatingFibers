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
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include <iostream>

using namespace std;

/// Class to determine timing resolution.

class SFTimingRes : public TObject{
 
private:
  int     fSeriesNo;		///< Number of analyzed series.
  TString fThreshold;		///< Falg to identify which data should be analyzed.
  TString fMethod;		///< Flag to determine type of analysis - with or without energy cut.
  SFData  *fData;		///< Experimental series to be analyzed.
  
  vector <TH1D*>  fT0Diff;
  vector <double> fTimeRes;
  vector <double> fTimeResErr;
  
  int GetIndex(double position);
  bool AnalyzeWithECut(void);
  bool AnalyzeNoECut(void);
  
public:
  SFTimingRes();
  SFTimingRes(int seriesNo, TString threshold, TString method);
  ~SFTimingRes();

  vector <TH1D*>  GetT0Diff(void);
  TH1D*           GetT0Diff(double position);
  vector <double> GetTimingResolution(double position);
  vector <double> GetTimingResolutions(void);
  vector <double> GetTimingResErrors(void);
  void            Print(void);
  
  ClassDef(SFTimingRes,1)
};

#endif