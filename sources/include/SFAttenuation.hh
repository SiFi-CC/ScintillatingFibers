// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFAttenuation.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFAttenuation_H_
#define __SFAttenuation_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include <iostream>

using namespace std;

class SFAttenuation : public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  
  double fAttnLen;		///< Attenuation length determined with averaged channels method [mm]
  double fAttnErr;		///< Error on attenuation length fAttnLen [mm]
  vector <TH1D*> fRatios;	///< Vector containing histograms of signal ratios from both channels, for whole series 
  TGraphErrors *fAttnGraph;	///< Attenuation graph i.e. ln(M_{FB}) vs. source position
  
  double fAttnLenCh0;
  double fAttnLenCh1;
  double fAttnErrCh0;
  double fAttnErrCh1;
  TGraphErrors *fAttnGraphCh0;
  TGraphErrors *fAttnGraphCh1;
  vector <TH1D*> fSpectraCh0;
  vector <TH1D*> fSpectraCh1;
  vector <TH1D*> fPeaksCh0;
  vector <TH1D*> fPeaksCh1;
  
public:
  SFAttenuation();
  SFAttenuation(int seriesNo);
  ~SFAttenuation();
  
  bool            AttAveragedCh(void);
  vector <TH1D*>  GetRatios(void);
  vector <double> GetAttenuation(void);
  TGraphErrors*   GetAttGraph(void);
  double          GetAttLength(void);
  double          GetAttError(void);
  
  bool            AttSeparateCh(int ch);
  vector <TH1D*>  GetSpectra(int ch);
  vector <TH1D*>  GetPeaks(int ch);
  vector <double> GetAttenuation(int ch);
  TGraphErrors*   GetAttGraph(int ch);
  double          GetAttLength(int ch);
  double          GetAttError(int ch);
  
  void Print(void);
  void Clear();
  
  ClassDef(SFAttenuation,1)
};

#endif