// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFEnergyResolution.cc         *
// *         	Jonas Kasper		   *
// *   	  kasper@physik.rwth-aachen.de	   *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFEnRes_H_
#define __SFEnRes_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

using namespace std;

/// Class to determine the energy resolution of the fibres for 511 keV. This class is suitable for experimental 
/// series with different positions of source, i.e. 1, 2, 3, 4 and 5. or for single positions of a series.

class SFEnergyResolution: public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  double fEnergyResAve;		///< EnergyResolution Averaged for all positions [%] 
  double fEnergyResErrAve;	///< Error on fEnergyResAve [%]
  double fEnergyResCh0;		///< EnergyResolution Ch0 for all positions [%] 
  double fEnergyResErrCh0;	///< Error on fEnergyResCh0 [%]
  double fEnergyResCh1;		///< EnergyResolution Ch1 for all positions [%] 
  double fEnergyResErrCh1;	///< Error on fEnergyResCh1 [%]
  TGraphErrors *fEnergyResGraphAve;	///< EnergyResolution graph for channel 0
  TGraphErrors *fEnergyResGraphCh0;	///< EnergyResolution graph for channel 0
  TGraphErrors *fEnergyResGraphCh1;	///< EnergyResolution graph for channel 1
  vector <TH1D*> fSpectraCh0;	///< Vector containing charge spectra from channel 0
  vector <TH1D*> fSpectraCorCh0;	///< Vector containing charge spectra from channel 0 corrected for attenuationlength
  vector <TH1D*> fSpectraCh1;	///< Vector containing charge spectra from channel 1 
  vector <TH1D*> fSpectraCorCh1;	///< Vector containing charge spectra from channel 1 corrected for attenuation length
  vector <TH1D*> fSpectraAve;	///< Vector containing charge spectra sum corrected for attenuation
  vector <SFPeakFinder*> fPFCh0;
  vector <SFPeakFinder*> fPFCh1;

  void CalculateER(int ch);
  void CalculateERAve();

  SFAttenuation *fAtt;		///< SFAttenuation object of the analyzed series
  vector <double> fAttLen; ///< Contains the averaged Attenuation Length and the Error

public:
  SFEnergyResolution();
  SFEnergyResolution(int seriesNo);
  ~SFEnergyResolution();
  
  vector <double> GetEnergyResolution(); ///< Returns the fEnergyResAve and fEnergyResErrAve in an vector 
  vector <double> GetEnergyResolution(int ch); ///< Returns the fEnergyRes and fEnergyResErr in an vector for the requested channel
  TGraphErrors*   GetEnergyResolutionGraph();///< Returns the fEnergyResolutionGraphAve 
  TGraphErrors*   GetEnergyResolutionGraph(int ch);///< Returns the fEnergyResolutionGraph for the requested channel 

  vector <TH1D*>  GetSpectra(int ch);
  vector <TH1D*>  GetAveSpectra();
  vector <TH1D*>  GetAveSpectra(int ch);
  
  void Print(void);
  void Reset();
  
  ClassDef(SFEnergyResolution,1)
};

#endif
