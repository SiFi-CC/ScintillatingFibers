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
#include <iostream>

using namespace std;

/// Class to determine the energy resolution of the fibres for 511 keV. This class is suitable for experimental 
/// series with different positions of source, i.e. 1, 2, 3, 4 and 5. or for single positions of a series.

class SFEnergyResolution: public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  double fEnergyResCh0;		///< EnergyResolution Ch0 for all positions [%] 
  double fEnergyResErrCh0;	///< Error on fEnergyResCh0 [%]
  double fEnergyResCh1;		///< EnergyResolution Ch1 for all positions [%] 
  double fEnergyResErrCh1;	///< Error on fEnergyResCh1 [%]
  TGraphErrors *fEnergyResGraphCh0;	///< EnergyResolution graph for channel 0
  TGraphErrors *fEnergyResGraphCh1;	///< EnergyResolution graph for channel 1
  vector <TH1D*> fSpectraCh0;	///< Vector containing charge spectra from channel 0
  vector <TH1D*> fSpectraCh1;	///< Vector containing charge spectra from channel 1
  vector <SFPeakFinder*> fPFCh0;
  vector <SFPeakFinder*> fPFCh1;

  void CalculateER(int ch);

public:
  SFEnergyResolution();
  SFEnergyResolution(int seriesNo);
  ~SFEnergyResolution();
  
  vector <double> GetEnergyResolution(int ch); ///< Returns the fEnergyRes and fEnergyResErr in an vector for the requested channel
  TGraphErrors*   GetEnergyResolutionGraph(int ch);///< Returns the fEnergyResolutionGraph for the requested channel 

  vector <TH1D*>  GetSpectra(int ch);
  
  void Print(void);
  void Reset();
  
  ClassDef(SFEnergyResolution,1)
};

#endif
