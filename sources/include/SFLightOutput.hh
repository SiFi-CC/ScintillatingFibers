// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFLightOutput.cc           *
// *         	Jonas Kasper		   *
// *   	  kasper@physik.rwth-aachen.de	   *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFLightOutput_H_
#define __SFLightOutput_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

using namespace std;

/// Class to determine the light output of the fibres. This class is suitable for experimental 
/// series with different positions of source, i.e. 1, 2, 3, 4 and 5. or for single positions of a series.
/// The amount of recorded ligth is corrected for the attenuation in the fibre and the sipm pde.

class SFLightOutput : public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  SFAttenuation *fAtt;		///< SFAttenuation object of the analyzed series
  double fPDE; ///< Photon Detection Efficiency of the SiPM
  double fLightOutAve;		///< Lightoutput averaged for all positions [#photons/MeV] 
  double fLightOutAveErr;		///< Error on lightoutput fLightOut [#photons/MeV]
  double fLightOutAveCh0;		///< Lightoutput Ch0 for all positions [#photons/MeV] 
  double fLightOutAveErrCh0;		///< Error on lightoutput fLightOutAveCh0 [#photons/MeV]
  double fLightOutAveCh1;		///< Lightoutput Ch1 for all positions [#photons/MeV] 
  double fLightOutAveErrCh1;		///< Error on lightoutput fLightOutAveCh1 [#photons/MeV]
  TGraphErrors *fLightOutAveGraph;	///< LightOutput graph i.e. LO vs. source position
  TGraphErrors *fLightOutAveGraphCh0;	///< LightOutput graph for channel 0
  TGraphErrors *fLightOutAveGraphCh1;	///< LightOutput graph for channel 1
  double fLightOutSep;		///< Lightoutput averaged for all positions [#photons/MeV] - cal wiht sep attlen 
  double fLightOutSepErr;		///< Error on lightoutput fLightOut [#photons/MeV] - cal wiht sep attlen
  double fLightOutSepCh0;		///< Lightoutput Ch0 for all positions [#photons/MeV]  - cal wiht sep attlen
  double fLightOutSepErrCh0;		///< Error on lightoutput fLightOutAveCh0 [#photons/MeV] - cal wiht sep attlen
  double fLightOutSepCh1;		///< Lightoutput Ch1 for all positions [#photons/MeV]  - cal wiht sep attlen
  double fLightOutSepErrCh1;		///< Error on lightoutput fLightOutAveCh1 [#photons/MeV] - cal wiht sep attlen
  TGraphErrors *fLightOutSepGraph;	///< LightOutput graph i.e. LO vs. source position - cal wiht sep attlen
  TGraphErrors *fLightOutSepGraphCh0;	///< LightOutput graph for channel 0 - cal wiht sep attlen
  TGraphErrors *fLightOutSepGraphCh1;	///< LightOutput graph for channel 1 - cal wiht sep attlen
  vector <TH1D*> fSpectraCh0;	///< Vector containing charge spectra from channel 0
  vector <TH1D*> fSpectraCh1;	///< Vector containing charge spectra from channel 1
  vector <TH1D*> fPeaksCh0;	///< Vector containing 511 keV peaks, channel 0
  vector <TH1D*> fPeaksCh1;	///< Vector containing 511 keV peaks, channel 1
  vector <double> fAttLen; ///< Contains the averaged Attenuation Length and the Error
  vector <double> fAttLenCh0; ///< Contains the Attenuation Length and the Error for Ch0
  vector <double> fAttLenCh1; ///< Contains the Attenuation Length and the Error for Ch1

  void CalculateLO(bool mode,int ch);
  void CalculateSLO(bool mode);

public:
  SFLightOutput();
  SFLightOutput(int seriesNo);
  ~SFLightOutput();
  
  vector <double> GetLightOutput(TString mode); ///< Returns the fLightOut and fLightOutErr in an vector
  TGraphErrors*   GetLightOutputGraph(TString mode);///< Returns the fLightOutGraph 

  vector <TH1D*>  GetSpectra(int ch);
  vector <TH1D*>  GetPeaks(int ch);
  vector <double> GetLightOutput(TString mode,int ch); ///< Returns the fLightOut and fLightOutErr in an vector
  TGraphErrors*   GetLightOutputGraph(TString mode,int ch);///< Returns the fLightOutGraph corresponding to one Channel
  
  void Print(void);
  void Reset();
  
  ClassDef(SFLightOutput,1)
};

#endif
