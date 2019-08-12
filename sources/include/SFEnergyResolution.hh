// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFEnergyResolution.hh         *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************
#ifndef __SFEnergyResolution_H_
#define __SFEnergyResolution_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

class SFEnergyResolution: public TObject{
 
private:
  int     fSeriesNo;         ///< Number of experimental series to be analyzed
  SFData *fData;             ///< SFData object of the analyzed series
  double  fEnergyResCh0;     ///< Energy resolution of Ch0 for all positions [%] 
  double  fEnergyResCh0Err;  ///< Error on fEnergyResCh0 [%]
  double  fEnergyResCh1;     ///< Energy resolution of Ch1 for all positions [%] 
  double  fEnergyResCh1Err;  ///< Error on fEnergyResCh1 [%]
  double  fEnergyResSum;     ///< Energy resolution of summed spectra for all positions [%] 
  double  fEnergyResSumErr;  ///< Error on fEnergyResSum [%]
  
  TGraphErrors *fEnergyResGraphCh0;  ///< Energy resolution graph for channel 0
  TGraphErrors *fEnergyResGraphCh1;  ///< Energy resolution graph for channel 1
  TGraphErrors *fEnergyResGraphSum;  ///< Energy resolution graph for summed spectra
  
  std::vector <TH1D*> fSpectraCh0;      ///< Vector containing charge spectra from channel 0
  std::vector <TH1D*> fSpectraCh1;      ///< Vector containing charge spectra from channel 1
  std::vector <TH1D*> fSpectraCorrCh0;  ///< Vector containing charge spectra from channel 0
                                        ///< corrected for attenuationlength 
  std::vector <TH1D*> fSpectraCorrCh1;  ///< Vector containing charge spectra from channel 1 
                                        ///< corrected for attenuation length
  std::vector <TH1D*> fSpectraSum;      ///< Vector containing charge spectra summed and corrected 
                                        ///< for attenuation length
  std::vector <TH1D*> fPeaksCh0;        ///< Vector containing spectra after background subtraction;
                                        ///< this vector is filled only if series measured with lead
                                        ///< collimator is analyzed
  std::vector <TH1D*> fPeaksCh1;
  std::vector <TH1D*> fPeaksSum;
  
public:
  SFEnergyResolution();
  SFEnergyResolution(int seriesNo);
  ~SFEnergyResolution();
    
  bool CalculateEnergyRes(int ch);
  bool CalculateEnergyRes(void);
  
  std::vector <double> GetEnergyResolution(void);
  std::vector <double> GetEnergyResolution(int ch);
  TGraphErrors*        GetEnergyResolutionGraph();
  TGraphErrors*        GetEnergyResolutionGraph(int ch);
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetSpectraCorrected(int ch);
  std::vector <TH1D*>  GetSpectraSum(void); 
  std::vector <TH1D*>  GetPeaks(void);
  std::vector <TH1D*>  GetPeaks(int ch);
  
  void Print(void);
  
  ClassDef(SFEnergyResolution,1)
};

#endif
