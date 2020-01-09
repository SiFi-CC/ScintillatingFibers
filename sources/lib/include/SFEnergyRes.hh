// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFEnergyRes.hh              *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************
#ifndef __SFEnergyRes_H_
#define __SFEnergyRes_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include <iostream>

struct EnergyResResults{
    
    double fEnergyResCh0    = -1;
    double fEnergyResCh0Err = -1;
    
    double fEnergyResCh1    = -1;
    double fEnergyResCh1Err = -1;
    
    double fEnergyResSum    = -1;
    double fEnergyResSumErr = -1;
};

class SFEnergyRes: public TObject{
 
private:
  int     fSeriesNo;         ///< Number of experimental series to be analyzed
  SFData *fData;             ///< SFData object of the analyzed series
  
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
                                        ///< collimator is analyzed [not used]
  std::vector <TH1D*> fPeaksCh1;
  std::vector <TH1D*> fPeaksSum;
  
  EnergyResResults fResults;
  
public:
  SFEnergyRes(int seriesNo);
  ~SFEnergyRes();
    
  bool CalculateEnergyRes(int ch);
  bool CalculateEnergyRes(void);
  
  EnergyResResults     GetResults(void) { return fResults; };
  TGraphErrors*        GetEnergyResolutionGraph(void);
  TGraphErrors*        GetEnergyResolutionGraph(int ch);
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetSpectraCorrected(int ch);
  std::vector <TH1D*>  GetSpectraSum(void); 
  std::vector <TH1D*>  GetPeaks(void);
  std::vector <TH1D*>  GetPeaks(int ch);
  
  void Print(void);
  
  ClassDef(SFEnergyRes,1)
};

#endif
