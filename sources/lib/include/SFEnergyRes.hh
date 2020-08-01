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

/// Structure containing numerical results of the 
/// energy resolution analysis.

struct SFEnergyResResults
{
    double fEnergyResCh0    = -1; ///< Energy resolution of channel 0
    double fEnergyResCh0Err = -1; ///< Uncertainty of energy resolution of channel 0
    
    double fEnergyResCh1    = -1; ///< Energy resolution of channel 1
    double fEnergyResCh1Err = -1; ///< Uncertainty of energy resolution of channel 1
    
    double fEnergyResAve    = -1; ///< Energy resolution for averaged channels
    double fEnergyResAveErr = -1; ///< Uncertainty of energy resolution for averaged channels
};

/// Class performing analysis of energy resolution for requested series. 
/// Energy resolution is determined separately for channels 0 and 1 as 
/// well as for the averaged spectra. 
/// 
/// Energy resolution for single measurement is calculated as follows:
/// \f[
/// ER = \frac{\sigma_{511}}{\mu_{511}} \cdot 100\%
/// \f]
/// where \f$\sigma_{511}\f$ and \f$\mu_{511}\f$ are parameters of the 511 keV
/// peak determined with the use of SFPeakFinder class. Uncertainty of \f$ER\f$ 
/// is calculated as follows:
/// \f[
/// \Delta ER = ER \cdot \sqrt{\frac{\Delta \mu_{511}^2}{\mu_{511}^2} + 
/// \frac{\Delta \sigma_{511}^2}{\sigma_{511}^2}} \cdot 100\%
/// \f]
/// Results of calculations are plotted on graphs: fEnergyResGraphCh0, 
/// fEnergyResCh1 and fEnergyResAve. Each of them is accessible via dedicated 
/// functions.
/// 
/// Final energy resolution values for the series are calculated as mean weighted 
/// with the uncertainties as follows:
/// \f[
/// \overline{ER} = \frac{\sum_{i} ER_i/\Delta ER_i^2}{\sum_i 1/\Delta ER_1^2}
/// \f]
/// And the uncertainty of the mean value is calculated as follows:
/// \f[
/// \Delta \overline{ER} = \sqrt{\frac{1}{\sum_i 1 / \Delta ER_i^2}}
/// \f]
/// Numerical results of the anaylis are stored in fResults structure (EnergyResResults 
/// type object) and are accesible via dedicated function. 

class SFEnergyRes: public TObject
{
 
private:
  int     fSeriesNo; ///< Number of experimental series to be analyzed
  SFData *fData;     ///< SFData object of the analyzed series
  
  TGraphErrors *fEnergyResGraphCh0;  ///< Energy resolution graph for channel 0
  TGraphErrors *fEnergyResGraphCh1;  ///< Energy resolution graph for channel 1
  TGraphErrors *fEnergyResGraphAve;  ///< Energy resolution graph for averaged spectra
  
  std::vector <TH1D*> fSpectraCh0;   ///< Vector containing charge spectra from channel 0
  std::vector <TH1D*> fSpectraCh1;   ///< Vector containing charge spectra from channel 1
  std::vector <TH1D*> fSpectraAve;   ///< Vector containing averaged charge spectra
  
  std::vector <TH1D*> fPeaksCh0;     ///< Vector containing spectra after background subtraction (currently not used)
  std::vector <TH1D*> fPeaksCh1;     ///< Vector containing spectra after background subtraction (currently not used)
  std::vector <TH1D*> fPeaksAve;     ///< Vector containing spectra after background subtraction (currently not used)
  
  SFEnergyResResults fResults;         ///< Structure containing numerical results od analysis
  
public:
  SFEnergyRes(int seriesNo);
  ~SFEnergyRes();
    
  bool CalculateEnergyRes(int ch);
  bool CalculateEnergyRes(void);
  
  /// Returns structure fResults numerical results of the energy resolution anaalysis 
  SFEnergyResResults   GetResults(void) { return fResults; };
  TGraphErrors*        GetEnergyResolutionGraph(void);
  TGraphErrors*        GetEnergyResolutionGraph(int ch);
  std::vector <TH1D*>  GetSpectra(int ch);
  std::vector <TH1D*>  GetSpectra(void); 
  std::vector <TH1D*>  GetPeaks(void);
  std::vector <TH1D*>  GetPeaks(int ch);
  
  void Print(void);
  
  ClassDef(SFEnergyRes,1)
};

#endif
