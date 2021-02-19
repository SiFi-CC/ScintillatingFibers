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

#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "SFResults.hh"

#include <TF1.h>
#include <TGraphErrors.h>
#include <TObject.h>

#include <iostream>

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

class SFEnergyRes : public TObject
{

  private:
    int     fSeriesNo; ///< Number of experimental series to be analyzed
    SFData* fData;     ///< SFData object of the analyzed series

    TGraphErrors* fEnergyResCh0Graph; ///< Energy resolution graph for channel 0
    TGraphErrors* fEnergyResCh1Graph; ///< Energy resolution graph for channel 1
    TGraphErrors* fEnergyResAveGraph; ///< Energy resolution graph for averaged spectra

    std::vector<TH1D*> fSpectraCh0; ///< Vector containing charge spectra from channel 0
    std::vector<TH1D*> fSpectraCh1; ///< Vector containing charge spectra from channel 1
    std::vector<TH1D*> fSpectraAve; ///< Vector containing averaged charge spectra

    std::vector<TH1D*> fPeaksCh0; ///< Vector containing spectra after background subtraction (currently not used)
    std::vector<TH1D*> fPeaksCh1; ///< Vector containing spectra after background subtraction (currently not used)
    std::vector<TH1D*> fPeaksAve; ///< Vector containing spectra after background subtraction (currently not used)

    SFResults* fResultsCh0; ///< Object containing numerical results of analysis for channel 0
    SFResults* fResultsCh1; ///< Object containing numerical results of analysis for channel 0
    SFResults* fResultsAve; ///< Object containing numerical results of analysis for channel 0

  public:
    SFEnergyRes(int seriesNo);
    ~SFEnergyRes();

    bool CalculateEnergyRes(int ch);
    bool CalculateEnergyRes(void);

    std::vector<SFResults*> GetResults(void);
    std::vector<TH1D*>      GetSpectra(int ch);
    std::vector<TH1D*>      GetSpectra(void);
    std::vector<TH1D*>      GetPeaks(void);
    std::vector<TH1D*>      GetPeaks(int ch);

    void Print(void);

    ClassDef(SFEnergyRes, 1)
};

#endif
