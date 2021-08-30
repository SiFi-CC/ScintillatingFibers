// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFDrawCommands.hh           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFDrawCommands_H_
#define __SFDrawCommands_H_ 1
#include <TObject.h>
#include <TString.h>
#include <iostream>

/// \file
/// Enumeration representing different types of selections
/// for the analyzed data:

/// Enumeration representing different types of selections
/// for the analyzed data:
enum class SFSelectionType
{
    ///----- selections for single fiber measurements (Desktop Digitizer)
    kPE,                   ///< shows PE (calibrated charge) spectrum
    kCharge,               ///< shows uncalibrated charge spectrum
    kAmplitude,            ///< shows amplitude spectrum
    kT0,                   ///< shows T0 spectrum
    kTOT,                  ///< shows TOT spectrum
    kLogSqrtPERatio,       ///< shows PE ratio spectrum:
                           ///< \f$ \ln\sqrt{Q_{ch1}/Q_{ch0}} \f$
    kT0Difference,         ///< shows T0 difference spectrum:
                           ///< \f$ T_{0 ch0} - T_{0 ch1} \f$
    kPEAverage,            ///< shows spectrum of PE goemetric mean:
                           ///< \f$ \sqrt{Q_{ch0} \cdot Q_{ch1}} \f$
    kAmplitudeAverage,     ///< shows spectrum of amplitude geometric mean:
                           ///< \f$ \sqrt{A_{ch0} \cdot A_{ch1}} \f$
    kPECorrelation,        ///< shows 2D PE correlation spectrum:
                           ///< \f$Q_{ch0}\f$ vs. \f$Q_{ch1}\f$
    kAmplitudeCorrelation, ///< shows 2D amiplitude correlation spectrum:
                           ///< \f$A_{ch0}\f$ vs. \f$A_{ch1}\f$
    kT0Correlation,        ///< shows 2D T0 correlation spectrum:
                           ///< \f$T_{0 ch0}\f$ vs. \f$T_{0 ch1}\f$
    kAmpPECorrelation,     ///< shows 2D distribution for chosen channel i:
                           ///< \f$A_{chi}\f$ vs. \f$Q_{chi}\f$
    kPEAttCorrected,       ///< shows PE spectrum corrected for attenuation length:
                           ///< \f$ Q_{chi} \cdot /\exp{(-z/\lambda_{att})}\f$
    kPEAttCorrectedSum,    ///< shows sum of Ch0 and Ch1 PE spectra, both corrected
                           ///< for attenuation length: \f$ Q_{ch0}/\exp{\frac{z}
                           ///< {\lambda_{att}}} + Q_{ch1}/\exp{\frac{(L-z)}
                           ///< {\lambda_{att}}}\f$
    kBL,                   ///< shows base line histogram
    kBLSigma,              ///< shows histogram of base line sigma
    
    //----- selections for PMI measurements
    kPMICharge,               ///< shows charge (photon count) spectrum (PMI data)
    kPMIChargeAverage,        ///< shows charge (photon count) spectrum calculated
                              ///< as geometric mean of two fiber ends: 
                              ///< \f$ \sqrt{Q_{ch0} \cdot Q_{ch1}} \f$ (PMI data)
    kPMIT0,                   ///< shows T0 spectrum (PMI data)
    kPMIChargeCorrelation,    ///< shows charge (photon count) correlation spectum:
                              ///< \f$Q_{ch0}\f$ vs. \f$Q_{ch1}\f$ (PMI data)
    kPMIT0Correlation,        ///< shows time T0 correlation spectrum: 
                              ///< \f$T_{0 ch0}\f$ vs. \f$T_{0 ch1}\f$ (PMI data)
    kPMIT0Difference,         ///< shows time T0 difference spectrum: 
                              ///< \f$ T_{0 ch0} - T_{0 ch1} \f$ (PMI data)
    kPMILogSqrtChargeRatio,   ///< shows charge (photon count) ratio spectrum:
                              ///< \f$ \ln\sqrt{Q_{ch1}/Q_{ch0}} \f$ (PMI data)
    kPMIChargeAttCorrected,   ///< shows charge (photon count) spectra corrected for 
                              ///< the attenuation length: 
                              ///< \f$ Q_{chi} \cdot /\exp{(-z/\lambda_{att})}\f$ (PMI data)
    kPMIChargeAttCorrectedSum ///< shows sum of Ch0 and Ch1 charge (photon count) spectra
                              ///< both corrected for the attenuation length: 
                              ///< \f$ Q_{ch0}/\exp{\frac{z}
                              ///< {\lambda_{att}}} + Q_{ch1}/\exp{\frac{(L-z)}
                              ///< {\lambda_{att}}}\f$ (PMI data)
}; 

/// \file
/// Enumaeration representing different types of cuts for
/// the analyzed data. All cuts include selection of correct
/// address (module) for the requested channel.

/// Enumaeration representing different types of cuts for
/// the analyzed data. All cuts include selection of correct
/// address (module) for the requested channel.
enum class SFCutType
{
    kSpecCh0,    ///< Cut for Ch0 spectra, includes: BLsig<max, PE>0, 0>T0<590,
                 ///< TOT>0 and Amp<maxAmp
    kSpecCh0A,   ///< Cut for Ch0 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kSpecCh1,    ///< Cut for Ch1 spectra, includes: BLsig<max, PE>0, 0>T0<590,
                 ///< TOT>0 and Amp<maxAmp
    kSpecCh1A,   ///< Cut for Ch1 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kSpecCh2,    ///< Cut for Ch2 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kCombCh0Ch1, ///< Cut for any combination of Ch0 and Ch1, includes: BLsig<max,
                 ///< PE>0, 0>T0< 590, TOT>0 and Amp<maxAmp for both channels
    kT0Diff,     ///< Cut for T0 difference spectra, includes: BLsig<max, PE>min,
                 ///< 0>T0<590, TOT>0, Amp<maxAmp for both channels and \f$M_{LR}\f$>
                 ///< min and \f$M_{LR}\f$<max
    kT0DiffECut, ///< Cut for T0 difference spectra with energy cut, includes:
                 ///< BLsig<max, min>PE<max, 0>T0<590, TOT>0 for both channels
                 ///< and \f$M_{LR}\f$>min and \f$M_{LR}\f$<max
    kBLCh0,      ///< Cut for Ch0 base line histogram, module selection only
    kBLCh1,      ///< Cut for Ch1 base line histogram, module selection only
    kBLCh2,      ///< Cut for Ch2 base line histogram, module selection only
    
    //----- cuts for PMI measurements
    kPMISpecCh0,    ///< cut for Ch0 charge (photon count) spectrum in PMI measurements
                    ///< (module selection only)
    kPMISpecCh1,    ///< cut for Ch1 charge (photon count) spectrum in PMI measurements
                    ///< (module selection only)
    kPMICombCh0Ch1, ///< cut for combination of Ch0 and Ch1 spectra in PMI measurements
                    ///< (module selection only)
    kPMIT0Diff,     ///< cut for T0 difference spectra in PMI measurements, includes: 
                    ///< \f$M_{LR}\f$>min and \f$M_{LR}\f$<max
    kPMIT0DiffECut  ///< cut for T0 difference spectra with energy cut in PMI measurements,
                    ///< includes: min>PE<max for both channels and \f$M_{LR}\f$>min and 
                    ///< \f$M_{LR}\f$<max
};

/// Structure representing full address of a channel
/// according to the convention from sifi-framework

struct ChannelAddress
{
    int  fAddress = 0x0000; ///< channel address
    int  fChID    = -1;     ///< channel ID/number
    int  fModule  = -1;     ///< module number
    int  fLayer   = -1;     ///< layer number
    int  fFiber   = -1;     ///< fiber number
    char fSide    = ' ';    ///< side: 'l' fot left (ch0/ch2) and 'r' for right (ch1)
};

/// Class providing standarized and uniform set of selections and cuts for the
/// analyzed data. Selections and cuts are returned as TString and are consistent
/// with ROOT's TTree style. All commands are consistent with trees produced with
/// the sifi-framework software.

class SFDrawCommands : public TObject
{

  private:
    static void CheckCommand(TString string);

  public:
    /// Default constructor.
    SFDrawCommands(){};
    /// Default destructor.
    ~SFDrawCommands(){};

    static TString        GetSelectionName(SFSelectionType selection);
    static TString        GetSelection(SFSelectionType selection, int unique, int ch,
                                       std::vector<double> customNum = {});
    static TString        GetSelection(SFSelectionType selection, int unique,
                                       std::vector<double> customNum = {});
    static TString        GetCut(SFCutType cut, std::vector<double> customNum = {});
    static ChannelAddress GetChannelAddress(int ch);

    void Print(void);

    
    ClassDef(SFDrawCommands, 1)
};

#endif /* __SFDrawCommands_H_ */
