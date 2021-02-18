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

enum class SFSelectionType
{
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
    kBLSigma               ///< shows histogram of base line sigma
};

/// \file
/// Enumaeration representing different types of cuts for
/// the analyzed data. All cuts include selection of correct
/// address (module) for the requested channel.

enum class SFCutType
{
    kSpecCh0,    ///< Cut for ch0 spectra, includes: BLsig<max, PE>0, 0>T0<590,
                 ///< TOT>0 and Amp<maxAmp
    kSpecCh0A,   ///< Cut for ch0 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kSpecCh1,    ///< Cut for ch1 spectra, includes: BLsig<max, PE>0, 0>T0<590,
                 ///< TOT>0 and Amp<maxAmp
    kSpecCh1A,   ///< Cut for ch1 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kSpecCh2,    ///< Cut for ch2 spectra, includes: BLsig<max, PE>0, 0>T0<590
                 ///< and TOT>0
    kCombCh0Ch1, ///< Cut for any combination of ch0 and ch1, includes: BLsig<max,
                 ///< PE>0, 0>T0< 590, TOT>0 and Amp<maxAmp for both channels
    kT0Diff,     ///< Cut for T0 difference spectra, includes: BLsig<max, PE>min,
                 ///< 0>T0<590, TOT>0, Amp<maxAmp for both channels and \f$M_{LR}\f$>
                 ///< min and \f$M_{LR}\f$<max
    kT0DiffECut, ///< Cut for T0 difference spectra with energy cut, includes:
                 ///< BLsig<max, min>PE<max, 0>T0<590, TOT>0 for both channels
                 ///< and \f$M_{LR}\f$>min and \f$M_{LR}\f$<max
    kBLCh0,      ///< Cut for ch0 base line histogram, module selection only
    kBLCh1,      ///< Cut for ch1 base line histogram, module selection only
    kBLCh2       ///< Cut for ch2 base line histogram, module selection only
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

#endif
