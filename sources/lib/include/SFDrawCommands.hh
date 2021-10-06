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
    kDDPE,                   ///< shows PE (calibrated charge) spectrum
    kDDCharge,               ///< shows uncalibrated charge spectrum
    kDDAmplitude,            ///< shows amplitude spectrum
    kDDT0,                   ///< shows T0 spectrum
    kDDTOT,                  ///< shows TOT spectrum
    kDDLogSqrtPERatio,       ///< shows PE ratio spectrum:
                             ///< \f$ \ln\sqrt{Q_{R}/Q_{L}} \f$
    kDDT0Difference,         ///< shows T0 difference spectrum:
                             ///< \f$ T_{0 L} - T_{0 R} \f$
    kDDPEAverage,            ///< shows spectrum of PE goemetric mean:
                             ///< \f$ \sqrt{Q_{L} \cdot Q_{R}} \f$
    kDDAmplitudeAverage,     ///< shows spectrum of amplitude geometric mean:
                             ///< \f$ \sqrt{A_{L} \cdot A_{R}} \f$
    kDDPECorrelation,        ///< shows 2D PE correlation spectrum:
                             ///< \f$Q_{L}\f$ vs. \f$Q_{R}\f$
    kDDAmplitudeCorrelation, ///< shows 2D amiplitude correlation spectrum:
                             ///< \f$A_{L}\f$ vs. \f$A_{R}\f$
    kDDT0Correlation,        ///< shows 2D T0 correlation spectrum:
                             ///< \f$T_{0 L}\f$ vs. \f$T_{0 R}\f$
    kDDAmpPECorrelation,     ///< shows 2D distribution for chosen channel i:
                             ///< \f$A_{chi}\f$ vs. \f$Q_{chi}\f$
    kDDPEAttCorrected,       ///< shows PE spectrum corrected for attenuation length:
                             ///< \f$ Q_{chi} \cdot /\exp{(-z/\lambda_{att})}\f$
    kDDPEAttCorrectedSum,    ///< shows sum of L and R PE spectra, both corrected
                             ///< for attenuation length: \f$ Q_{L}/\exp{\frac{z}
                             ///< {\lambda_{att}}} + Q_{R}/\exp{\frac{(L-z)}
                             ///< {\lambda_{att}}}\f$
    kDDBL,                   ///< shows base line histogram
    kDDBLSigma,              ///< shows histogram of base line sigma
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
    ///----- cuts for single fiber measurements (Desktop Digitizer)
    kDDSpec,       ///< Cut for spectra, includes: BLsig<max, PE>0, 0>T0<590,
                   ///< TOT>0 and Amp<maxAmp
    kDDSpecA,      ///< Cut for spectra, includes: BLsig<max, PE>0, 0>T0<590
                   ///< and TOT>0
    kDDCombLR,     ///< Cut for any combination of L and R, includes: BLsig<max,
                   ///< PE>0, 0>T0< 590, TOT>0 and Amp<maxAmp for both channels
    kDDT0Diff,     ///< Cut for T0 difference spectra, includes: BLsig<max, PE>min,
                   ///< 0>T0<590, TOT>0, Amp<maxAmp for both channels and \f$M_{LR}\f$>
                   ///< min and \f$M_{LR}\f$<max
    kDDT0DiffECut, ///< Cut for T0 difference spectra with energy cut, includes:
                   ///< BLsig<max, min>PE<max, 0>T0<590, TOT>0 for both channels
                   ///< and \f$M_{LR}\f$>min and \f$M_{LR}\f$<max
    kDDBL,         ///< Cut for Ch0 base line histogram, module selection only
};

/// Structure representing full address of a channel
/// according to the convention from sifi-framework

struct SFChAddr
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

namespace SFDrawCommands
{
    void     CheckCommand(TString string);
    TString  GetSelectionName(SFSelectionType selection);
    
    TString  GetSelection(SFSelectionType selection,
                          int unique,
                          SFChAddr chAddr,
                          std::vector<double> customNum = {});
    TString  GetSelection(SFSelectionType selection,
                          int unique,
                          std::vector<double> customNum = {});
    TString  GetCut(SFCutType cut,
                    SFChAddr chAddr,
                    std::vector<double> customNum = {});
    
};

#endif /* __SFDrawCommands_H_ */
