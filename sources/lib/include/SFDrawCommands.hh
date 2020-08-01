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
#include "TObject.h"
#include "TString.h"
#include <iostream>

/// \file
/// Enumeration representing different types of selections
/// for the analyzed data:

enum class SFSelectionType
{
     PE,                    ///< shows PE (calibrated charge) spectrum
     Charge,                ///< shows uncalibrated charge spectrum
     Amplitude,             ///< shows amplitude spectrum
     T0,                    ///< shows T0 spectrum
     TOT,                   ///< shows TOT spectrum
     LogSqrtPERatio,        ///< shows PE ratio spectrum: \f$ \ln\sqrt{Q_{ch1}/Q_{ch0}} \f$
     T0Difference,          ///< shows T0 difference spectrum: \f$ T_{0 ch0} - T_{0 ch1} \f$
     PEAverage,             ///< shows spectrum of PE goemetric mean: \f$ \sqrt{Q_{ch0} \cdot Q_{ch1}} \f$
     AmplitudeAverage,      ///< shows spectrum of amplitude geometric mean: \f$ \sqrt{A_{ch0} \cdot A_{ch1}} \f$
     PECorrelation,         ///< shows 2D PE correlation spectrum: \f$Q_{ch0}\f$ vs. \f$Q_{ch1}\f$
     AmplitudeCorrelation,  ///< shows 2D amiplitude correlation spectrum: \f$A_{ch0}\f$ vs. \f$A_{ch1}\f$
     T0Correlation,         ///< shows 2D T0 correlation spectrum: \f$T_{0 ch0}\f$ vs. \f$T_{0 ch1}\f$
     AmpPECorrelation,      ///< shows 2D distribution for chosen channel i: \f$A_{chi}\f$ vs. \f$Q_{chi}\f$
     PEAttCorrected,        ///< shows PE spectrum corrected for attenuation length: \f$ Q_{chi} \cdot /\exp{(-z/\lambda_{att})}\f$
     PEAttCorrectedSum,     ///< shows sum of Ch0 and Ch1 PE spectra, both corrected 
                            ///< for attenuation length: \f$ Q_{ch0}/\exp{\frac{z}{\lambda_{att}}} + Q_{ch1}/\exp{\frac{(L-z)}{\lambda_{att}}}\f$
     BL,
     BLSigma
};

/// \file
/// Enumaeration representing different types of cuts for
/// the analyzed data. All cuts, except from custom expressions,
/// include selection of correct module and cut of base line sigma. 

enum class SFCutType
{
    SpecCh0,    ///< Cut for channel 0 spectra, includes PE>0, T0>0, T0<590 and TOT>0
    SpecCh0A,
    SpecCh1,    ///< Cut for channel 1 spectra, includes, PE>0, T0>0, T0<590 and TOT>0
    SpecCh1A,
    SpecCh2,    ///< Cut for channel 2 spectra, includes, PE>0, T0>0, T0<590 and TOT>0
    CombCh0Ch1, ///< Cut for any combination of channels 0 and 1, includes PE>0, T0>0, T0<590 and TOT>0 for both channels
    T0Diff,     ///< Cut for T0 difference spectra, includes PE>min, T0>0, T0<590 and TOT>0 for both channels and \f$ M_{LR} \f$  > min and \f$ M_{LR} \f$  < max
    T0DiffECut,  ///< Cut for T0 difference spectra including energy cut, includes PE>min, PE<max, T0>0, T0<590 and TOT>0 for both channels and \f$ M_{LR} \f$  > min and \f$ M_{LR} \f$  < max
    BLCh0,
    BLCh1,
    BLCh2
};

/// Structure representing full address of channel 
/// according to the convention from sifi-framework

struct ChannelAddress
{
    int  fAddress = 0x0000; ///< channel address
    int  fChID    = -1;     ///< channel ID/number
    int  fModule  = -1;     ///< module number
    int  fLayer   = -1;     ///< layer number
    int  fFiber   = -1;     ///< fiber number
    char fSide    = ' ';    ///< side: 'l' fot left (ch0) and 'r' for right (ch1)
};

/// Class providing standarized and uniform set of selections and cuts for the 
/// analyzed data. Selections and cuts are returned as TString and are consistent 
/// with ROOT's TTree style. All commands are consistent with trees produced with
/// the sifi-framework software.

class SFDrawCommands : public TObject{
    
private:
    static void CheckCommand(TString string);
    
public:
    /// Default constructor.
    SFDrawCommands() {};
    /// Default destructor.
    ~SFDrawCommands() {};
    
    static TString        GetSelectionName(SFSelectionType selection);
    static TString        GetSelection(SFSelectionType selection, int unique, 
                                       int ch, std::vector <double> customNum={});
    static TString        GetSelection(SFSelectionType selection, int unique, 
                                       std::vector <double> customNum={});
    static TString        GetCut(SFCutType cut, std::vector <double> customNum={});
    static ChannelAddress GetChannelAddress(int ch);
    
    void Print(void);
    
    ClassDef(SFDrawCommands,1)
};

#endif
