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
enum class SFSelectionType{
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
     PEvsPEch2Correlation,  ///< shows 2D distribution for chosen channel i: \f$Q_{chi}\f$ vs. \f$Q_{ch2}\f$
     PEAttCorrected,        ///< shows PE spectrum corrected for attenuation length: \f$ Q_{chi} \cdot /\exp{(-z/\lambda_{att})}\f$
     PEAttCorrectedSum      ///< shows sum of Ch0 and Ch1 PE spectra, both corrected 
                            ///< for attenuation length: \f$ Q_{ch0}/\exp{\frac{z}{\lambda_{att}}} + Q_{ch1}/\exp{\frac{(L-z)}{\lambda_{att}}}\f$
};

/// Class providing standarized and uniform set of selections for analyzed 
/// data. Selections are returned as TString and are consistent with 
/// ROOT's TTree style.

class SFDrawCommands : public TObject{
    
private:
    static void CheckSelection(TString string);
    
public:
    /// Default constructor.
    SFDrawCommands() {};
    /// Default destructor.
    ~SFDrawCommands() {};
    
    static TString GetSelectionName(SFSelectionType selection);
    static TString GetSelection(SFSelectionType selection, int unique, 
                                int ch, std::vector <double> customNum={});
    static TString GetSelection(SFSelectionType selection, int unique, 
                                std::vector <double> customNum={});
    
    void Print(void);
    
    ClassDef(SFDrawCommands,1)
};

#endif
