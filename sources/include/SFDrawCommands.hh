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

/// Enumeration representing different types of selections
/// for the analyzed data:
enum class SFSelectionType{
     PE,                    //< shows PE (calibrated charge) spectrum
     Charge,                //< shows uncalibrated charge spectrum
     Amplitude,             //< shows amplitude spectrum
     T0,                    //< shows T0 spectrum
     TOT,                   //< shows TOT spectrum
     LogSqrtPERatio,        //< shows PE ratio spectrum
     T0Difference,          //< shows T0 difference spectrum
     PEAverage,             //< shows spectrum of PE goemetric mean
     AmplitudeAverage,      //< shows spectrum of amplitude geometric mean
     PECorrelation,         //< shows 2D PE correlation spectrum
     AmplitudeCorrelation,  //< shows 2D amiplitude correlation spectrum
     T0Correlation,         //< shows 2D T0 correlation spectrum
     PEAttCorrected,        //< shows PE spectrum corrected for attenuation length
     PEAttCorrectedSum      //< shows sum of Ch0 and Ch1 PE spectra, both corrected 
                            //< for attenuation length
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
