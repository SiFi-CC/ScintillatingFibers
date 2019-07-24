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
     T0Correlation          //< shows 2D T0 correlation spectrum
};

/// Enumeration representing different types of cuts for the
/// analyzed data:
enum class SFCutType{
     //RefOnly511Tight,   //< tight cut on 511 keV peak in the reference detector
     //RefOnly511Loose,   //< loose cut on 511 keV peak in the reference detector
     //RefBelow511,       //< events only beleow 511 keV peak in the reference detector
     PEAbove0,          //< events with the PE value larger than 0
     T0Valid,           //< events with T0 larger than -100
     T0Above0,          //< events with the T0 value above 0
     TOTValid,          //< events with the TOT value above -100
     T0ChBelowMax,      //< events with the T0 value below some maximum value, chosen channel
     T0BelowMax,        //< events with the T0 value below some maximum value, both channels
     ScatteredEvents,   //< rejection of scattered events
     Only511,           //< events only from 511 keV peak
     Empty              //< empty cut
};

/// Class providing standarized and uniform set of selections and cuts for analyzed 
/// data. Both, selections and cuts, are returned as TString and are consistent with 
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
    static TString GetSelection(SFSelectionType selection, int unique, int ch);
    static TString GetSelection(SFSelectionType selection, int unique);
    //static TString GetSelection(SFSelectionType selection, int unique, 
    //                            std::vector <double> customNumbers);
    
    static TString GetCut(SFCutType cut);
    static TString GetCut(SFCutType cut, std::vector <double> customNumbers, int ch);
    static TString GetCut(SFCutType cut, std::vector <double> customNumbers);
    
    void Print(void);
    
    ClassDef(SFDrawCommands,1)
};

#endif
