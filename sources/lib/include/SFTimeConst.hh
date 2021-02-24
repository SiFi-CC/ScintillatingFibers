// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimeConst.hh             *
// *       J. Kasper, K. Rusiecka          *
// *     kasper@physik.rwth-aachen.de      *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *           Created in 2018             *
// *                                       *
// *****************************************

#ifndef __SFTimeConst_H_
#define __SFTimeConst_H_ 1

#include "SFData.hh"
#include "SFFitResults.hh"
#include "SFTools.hh"
#include "SFResults.hh"

#include <Math/MinimizerOptions.h>
#include <TF1.h>
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>

#include <iostream>
#include <stdlib.h>
#include <string>

/// Class for determination of decay time constants from the averaged signals.
/// The accessed histograms are averaging of maximum of 50 signals.
/// In order to determine the time constants a sum of two exponential functions
/// is fitted to the signal TProfile histograms. Details of the fitting and the
/// function can be found in the presenation of KR posted on wiki
/// [LINK](http://bragg.if.uj.edu.pl/gccbwiki/images/5/5e/KR_20180604_TimeConstSummary.pdf) , on
/// slide 15. Double decay mode is assumed, i.e. we determine fast and slow decay time constants.
/// Additionally, only falling slope of the signal is fitted, i.e. the rise time is not determined.
/// Fitting results are stored in SFFitResults class objects. If function FitAllSignals() is called,
/// average values of time constants and intensities for the whole series are calculated.

// struct TimeConstResults
// {
// 
//     double fFastDecAv    = -1;
//     double fFastDecAvErr = -1;
// 
//     double fSlowDecAv    = -1;
//     double fSlowDecAvErr = -1;
// 
//     double fIfastAv = -1;
//     double fIslowAv = -1;
// 
//     std::vector<SFFitResults*> fResultsCh0;
//     std::vector<SFFitResults*> fResultsCh1;
// };

class SFTimeConst : public TObject
{

  private:
    int     fSeriesNo; ///< Number of analyzed fSeriesNo
    SFData* fData;     ///< Data of the measurement series
    double  fPE;       ///< Value of signals PE
    bool    fVerb;     ///< Verbose level: false - quiet, true - verbose

    std::vector<TProfile*>     fSignalsCh0; ///< Vector containing all signals from channel 0
    std::vector<TProfile*>     fSignalsCh1; ///< Vector containing all signals from channel 1
    std::vector<SFFitResults*> fFitResultsCh0;
    std::vector<SFFitResults*> fFitResultsCh1;
    
    SFResults* fResults;

  public:
    SFTimeConst();
    SFTimeConst(int seriesNo, double PE, bool verb);
    ~SFTimeConst();

    bool                       SetDetails(int seriesNo, double PE, bool verb);
    bool                       FitDecayTimeDouble(TProfile* signal, int ID);
    bool                       FitDecayTimeSingle(TProfile* signal, int ID);
    bool                       FitAllSignals(void);
    bool                       FitAllSignals(int ch);
    void                       Print(void);
    std::vector<TProfile*>     GetSignals(int ch);
    std::vector<SFFitResults*> GetFitResults(int ch);
    SFResults*                 GetResults(void) { return fResults; };

    ClassDef(SFTimeConst, 1)
};

#endif
