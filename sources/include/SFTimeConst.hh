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
#include "TObject.h"
#include "TString.h"
#include "TProfile.h"
#include "TF1.h"
#include "Math/MinimizerOptions.h"

#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

///Class for determination of decay time constants from the averaged signals.
///The accessed histograms are averaging of maximum of 50 signals. 
///In order to determine the time constants a sum of two exponential functions 
///is fitted to the signal TProfile histograms. Details of the fitting and the
///function can be found in the presenation of KR posted on wiki [LINK](http://bragg.if.uj.edu.pl/gccbwiki/images/5/5e/KR_20180604_TimeConstSummary.pdf)
///, on slide 15. Double decay mode is assumed, i.e. we determine fast and slow 
///decay time constants. Additionally, only falling slope of the signal is fitted, 
///i.e. the rise time is not determined. Fitting results are stored in SFFitResults 
///class objects. If function FitAllSignals() is called, average values of time constants 
///and intensities for the whole series are calculated.

class SFTimeConst : public TObject{
  
private:
  int     fSeriesNo;	///< Number of analyzed fSeriesNo
  SFData  *fData;	///< Data of the measurement series  
  double  fPE;		///< Value of signals PE
  bool    fVerb;	///< Verbose level: false - quiet, true - verbose
  
  vector <TProfile*>     fSignalsCh0;	///< Vector containing all signals from channel 0
  vector <TProfile*>     fSignalsCh1;	///< Vector containing all signals from channel 1 
  vector <SFFitResults*> fResultsCh0;	///< Vector containing fit results for signals of channel 0
  vector <SFFitResults*> fResultsCh1;	///< Vector containing fit results for signals of channel 1

  int    GetIndex(double position);
  
public:
  SFTimeConst();
  SFTimeConst(int seriesNo, double PE, bool verb);
  ~SFTimeConst();
  
  bool          SetDetails(int seriesNo, double PE, bool verb);
  bool          FitDecayTime(TProfile *signal, double position);
  bool          FitAllSignals(void);
  bool          FitAllSignals(int ch);
  bool          FitSingleSignal(int ch, double position);
  TProfile*     GetSingleSignal(int ch, double position);
  SFFitResults* GetSingleResult(int ch, double position);
  void          Reset(void);
  void          Print(void);
  vector <TProfile*>     GetAllSignals(int ch);
  vector <SFFitResults*> GetAllResults(int ch);
  
  ClassDef(SFTimeConst,1)
  
};

#endif