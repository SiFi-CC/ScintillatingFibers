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

class SFTimeConst : public TObject{
  
private:
  int     fSeriesNo;	///< Number of analyzed fSeriesNo
  SFData  *fData;	///< Data of the measurement series  
  double  fPE;		///< Value of signals PE
  bool    fVerb;	///< Verbose level: false - quiet, true - verbose
  
  vector <TProfile*>     fSignals;
  vector <SFFitResults*> fResultsSingle;
  vector <SFFitResults*> fResultsDouble;
  
  int    GetIndex(double position);
  
public:
  SFTimeConst();
  SFTimeConst(int seriesNo, double PE, bool verb);
  ~SFTimeConst();
  
  bool          SetDetails(int seriesNo, double PE, bool verb);
  bool          FitDecTimeSingle(TProfile *signal, double position);
  bool          FitDecTimeDouble(TProfile *signal, double position);
  bool          FitAllSignals(void);
  bool          FitAllSignals(TString option);
  bool          FitSingleSignal(double position);
  bool          FitSingleSignal(double position, TString option);
  TProfile*     GetSingleSignal(double position);
  SFFitResults* GetSingleResult(double position, TString opt);
  void          Reset(void);
  void          Print(void);
  
  vector <TProfile*>     GetAllSignals(void) { return fSignals; };
  vector <SFFitResults*> GetAllResults(TString opt);
  
  ClassDef(SFTimeConst,1)
  
};

#endif