// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFFitResults.hh            *
// *              K. Rusiecka              *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *           Created in 2018             *
// *                                       *
// *****************************************

#ifndef __SFFitResults_H_
#define __SFFitResults_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TPaveText.h"

using namespace std;

/// This class stores results of the fit and the fitted function.

class SFFitResults : public TObject{
  
private:
  int     fStat;
  int     fNpar;
  double  fFastDecTime;
  double  fFastDecTimeErr;
  double  fSlowDecTime;
  double  fSlowDecTimeErr;
  double  fAmpFast;
  double  fAmpFastErr;
  double  fAmpSlow;
  double  fAmpSlowErr;
  double  fT0;
  double  fT0Err;
  double  fConst;
  double  fIslow;
  double  fIfast;
  double  fChi2;
  double  fNDF;
  double  fFitXmin;
  double  fFitXmax;
  TString fFormula;
  TString fName;
  TF1     *fFunction;
  vector <double> fParameters;
  vector <double> fParErrors;
  
public:
  SFFitResults();
  SFFitResults(TString name);
  SFFitResults(TString name, TF1 *fun);
  ~SFFitResults();
  
  void Reset(void);
  void Print(void);
  
  bool SetFromFunction(TF1 *fun);
  void SetFastDecTime(double t, double err);
  void SetSlowDecTime(double t, double err);
  void SetAmpFast(double amp, double err);
  void SetAmpSlow(double amp, double err);
  void SetT0(double t0, double err);
  void SetIntensities(double iSlow, double iFast);
  void SetFitRange(double xmin, double xmax);
  bool SetParameters(vector <double> par);
  bool SetParErrors(vector <double> parErr);
  void SetStat(int stat)           { fStat = stat; };
  void SetNpar(int npar)           { fNpar = npar; };
  void SetConst(double BL)         { fConst = BL; };
  void SetChi2(double chi2)        { fChi2 = chi2; };
  void SetNDF(double ndf)          { fNDF = ndf; };
  void SetFormula(TString formula) { fFormula = formula; };
  void SetName(TString name)       { fName = name; };
  
  bool            GetFastDecTime(double &t, double &err);
  bool            GetSlowDecTime(double &t, double &err);
  bool            GetAmpFast(double &ampFast, double &ampFastErr);
  bool            GetAmpSlow(double &ampSlow, double &ampSlowErr);
  bool            GetT0(double &t0, double &t0Err);
  bool            GetIntensities(double &iSlow, double &iFast);
  bool            GetFitRange(double &xmin, double &xmax);
  vector <TF1*>   GetCompFunctions(void);
  TPaveText*      GetResultsPave(void);
  int             GetStat(void)       { return fStat; };
  int             GetNpar(void)       { return fNpar; };
  double          GetConst(void)      { return fConst; };
  double          GetChi2(void)       { return fChi2; };
  double          GetNDF(void)        { return fNDF; };
  TString         GetFormula(void)    { return fFormula; };
  TString         GetName(void)       { return fName; };
  TF1*            GetFunction(void)   { return fFunction; };
  vector <double> GetParameters(void) { return fParameters; };
  vector <double> GetParErrors(void)  { return fParErrors; };
  
  ClassDef(SFFitResults,1)
};

#endif