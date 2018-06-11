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
#include "TMinuit.h"

using namespace std;

/// This class stores results of the fit and the fitted function.

class SFFitResults : public TObject{
  
private:
  int     fNpar;
  double  fDecayTime;
  double  fDecayTimeErr;
  double  fFastDecTime;
  double  fFastDecTimeErr;
  double  fSlowDecTime;
  double  fSlowDecTimeErr;
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
  void SetDecayTime(double t, double err);
  void SetFastDecTime(double t, double err);
  void SetSlowDecTime(double t, double err);
  void SetIntensities(double iSlow, double iFast);
  void SetFitRange(double xmin, double xmax);
  bool SetParameters(vector <double> par);
  bool SetParErrors(vector <double> parErr);
  void SetNpar(int npar)           { fNpar = npar; };
  void SetChi2(double chi2)        { fChi2 = chi2; };
  void SetNDF(double ndf)          { fNDF = ndf; };
  void SetFormula(TString formula) { fFormula = formula; };
  void SetName(TString name)       { fName = name; };
  
  bool            GetDecayTime(double &t, double &err);
  bool            GetFastDecTime(double &t, double &err);
  bool            GetSlowDecTime(double &t, double &err);
  bool            GetIntensities(double &iSlow, double &iFast);
  bool            GetFitRange(double &xmin, double &xmax);
  vector <TF1*>   GetCompFunctions(void);
  int             GetNpar(void)       { return fNpar; };
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