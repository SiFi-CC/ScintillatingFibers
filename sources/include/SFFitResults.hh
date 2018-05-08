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

using namespace std;

class SFFitResults : public TObject{
  
private:
  double  fDecayTime;
  double  fDecayTimeErr;
  double  fFastDec;
  double  fFastDecErr;
  double  fSlowDec;
  double  fSlowDecErr;
  double  fRiseTime;
  double  fRiseTimeErr;
  double  fIslow;
  double  fIfast;
  double  fChi2NDFRise;
  double  fChi2NDFDec;
  double  fXminRise;
  double  fXmaxRise;
  double  fXminDec;
  double  fXmaxDec;
  double  fTsplit;
  TString fFormulaRise;
  TString fFormulaDecay;
  TString fName;
  TF1     *fDecayFun;
  TF1     *fRiseFun;
  vector <double> fDecayPar;
  vector <double> fDecayParErr;
  vector <double> fRisePar;
  vector <double> fRiseParErr;
  
public:
  SFFitResults();
  SFFitResults(TString name);
  SFFitResults(TString name, TF1 *riseFun, TF1 *decayFun);
  ~SFFitResults();
  
  void Reset();
  void Print();
  
  bool SetFromFunctions(TF1 *riseFun, TF1 *decayFun);
  void SetDecayTime(double t, double err);
  void SetFastDecTime(double t, double err);
  void SetSlowDecTime(double t, double err);
  void SetRiseTime(double t, double err);
  void SetIntensities(double iSlow, double iFast);
  void SetRiseChi2NDF(double Chi2NDF)   { fChi2NDFRise = Chi2NDF; }; 
  void SetDecayChi2NDF(double Chi2NDF)  { fChi2NDFDec = Chi2NDF; };
  void SetRangesRise(double xmin, double xmax);
  void SetRangesDecay(double xmin, double xmax);
  void SetSplitTime(double tSplit)      { fTsplit = tSplit; };
  void SetFormulaRise(TString formula)  { fFormulaRise = formula; };
  void SetFormulaDecay(TString formula) { fFormulaDecay = formula; };
  void SetName(TString name)            { fName = name; };
  bool SetDecayFun(TF1 *decFun);
  bool SetRiseFun(TF1 *riseFun);
  bool SetDecayPar(vector <double> par);
  bool SetDecayParErrors(vector <double> parErr);
  bool SetRisePar(vector <double> par);
  bool SetRiseParErrors(vector <double> parErr);
  
  bool    GetDecayTime(double &t, double &err);
  bool    GetFastDecTime(double &t, double &err);
  bool    GetSlowDecTime(double &t, double &err);
  bool    GetRiseTime(double &t, double &err);
  bool    GetIntensities(double &iSlow, double &iFast);
  bool    GetChi2NDF(double Chi2NDFRise, double Chi2NDFDec);
  bool    GetRangesRise(double &xmin, double &xmax);
  bool    GetRangesDecay(double &xmin, double &xmax);
  double  GetSplitTime(void)           {return fTsplit; };
  TString GetFormulaRise(void)         { return fFormulaRise; };
  TString GetFormulaDecay(void)        { return fFormulaDecay; };
  TString GetName(void)                { return fName; };
  TF1*    GetDecayFunction(void)       { return fDecayFun; };
  TF1*    GetRiseFunction(void)        { return fRiseFun; };
  vector <double> GetDecayPar(void)    { return fDecayPar; };
  vector <double> GetDecayParErr(void) { return fDecayParErr; };
  vector <double> GetRisePar(void)     { return fRisePar; };
  vector <double> GetRiseParErr(void)  { return fRiseParErr; };
  
  ClassDef(SFFitResults,1)
};

#endif