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

/// This class stores results of the fit and the fitted function itself.
/// This is assistant class for SFTimeConst.

class SFFitResults : public TObject{
  
private:
  int     fStat;		///< Fit status. 0 if fit is succesfull, -1 if it failed
  int     fNpar;		///< Number of parameters
  double  fFastDecTime;		///< Fast decay time [ns]
  double  fFastDecTimeErr;	///< Uncertainty of the fast decay time [ns] 
  double  fSlowDecTime;		///< SLow decay time [ns]
  double  fSlowDecTimeErr;	///< Uncertainty of the slow decay time [n]
  double  fAmpFast;		///< Amplitude of the fast decay component
  double  fAmpFastErr;		///< Uncertainty of the fast component amplitude
  double  fAmpSlow;		///< Amplitude of the slow decay component
  double  fAmpSlowErr;		///< Uncertainty of the slow component amplitude
  double  fT0;			///< Time offset
  double  fT0Err;		///< Uncertainty of the time offset
  double  fConst;		///< Constant i.e. base line
  double  fIslow;		///< Intensity of the slow component [%]
  double  fIfast;		///< Intensity of the fast component [%]
  double  fChi2;		///< Chi square
  double  fNDF;			///< Number of degrees of freedom
  double  fFitXmin;		///< Minimum of the fitting range
  double  fFitXmax;		///< Maximum of the fitting range
  TString fFormula;		///< Formula of the fitted function
  TString fName;		///< Name of the SFFitResults object
  TF1     *fFunction;		///< Fitted function
  vector <double> fParameters;	///< Vector containing all parameters
  vector <double> fParErrors;	///< Vector containing errors of all parameters
  
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
  ///Sets fitting status. Fit status can be 0 
  ///if it was successfull or -1 if it failed
  void SetStat(int stat)           { fStat = stat; };
  ///Sets number of parameters.
  void SetNpar(int npar)           { fNpar = npar; };
  ///Sets constant (base line).
  void SetConst(double BL)         { fConst = BL; };
  ///Sets Chi square
  void SetChi2(double chi2)        { fChi2 = chi2; };
  ///Sets number of degrees of freedom
  void SetNDF(double ndf)          { fNDF = ndf; };
  ///Sets formula of the fitted function.
  void SetFormula(TString formula) { fFormula = formula; };
  ///Sets name of the SFFitResults class object.
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
  ///Returns fit status. If 0 is returned - fit was succesfull,
  ///if -1 was returned - fit failed.
  int             GetStat(void)       { return fStat; };
  ///Returns number of parameters.
  int             GetNpar(void)       { return fNpar; };
  ///Returns constant (base line).
  double          GetConst(void)      { return fConst; };
  ///Returns Chi square.
  double          GetChi2(void)       { return fChi2; };
  ///Returns number of degrees of freedom.
  double          GetNDF(void)        { return fNDF; };
  ///Returns formula of the fitted function.
  TString         GetFormula(void)    { return fFormula; };
  ///Returns name of the SFFitResults class object.
  TString         GetName(void)       { return fName; };
  ///Returns pointer to the fitted function.
  TF1*            GetFunction(void)   { return fFunction; };
  ///Returns vector containing fitted parameters.
  vector <double> GetParameters(void) { return fParameters; };
  ///Returns vector containing errors of parameters.
  vector <double> GetParErrors(void)  { return fParErrors; };
  
  ClassDef(SFFitResults,1)
};

#endif