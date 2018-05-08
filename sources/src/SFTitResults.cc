// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFFitResults.cc            *
// *              K. Rusiecka              *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "SFFitResults.hh"

ClassImp(SFFitResults);

//------------------------------------------------------------------
///Default constructor.
SFFitResults::SFFitResults(){
  cout << "##### Warning in SFFitResults constructor! You are using default constructor!" << endl;
  Reset();
}
//------------------------------------------------------------------
///Standard constructor (recommended).
SFFitResults::SFFitResults(TString name){
  Reset();
  fName = name;
}
//------------------------------------------------------------------
///Standasd constructor (recommended).
SFFitResults::SFFitResults(TString name, TF1 *riseFun, TF1 *decayFun){
  fName = name;
  SetFromFunctions(riseFun,decayFun);
}
//------------------------------------------------------------------
///Default destructor.
SFFitResults::~SFFitResults(){
}
//------------------------------------------------------------------
void SFFitResults::Reset(void){
  fDecayTime    = -1;
  fDecayTimeErr = -1;
  fFastDec      = -1;
  fFastDecErr   = -1;
  fSlowDec      = -1;
  fSlowDecErr   = -1;
  fRiseTime     = -1;
  fRiseTimeErr  = -1;
  fIslow        = -1;
  fIfast        = -1;
  fChi2NDFRise  = -1;
  fChi2NDFDec   = -1;
  fXminRise     = -1;
  fXmaxRise     = -1;
  fXminDec      = -1;
  fXmaxDec      = -1;
  fTsplit       = -1;
  fFormulaRise  = "dummy";
  fFormulaDecay = "dummy";
  fName         = "dummy";
  fDecayFun     = NULL;
  fRiseFun      = NULL;
  if(!fDecayPar.empty()) fDecayPar.clear();
  if(!fDecayParErr.empty()) fDecayParErr.clear();
  if(!fRisePar.empty()) fRisePar.clear();
  if(!fRiseParErr.empty()) fRiseParErr.clear();
}
//------------------------------------------------------------------
bool SFFitResults::SetFromFunctions(TF1 *riseFun, TF1 *decayFun){
 
  return true;
}
//------------------------------------------------------------------
void SFFitResults::SetDecayTime(double t, double err){
  fDecayTime    = t;
  fDecayTimeErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetFastDecTime(double t, double err){
  fFastDec    = t;
  fFastDecErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetSlowDecTime(double t, double err){
   fSlowDec    = t;
   fSlowDecErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetRiseTime(double t, double err){
  fRiseTime    = t;
  fRiseTimeErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetIntensities(double iSlow, double iFast){
  fIslow = iSlow;
  fIfast = iFast;
}
//------------------------------------------------------------------
void SFFitResults::SetRangesRise(double xmin, double xmax){
  fXminRise = xmin;
  fXmaxRise = xmax;
}
//------------------------------------------------------------------
void SFFitResults::SetRangesDecay(double xmin, double xmax){
  fXminDec = xmin;
  fXmaxDec = xmax;
}
//------------------------------------------------------------------
bool SFFitResults::SetDecayFun(TF1 *decFun){
  if(decFun==NULL) return false;
  fDecayFun = decFun;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetRiseFun(TF1 *riseFun){
  if(riseFun==NULL) return false;
  fRiseFun = riseFun;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetDecayPar(vector <double> par){
  if(par.empty()) return false;
  int size = par.size();
  for(int i=0; i<size; i++){
    fDecayPar.push_back(par[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetDecayParErrors(vector <double> parErr){
  if(parErr.empty()) return false;
  int size = parErr.size();
  for(int i=0; i<size; i++){
    fDecayParErr.push_back(parErr[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetRisePar(vector <double> par){
  if(par.empty()) return false;
  int size = par.size();
  for(int i=0; i<size; i++){
    fRisePar.push_back(par[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetRiseParErrors(vector <double> parErr){
  if(parErr.empty()) return false;
  int size = parErr.size();
  for(int i=0; i<size; i++){
    fRiseParErr.push_back(parErr[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetDecayTime(double &t, double &err){
  if(fDecayTime==-1 || fDecayTimeErr==-1){
    cout << "##### Error in SFFitResults::GetDecayTime()!" << endl;
    return false;
  }
  t   = fDecayTime;
  err = fDecayTimeErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetFastDecTime(double &t, double &err){
  if(fFastDec==-1 || fFastDecErr==-1){
    cout << "##### Error in SFFitResults::GetFastDecTime()!" << endl;
    return false;
  }
  t   = fFastDec; 
  err = fFastDecErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetSlowDecTime(double &t, double &err){
  if(fSlowDec==-1 || fSlowDecErr==1){
  cout << "##### Error in SFFitResults::GetSlowDecTime()!" << endl;
  return false;
  }
  t   = fSlowDec;
  err = fSlowDecErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetRiseTime(double &t, double &err){
  if(fRiseTime==-1 || fRiseTimeErr==-1){
    cout << "##### Error in SFFitResults::GetRiseTime()!" << endl;
    return false;
  }
  t   = fRiseTime;
  err = fRiseTimeErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetIntensities(double &iSlow, double &iFast){
  if(fIslow==-1 || fIfast==-1){
    cout << "##### Error in SFFitResults::GetIntensities()" << endl;
    return false;
  }
  iSlow = fIslow;
  iFast = fIfast;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetChi2NDF(double Chi2NDFRise, double Chi2NDFDec){
  if(fChi2NDFDec==-1 || fChi2NDFRise==-1){
    cout << "##### Error in SFFitResults::GetChi2NDF()" << endl;
    return false;
  }
  Chi2NDFRise = fChi2NDFRise;
  Chi2NDFDec  = fChi2NDFDec;
  return false;
}
//------------------------------------------------------------------
bool SFFitResults::GetRangesRise(double &xmin, double &xmax){
  if(fXminRise==-1 || fXmaxRise==-1){
    cout << "##### Error in SFFitResults::GetRangesRise()" << endl;
    return false;
  }
  xmin = fXminRise;
  xmax = fXmaxRise;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetRangesDecay(double &xmin, double &xmax){
  if(fXminDec==-1 || fXmaxDec==-1){
    cout << "##### Error in SFFitResults::GetRangesDecay()" << endl;
    return false;
  }
  xmin = fXminDec;
  xmax = fXmaxDec;
  return true;
}
//------------------------------------------------------------------
void SFFitResults::Print(void){
  
  cout << "\n------------------------------------------------------" << endl; 
  cout << "------------------------------------------------------\n" << endl; 
  
}
//------------------------------------------------------------------