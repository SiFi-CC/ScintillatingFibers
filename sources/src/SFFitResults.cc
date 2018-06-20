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
///Standard constructor (recommended).
SFFitResults::SFFitResults(TString name, TF1 *fun){
  fName = name;
  SetFromFunction(fun);
}
//------------------------------------------------------------------
///Default destructor.
SFFitResults::~SFFitResults(){
}
//------------------------------------------------------------------
void SFFitResults::Reset(void){
  fStat           = -1;
  fNpar           = -1;
  fFastDecTime    = -1;
  fFastDecTimeErr = -1;
  fSlowDecTime    = -1;
  fSlowDecTimeErr = -1;
  fAmpFast        = -1;
  fAmpFastErr     = -1;
  fAmpSlow        = -1;
  fAmpSlowErr     = -1;
  fT0             = -1;
  fT0Err          = -1;
  fConst          = -1;
  fIslow          = -1;
  fIfast          = -1;
  fChi2           = -1;
  fNDF            = -1;
  fFitXmin        = -1;
  fFitXmax        = -1;
  fFormula        = "dummy";
  fName           = "dummy";
  fFunction       = NULL;
  if(!fParameters.empty()) fParameters.clear();
  if(!fParErrors.empty())  fParErrors.clear();
}
//------------------------------------------------------------------
bool SFFitResults::SetFromFunction(TF1 *fun){
 
  if(fun==NULL){
   cout << "##### Error in SFFitResults::SetFromFunction()! Empty pointer!" << endl;
   return false;
  }
  
  //----- Setting unambiguous parameters 
  fFunction = fun;
  fNpar     = fFunction->GetNpar();
  fChi2     = fFunction->GetChisquare();
  fNDF      = fFunction->GetNDF();
  fFormula  = fFunction->GetExpFormula();
  fFunction->GetRange(fFitXmin,fFitXmax);
  
  for(int i=0; i<fNpar; i++){
   fParameters.push_back(fFunction->GetParameter(i));
   fParErrors.push_back(fFunction->GetParError(i));
  }

  fT0       = fParameters[1];
  fT0Err    = fParErrors[1];
  fConst    = fParameters[fNpar-1];
  
  //----- Recognizing fast and slow component
  //----- Assigning ambiguous parameters
  if(fParameters[2]<fParameters[4]){
    fFastDecTime    = fParameters[2];
    fSlowDecTime    = fParameters[4];
    fAmpFast        = fParameters[0];
    fAmpSlow        = fParameters[3];
    fFastDecTimeErr = fParErrors[2];
    fSlowDecTimeErr = fParErrors[4];
    fAmpFastErr     = fParErrors[0];
    fAmpSlowErr     = fParErrors[3];
  }
  else{
    fFastDecTime    = fParameters[4];
    fSlowDecTime    = fParameters[2];
    fAmpFast        = fParameters[3];
    fAmpSlow        = fParameters[0];
    fFastDecTimeErr = fParErrors[4];
    fSlowDecTimeErr = fParErrors[2];
    fAmpFastErr     = fParErrors[3];
    fAmpSlowErr     = fParErrors[0];
  }
  
  //----- Calculating intensities
  double denom = fAmpFast*fFastDecTime + fAmpSlow*fSlowDecTime;
  fIfast = (fAmpFast*fFastDecTime/denom)*100;
  fIslow = (fAmpSlow*fSlowDecTime/denom)*100;
    
  //----- Determining fit validity
  if(fAmpFast<0 || fAmpSlow<0){
    cout << "##### Warning in SFFitResults::SetFromFunction()!" << endl;
    cout << "\t Negative amplitude!" << endl;
    fStat = -1;
  }
  else if(fFastDecTime<0 || fFastDecTime>1E3){
    cout << "##### Warning in SFFitResults::SetFromFunction()!" << endl;
    cout << "\t Fast decay time out of range!" << endl;
    fStat = -1;
  }
  else if(fSlowDecTime<0 || fSlowDecTime>1E4){
    cout << "##### Warning in SFFitResults::SetFromFunction()!" << endl;
    cout << "\t Slow decay time out of range!" << endl;
    fStat = -1;
  }
  else
    fStat = 0;
  
  return true;
}
//------------------------------------------------------------------
void SFFitResults::SetFastDecTime(double t, double err){
  fFastDecTime    = t;
  fFastDecTimeErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetSlowDecTime(double t, double err){
   fSlowDecTime    = t;
   fSlowDecTimeErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetAmpFast(double amp, double err){
   fAmpFast    = amp;
   fAmpFastErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetAmpSlow(double amp, double err){
   fAmpSlow    = amp;
   fAmpSlowErr = err;
}
//------------------------------------------------------------------
void SFFitResults::SetT0(double t0, double err){
   fT0    = t0;
   fT0Err = err;
}
//------------------------------------------------------------------
void SFFitResults::SetIntensities(double iSlow, double iFast){
  fIslow = iSlow;
  fIfast = iFast;
}
//------------------------------------------------------------------
void SFFitResults::SetFitRange(double xmin, double xmax){
  fFitXmin = xmin;
  fFitXmax = xmax;
}
//------------------------------------------------------------------
bool SFFitResults::SetParameters(vector <double> par){
  if(par.empty()) return false;
  int size = par.size();
  for(int i=0; i<size; i++){
    fParameters.push_back(par[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::SetParErrors(vector <double> parErr){
  if(parErr.empty()) return false;
  int size = parErr.size();
  for(int i=0; i<size; i++){
    fParErrors.push_back(parErr[i]);
  }
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetFastDecTime(double &t, double &err){
  if(fFastDecTime==-1 || fFastDecTimeErr==-1){
    cout << "##### Error in SFFitResults::GetFastDecTime()!" << endl;
    return false;
  }
  t   = fFastDecTime; 
  err = fFastDecTimeErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetSlowDecTime(double &t, double &err){
  if(fSlowDecTime==-1 || fSlowDecTimeErr==1){
  cout << "##### Error in SFFitResults::GetSlowDecTime()!" << endl;
  return false;
  }
  t   = fSlowDecTime;
  err = fSlowDecTimeErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetAmpFast(double &ampFast, double &ampFastErr){
  if(fAmpFast==-1 || fAmpFastErr==-1){
    cout << "##### Error in SFFitResults::GetAmpFast()" << endl;
    return false;
  }
  ampFast    = fAmpFast;
  ampFastErr = fAmpFastErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetAmpSlow(double &ampSlow, double &ampSlowErr){
  if(fAmpSlow==-1 || fAmpSlowErr==-1){
    cout << "##### Error in SFFitResults::GetAmpSlow()" << endl;
    return false;
  }
  ampSlow    = fAmpSlow;
  ampSlowErr = fAmpSlowErr;
  return true;
}
//------------------------------------------------------------------
bool SFFitResults::GetT0(double &t0, double &t0Err){
  if(fT0==-1 || fT0Err==-1){
    cout << "##### Error in SFFitResults::GetT0()" << endl;
    return false;
  }
  t0    = fT0;
  t0Err = fT0Err;
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
bool SFFitResults::GetFitRange(double &xmin, double &xmax){
  if(fFitXmin==-1 || fFitXmax==-1){
    cout << "##### Error in SFFitResults::GetFitRange()" << endl;
    return false;
  }
  xmin = fFitXmin;
  xmax = fFitXmax;
  return true;
}
//------------------------------------------------------------------
vector <TF1*> SFFitResults::GetCompFunctions(void){
  
 vector <TF1*> functions;

 TF1 *funBL = new TF1("funBL","pol0",0,1024);
 funBL->FixParameter(0,fConst);
 funBL->SetLineColor(kMagenta+3);
 
 TF1 *funFastDecay = new TF1("funFastDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
 funFastDecay->FixParameter(0,fAmpFast);
 funFastDecay->FixParameter(1,fT0);
 funFastDecay->FixParameter(2,fFastDecTime);
 funFastDecay->SetLineColor(kMagenta);
 
 TF1 *funSlowDecay = new TF1("funSlowDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
 funSlowDecay->FixParameter(0,fAmpSlow);
 funSlowDecay->FixParameter(1,fT0);
 funSlowDecay->FixParameter(2,fSlowDecTime);
 funSlowDecay->SetLineColor(kMagenta-8);

 functions.push_back(funBL);
 functions.push_back(funFastDecay);
 functions.push_back(funSlowDecay);
 
 return functions;  
}
//------------------------------------------------------------------
TPaveText *SFFitResults::GetResultsPave(void){
 
  TPaveText *pave = new TPaveText(0.456,0.436,0.903,0.901,"NBNDC");
  pave->SetTextFont(42);

  pave->AddText(Form("A_fast = %.2f +/- %.2f",fAmpFast,fAmpFastErr));
  pave->AddText(Form("t_0 = %.2f +/- %.2f",fT0,fT0Err));
  pave->AddText(Form("tau_fast = %.2f +/- %.2f",fFastDecTime,fFastDecTimeErr));
  pave->AddText(Form("A_slow = %.2f +/- %.2f",fAmpSlow,fAmpSlowErr));
  pave->AddText(Form("tau_slow = %.2f +/- %.2f",fSlowDecTime,fSlowDecTimeErr));
  pave->AddText(Form("const = %.2f",fConst));
  pave->AddText(Form("Chi2/NDF = %.2f",fChi2/fNDF));
  pave->AddText(Form("I_fast = %.2f perc.",fIfast));
  pave->AddText(Form("I_slow = %.2f perc.",fIslow));
  return pave;
}
//------------------------------------------------------------------
void SFFitResults::Print(void){
  cout << "\n------------------------------------------------------" << endl;
  cout << "This is Print of SFFitResults class object " << fName << endl;
  cout << "Fitted function: " << fFormula << endl;
  cout << "Fit status: " << fStat << endl;
  cout << "Number of parameters: " << fNpar << endl;
  cout << "Chi2 = " << fChi2 << endl;
  cout << "Chi2/NDF = " << fChi2/fNDF << endl;
  cout << "Fitting range: " << fFitXmin << " - " << fFitXmax << endl;
  cout << "Fast decay time: " << fFastDecTime << " +/- " << fFastDecTimeErr << " ns" << endl;
  cout << "Slow decay time: " << fSlowDecTime << " +/- " << fSlowDecTimeErr << " ns" << endl;
  cout << "Fast component intensity: " << fIfast << " %" << endl;
  cout << "Slow component intensity: " << fIslow << " %" << endl;
  cout << "\nAll parameters: " << endl;
  for(int i=0; i<fNpar; i++){
   cout << setw(20);
   cout << "\t" << fParameters[i] << " +/- " << fParErrors[i] << endl;
  }
  cout << "------------------------------------------------------\n" << endl; 
  return;
}
//------------------------------------------------------------------