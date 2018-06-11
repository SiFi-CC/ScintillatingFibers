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
  fNpar           = -1;
  fDecayTime      = -1;
  fDecayTimeErr   = -1;
  fFastDecTime    = -1;
  fFastDecTimeErr = -1;
  fSlowDecTime    = -1;
  fSlowDecTimeErr = -1;
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
  
  double denom = 0;
  
  if(fNpar==4){
   fDecayTime    = fParameters[2];
   fDecayTimeErr = fParErrors[2];
  }
  else if(fNpar==6){
   fFastDecTime    = fParameters[2];
   fFastDecTimeErr = fParErrors[2];
   fSlowDecTime    = fParameters[4];
   fSlowDecTimeErr = fParErrors[4];
   denom = fParameters[0]*fParameters[2] + fParameters[3]*fParameters[4];
   fIfast = (fParameters[0]*fParameters[2]/denom)*100;
   fIslow = (fParameters[3]*fParameters[4]/denom)*100;
  }
  else{
    cout << "##### Error in SFFitResults::SetFromFunction()! Incorrect number of parameters!" << endl;
    cout << "\t fNpar for single decay = 4" << endl;
    cout << "\t fNpar for double decay = 6" << endl;
    cout << "\t fNpar = " << fNpar << endl;
    return false;
  }
  
  return true;
}
//------------------------------------------------------------------
void SFFitResults::SetDecayTime(double t, double err){
  fDecayTime    = t;
  fDecayTimeErr = err;
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
 
 if(fNpar==4){
  functions.push_back(fFunction);
  cout << "##### Warning in SFFitResults::GetCompFunctions()!" << endl;
  cout << "This is SFFitResults for single decay i.e. only one function available!" << endl;
 }
 else if(fNpar==6){
   TF1 *funFastDecay = new TF1("funFastDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
   funFastDecay->FixParameter(0,fParameters[0]);
   funFastDecay->FixParameter(1,fParameters[1]);
   funFastDecay->FixParameter(2,fParameters[2]);
   funFastDecay->SetLineColor(kMagenta);
   functions.push_back(funFastDecay);
   
   TF1 *funSlowDecay = new TF1("funSlowDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
   funSlowDecay->FixParameter(0,fParameters[3]);
   funSlowDecay->FixParameter(1,fParameters[1]);
   funSlowDecay->FixParameter(2,fParameters[4]);
   funSlowDecay->SetLineColor(kMagenta-8);
   functions.push_back(funSlowDecay);
 }
 else{
   cout << "##### Error in SFFitResults::GetCompFunctions()! Incorrect number of parameters!" <<endl;
 }
 
 return functions;  
}
//------------------------------------------------------------------
void SFFitResults::Print(void){
  cout << "\n------------------------------------------------------" << endl;
  cout << "This is Print of SFFitResults class object " << fName << endl;
  cout << "Fitted function: " << fFormula << endl;
  cout << "Number of parameters: " << fNpar << endl;
  cout << "Chi2 = " << fChi2 << endl;
  cout << "Chi2/NDF = " << fChi2/fNDF << endl;
  cout << "Fitting range: " << fFitXmin << " - " << fFitXmax << endl;
  if(fNpar==4) 
    cout << "Decay time: " << fDecayTime << " +/- " << fDecayTimeErr << " ns" << endl;
  if(fNpar==6){
    cout << "Fast decay time: " << fFastDecTime << " +/- " << fFastDecTimeErr << " ns" << endl;
    cout << "Slow decay time: " << fSlowDecTime << " +/- " << fSlowDecTimeErr << " ns" << endl;
    cout << "Fast component intensity: " << fIfast << " %" << endl;
    cout << "Slow component intensity: " << fIslow << " %" << endl;
  }
  cout << "\nAll parameters: " << endl;
  for(int i=0; i<fNpar; i++){
   cout << setw(20);
   cout << "\t" << fParameters[i] << " +/- " << fParErrors[i] << endl;
  }
  cout << "------------------------------------------------------\n" << endl; 
  return;
}
//------------------------------------------------------------------