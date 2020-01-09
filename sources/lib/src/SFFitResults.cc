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
/// Default constructor.
SFFitResults::SFFitResults(): fComponents(-1),
                              fStat(100),
                              fNpar(-1),
                              fDecTime(-1),
                              fDecTimeErr(-1),
                              fFastDecTime(-1),
                              fFastDecTimeErr(-1),
                              fSlowDecTime(-1),
                              fSlowDecTimeErr(-1),
                              fAmp(-1),
                              fAmpErr(-1),
                              fAmpFast(-1),
                              fAmpFastErr(-1),
                              fAmpSlow(-1),
                              fAmpSlowErr(-1),
                              fT0(-1),
                              fT0Err(-1),
                              fConst(-1),
                              fIslow(-1),
                              fIfast(-1),
                              fChi2(-1),
                              fNDF(-1),
                              fFitXmin(-1),
                              fFitXmax(-1),
                              fFormula("dummy"),
                              fName("dummy"),
                              fFunction(nullptr) {
                                  
  std::cout << "##### Warning in SFFitResults constructor!" << std::endl;
  std::cout << "You are using default constructor!" << std::endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
SFFitResults::SFFitResults(TString name): fComponents(-1),
                                          fStat(100),
                                          fNpar(-1),
                                          fDecTime(-1),
                                          fDecTimeErr(-1),
                                          fFastDecTime(-1),
                                          fFastDecTimeErr(-1),
                                          fSlowDecTime(-1),
                                          fSlowDecTimeErr(-1),
                                          fAmp(-1),
                                          fAmpErr(-1),
                                          fAmpFast(-1),
                                          fAmpFastErr(-1),
                                          fAmpSlow(-1),
                                          fAmpSlowErr(-1),
                                          fT0(-1),
                                          fT0Err(-1),
                                          fConst(-1),
                                          fIslow(-1),
                                          fIfast(-1),
                                          fChi2(-1),
                                          fNDF(-1),
                                          fFitXmin(-1),
                                          fFitXmax(-1),
                                          fFormula("dummy"),
                                          fName(name),
                                          fFunction(nullptr) {

}
//------------------------------------------------------------------
/// Standard constructor (recommended).
SFFitResults::SFFitResults(TString name, TF1 *fun): fComponents(-1),
                                                    fStat(100),
                                                    fNpar(-1),
                                                    fDecTime(-1),
                                                    fDecTimeErr(-1),
                                                    fFastDecTime(-1),
                                                    fFastDecTimeErr(-1),
                                                    fSlowDecTime(-1),
                                                    fSlowDecTimeErr(-1),
                                                    fAmp(-1),
                                                    fAmpErr(-1),
                                                    fAmpFast(-1),
                                                    fAmpFastErr(-1),
                                                    fAmpSlow(-1),
                                                    fAmpSlowErr(-1),
                                                    fT0(-1),
                                                    fT0Err(-1),
                                                    fConst(-1),
                                                    fIslow(-1),
                                                    fIfast(-1),
                                                    fChi2(-1),
                                                    fNDF(-1),
                                                    fFitXmin(-1),
                                                    fFitXmax(-1),
                                                    fFormula("dummy"),
                                                    fName(name),
                                                    fFunction(fun) {

  int stat = SetFromFunction(fun);
  
  if(!stat){
      throw "##### Exception in SFFitResults constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
SFFitResults::~SFFitResults(){
}
//------------------------------------------------------------------
/// Sets values of the private members of the class based on the passed
/// pointer to the fitted function. Inside this function fast and slow
/// decay components are recognized, their intensities are calculculated
/// and fit validity is chacked.
/// \param fun - fitted function
bool SFFitResults::SetFromFunction(TF1 *fun){
 
  if(fun==nullptr){
   std::cerr << "##### Error in SFFitResults::SetFromFunction()! Empty pointer!" << std::endl;
   return false;
  }
  
  if(fFunction==nullptr)
    fFunction = fun;
  
  //----- Setting unambiguous parameters 
  fChi2     = fFunction->GetChisquare();
  fNDF      = fFunction->GetNDF();
  fNpar     = fFunction->GetNpar();
  fFunction->GetRange(fFitXmin, fFitXmax);
  
  if(fNpar==4){
    fComponents = 1;
    fFormula = "A*exp(-(x-t0)/tau) + const";
  }
  else if(fNpar==6){
    fComponents = 2;
    fFormula    = "A_fast*exp(-(x-t0)/tau_fast) + A_slow*exp(-(x-t0)/tau_slow) + const";
  }
  else{
    std::cerr << "##### Error in SFFitResults::SetFromFunction()!" << std::endl;
    std::cerr << "Unknown function!" << std::endl;
    return false;
  }
  
  for(int i=0; i<fNpar; i++){
   fParameters.push_back(fFunction->GetParameter(i));
   fParErrors.push_back(fFunction->GetParError(i));
  }

  SetT0(fParameters[1], fParErrors[1]);
  SetConst(fParameters[fNpar-1]);
  
  //----- Assigning ambiguous parameters
  //----- In double decay: Recognizing fast and slow component
  if(fComponents==1){
    SetAmp(fParameters[0], fParErrors[0]);
    SetDecTime(fParameters[2], fParErrors[2]);
  }
  else if(fComponents==2){
    if(fParameters[2]<fParameters[4]){
      SetFastDecTime(fParameters[2],fParErrors[2]);
      SetSlowDecTime(fParameters[4],fParErrors[4]);
      SetAmpFast(fParameters[0],fParErrors[0]);
      SetAmpSlow(fParameters[3],fParErrors[3]);
    }
    else{
      SetFastDecTime(fParameters[4],fParErrors[4]);
      SetSlowDecTime(fParameters[2],fParErrors[2]);
      SetAmpFast(fParameters[3],fParErrors[3]);
      SetAmpSlow(fParameters[0],fParErrors[0]);
    }
  
  //----- Calculating intensities
    double denom = fAmpFast*fFastDecTime + fAmpSlow*fSlowDecTime;
    fIfast = (fAmpFast*fFastDecTime/denom)*100;
    fIslow = (fAmpSlow*fSlowDecTime/denom)*100;
  }
  
  //----- Determining fit validity
  if(fComponents==1){
    if(fAmp<0){
      std::cout << "##### Warning in SFFitResults::SetFromFunction()!" << std::endl;
      std::cout << "\t Negative amplitude!" << std::endl;
      fStat = -1;
    }
    else if(fDecTime<0 || fDecTime>1E2){
      std::cout << "##### Warning in SFFitResults::SetFromFunction()!" << std::endl;
      std::cout << "\t Decay time out of range!" << std::endl;
      fStat = -1;
    }
    else fStat = 0;
  }
  
  if(fComponents==2){
    if(fAmpFast<0 || fAmpSlow<0){
      std::cout << "##### Warning in SFFitResults::SetFromFunction()!" << std::endl;
      std::cout << "\t Negative amplitude!" << std::endl;
      fStat = -1;
    }
    else if(fFastDecTime<0 || fFastDecTime>1E3){
      std::cout << "##### Warning in SFFitResults::SetFromFunction()!" << std::endl;
      std::cout << "\t Fast decay time out of range!" << std::endl;
      fStat = -1;
    }
    else if(fSlowDecTime<0 || fSlowDecTime>1E4){
      std::cout << "##### Warning in SFFitResults::SetFromFunction()!" << std::endl;
      std::cout << "\t Slow decay time out of range!" << std::endl;
      fStat = -1;
    }
    else
      fStat = 0;
  }
  
  return true;
}
//------------------------------------------------------------------
///Sets decay time and its uncertainty (single decay mode).
///\param t - decay time [ns]
///\param err - uncertainty
void SFFitResults::SetDecTime(double t, double err){
  fDecTime = t;
  fDecTimeErr = err;
  return;
}
//------------------------------------------------------------------
///Sets fast decay time and its uncertainty.
///\param t - fast decay time [ns]
///\param err - uncertainty
void SFFitResults::SetFastDecTime(double t, double err){
  fFastDecTime    = t;
  fFastDecTimeErr = err;
  return;
}
//------------------------------------------------------------------
///Sets slow decay time and its uncertainty.
///\param t - slow decay time [ns]
///\param err - uncertainty
void SFFitResults::SetSlowDecTime(double t, double err){
   fSlowDecTime    = t;
   fSlowDecTimeErr = err;
   return;
}
//------------------------------------------------------------------
///Sets amplitude of the decay function (single decay mode).
///\param amp - amplitude
///\param err - uncertainty
void SFFitResults::SetAmp(double amp, double err){
  fAmp    = amp;
  fAmpErr = err;
  return;
}
//------------------------------------------------------------------
///Sets amplitude of the of the fast decay component and its
///uncertainty.
///\param amp - amplitude
///\param err - uncertainty
void SFFitResults::SetAmpFast(double amp, double err){
   fAmpFast    = amp;
   fAmpFastErr = err;
   return;
}
//------------------------------------------------------------------
///Sets amplitude of the slow decay component and its 
///uncertainty.
///\param amp - amplitude 
///\param err - uncertainty 
void SFFitResults::SetAmpSlow(double amp, double err){
   fAmpSlow    = amp;
   fAmpSlowErr = err;
   return;
}
//------------------------------------------------------------------
///Sets time offset and its uncertainty.
///\param t0 - time offset [ns]
///\param err - uncertainty 
void SFFitResults::SetT0(double t0, double err){
   fT0    = t0;
   fT0Err = err;
   return;
}
//------------------------------------------------------------------
///Sets intensities of decay components.
///\param iSlow - intensity of the slow component [%]
///\param iFast - intensity of the fast component [%]
void SFFitResults::SetIntensities(double iSlow, double iFast){
  fIslow = iSlow;
  fIfast = iFast;
  return;
}
//------------------------------------------------------------------
///Sets fitting range. 
///\param xmin - lower fitting range
///\param xmax - upper fitting range 
void SFFitResults::SetFitRange(double xmin, double xmax){
  fFitXmin = xmin;
  fFitXmax = xmax;
  return;
}
//------------------------------------------------------------------
///Fills vector fParameters containing fitted parameters. 
///\param par - vector of parameters
bool SFFitResults::SetParameters(std::vector <double> par){
  if(par.empty()) return false;
  int size = par.size();
  for(int i=0; i<size; i++){
    fParameters.push_back(par[i]);
  }
  return true;
}
//------------------------------------------------------------------
///Fills vector fParErrors containing errors of parameters.
///\param parErr - vector of parameters errors
bool SFFitResults::SetParErrors(std::vector <double> parErr){
  if(parErr.empty()) return false;
  int size = parErr.size();
  for(int i=0; i<size; i++){
    fParErrors.push_back(parErr[i]);
  }
  return true;
}
//------------------------------------------------------------------
///Returns references to the decay time and its uncertainty (single decay mode)
///\param t - decay time [ns]
///\param err - uncertainty
bool SFFitResults::GetDecTime(double &t, double &err){
  if(fDecTime==-1 || fDecTimeErr==-1){
    std::cerr << "##### Error in SFFitResults::GetDecTime()!" << std::endl;
    return false;
  }
  t   = fDecTime;
  err = fDecTimeErr; 
}
//------------------------------------------------------------------
///Returns references to the fast decay time and its uncertainty.
///\param t - fast decay time [ns]
///\param err - uncertainty
bool SFFitResults::GetFastDecTime(double &t, double &err){
  if(fFastDecTime==-1 || fFastDecTimeErr==-1){
    std::cerr << "##### Error in SFFitResults::GetFastDecTime()!" << std::endl;
    return false;
  }
  t   = fFastDecTime; 
  err = fFastDecTimeErr;
  return true;
}
//------------------------------------------------------------------
///Returns references to the slow decay time and its uncertainty.
///\param t - slow decay time [ns]
///\param err - uncertainty
bool SFFitResults::GetSlowDecTime(double &t, double &err){
  if(fSlowDecTime==-1 || fSlowDecTimeErr==1){
  std::cerr << "##### Error in SFFitResults::GetSlowDecTime()!" << std::endl;
  return false;
  }
  t   = fSlowDecTime;
  err = fSlowDecTimeErr;
  return true;
}
//------------------------------------------------------------------
///Returns references to the decay amplitude and its uncertainty
///(single decay mode).
///\param amp - amplitude 
///\param ampErr - uncertainty
bool SFFitResults::GetAmp(double &amp, double &ampErr){
  if(fAmp==-1 || fAmpErr==-1){
    std::cout << "##### Error in SFFitResults::GetAmp()!" << std::endl;
    return false;
  }
  amp    = fAmp;
  ampErr = fAmpErr;
}
//------------------------------------------------------------------
///Returns references to the amplitude of the fast decay component.
///\param ampFast - amplitude
///\param ampFastErr - uncertainty
bool SFFitResults::GetAmpFast(double &ampFast, double &ampFastErr){
  if(fAmpFast==-1 || fAmpFastErr==-1){
    std::cerr << "##### Error in SFFitResults::GetAmpFast()" << std::endl;
    return false;
  }
  ampFast    = fAmpFast;
  ampFastErr = fAmpFastErr;
  return true;
}
//------------------------------------------------------------------
///Returns references to the amplitude of the slow decay component. 
///\param ampSlow - amplitude
///\param ampSlowErr - uncertainty
bool SFFitResults::GetAmpSlow(double &ampSlow, double &ampSlowErr){
  if(fAmpSlow==-1 || fAmpSlowErr==-1){
    std::cerr << "##### Error in SFFitResults::GetAmpSlow()" << std::endl;
    return false;
  }
  ampSlow    = fAmpSlow;
  ampSlowErr = fAmpSlowErr;
  return true;
}
//------------------------------------------------------------------
///Returns references to the time offset and its uncertainty.
///\param t0 - time offset
///\param t0Err - uncertainty
bool SFFitResults::GetT0(double &t0, double &t0Err){
  if(fT0==-1 || fT0Err==-1){
    std::cerr << "##### Error in SFFitResults::GetT0()" << std::endl;
    return false;
  }
  t0    = fT0;
  t0Err = fT0Err;
  return true;
}
//------------------------------------------------------------------
///Returns references to the intensities of decay components.
///\param iSlow - intensity of slow component [%]
///\param iFast - intensity of fast component [%]
bool SFFitResults::GetIntensities(double &iSlow, double &iFast){
  if(fIslow==-1 || fIfast==-1){
    std::cerr << "##### Error in SFFitResults::GetIntensities()" << std::endl;
    return false;
  }
  iSlow = fIslow;
  iFast = fIfast;
  return true;
}
//------------------------------------------------------------------
///Returns references to the fitting range. 
///\param xmin - lower fitting range
///\param xmax - upper fitting range
bool SFFitResults::GetFitRange(double &xmin, double &xmax){
  if(fFitXmin==-1 || fFitXmax==-1){
    std::cerr << "##### Error in SFFitResults::GetFitRange()" << std::endl;
    return false;
  }
  xmin = fFitXmin;
  xmax = fFitXmax;
  return true;
}
//------------------------------------------------------------------
///Returns vector containing functions representing each component 
///of the decay process. Order in the vector (for double decay mode): 
///constant, fast decay, slow component or for single decay mode:
///constant, decay.
std::vector <TF1*> SFFitResults::GetCompFunctions(void){
  
 std::vector <TF1*> functions;

 TF1 *funBL = new TF1("funBL","pol0",0,1024);
 funBL->FixParameter(0,fConst);
 funBL->SetLineColor(kMagenta+3);
 
 TF1 *funDecay;
 TF1 *funFastDecay;
 TF1 *funSlowDecay;
 
 if(fComponents==1){
   funDecay = new TF1("funDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
   funDecay->FixParameter(0,fAmp);
   funDecay->FixParameter(1,fT0);
   funDecay->FixParameter(2,fDecTime);
   funDecay->SetLineColor(kMagenta);
 }
 else if(fComponents==2){
  funFastDecay = new TF1("funFastDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
  funFastDecay->FixParameter(0,fAmpFast);
  funFastDecay->FixParameter(1,fT0);
  funFastDecay->FixParameter(2,fFastDecTime);
  funFastDecay->SetLineColor(kMagenta);
 
  funSlowDecay = new TF1("funSlowDecay","[0]*(exp(-(x-[1])/[2]))",fFitXmin,1024);
  funSlowDecay->FixParameter(0,fAmpSlow);
  funSlowDecay->FixParameter(1,fT0);
  funSlowDecay->FixParameter(2,fSlowDecTime);
  funSlowDecay->SetLineColor(kMagenta-8);
 }

 functions.push_back(funBL);
 
 if(fComponents==1){
  functions.push_back(funDecay);
 }
 else if(fComponents==2){
  functions.push_back(funFastDecay);
  functions.push_back(funSlowDecay);
 }
 
 return functions;  
}
//------------------------------------------------------------------
///Returns pointer to the pave with results of the fit.
TPaveText *SFFitResults::GetResultsPave(void){
 
  TPaveText *pave = new TPaveText(0.456,0.436,0.903,0.901,"NBNDC");
  pave->SetTextFont(42);
  
  if(fComponents==1){
   pave->AddText(Form("A = %.2f +/- %.2f",fAmp,fAmpErr));
   pave->AddText(Form("t_0 = %.2f +/- %.2f",fT0,fT0Err));
   pave->AddText(Form("tau = %.2f +/- %.2f",fDecTime,fDecTimeErr));
   pave->AddText(Form("const = %.2f",fConst));
   pave->AddText(Form("Chi2/NDF = %.2f",fChi2/fNDF));
  }
  else if(fComponents==2){
   pave->AddText(Form("A_fast = %.2f +/- %.2f",fAmpFast,fAmpFastErr));
   pave->AddText(Form("t_0 = %.2f +/- %.2f",fT0,fT0Err));
   pave->AddText(Form("tau_fast = %.2f +/- %.2f",fFastDecTime,fFastDecTimeErr));
   pave->AddText(Form("A_slow = %.2f +/- %.2f",fAmpSlow,fAmpSlowErr));
   pave->AddText(Form("tau_slow = %.2f +/- %.2f",fSlowDecTime,fSlowDecTimeErr));
   pave->AddText(Form("const = %.2f",fConst));
   pave->AddText(Form("Chi2/NDF = %.2f",fChi2/fNDF));
   pave->AddText(Form("I_fast = %.2f perc.",fIfast));
   pave->AddText(Form("I_slow = %.2f perc.",fIslow));
  }
  
  return pave;
}
//------------------------------------------------------------------
///Prints details of the SFFitResults class object.
void SFFitResults::Print(void){
  std::cout << "\n------------------------------------------------------" << std::endl;
  std::cout << "This is Print of SFFitResults class object " << fName << std::endl;
  std::cout << "Fitted function: " << fFormula << std::endl;
  std::cout << "Fit status: " << fStat << std::endl;
  std::cout << "Number of parameters: " << fNpar << std::endl;
  std::cout << "Chi2 = " << fChi2 << std::endl;
  std::cout << "Chi2/NDF = " << fChi2/fNDF << std::endl;
  std::cout << "Fitting range: " << fFitXmin << " - " << fFitXmax << std::endl;
  if(fComponents==1){
    std::cout << "Decay time: " << fDecTime << " +/- " << fDecTimeErr << " ns" << std::endl;
  }
  else if(fComponents==2){
   std::cout << "Fast decay time: " << fFastDecTime << " +/- " << fFastDecTimeErr << " ns" << std::endl;
   std::cout << "Slow decay time: " << fSlowDecTime << " +/- " << fSlowDecTimeErr << " ns" << std::endl;
   std::cout << "Fast component intensity: " << fIfast << " %" << std::endl;
   std::cout << "Slow component intensity: " << fIslow << " %" << std::endl;
  }
  std::cout << "\nAll parameters: " << std::endl;
  for(int i=0; i<fNpar; i++){
   std::cout << std::setw(20);
   std::cout << "\t" << fParameters[i] << " +/- " << fParErrors[i] << std::endl;
  }
  std::cout << "------------------------------------------------------\n" << std::endl; 
  return;
}
//------------------------------------------------------------------
