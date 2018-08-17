// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimeConst.cc             *
// *       J. Kasper, K. Rusiecka          *
// *     kasper@physik.rwth-aachen.de      *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "SFTimeConst.hh"

ClassImp(SFTimeConst);

//------------------------------------------------------------------
/// Default constructor.
SFTimeConst::SFTimeConst(){
  cout << "##### Warning in SFTimeConst constructor! You are using default constructor!" << endl;
  cout << "Set object attributes via SetDetails()" << endl;
  Reset();
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
///\param seriesNo - number of analyzed series
///\param PE - PE value for signals selection
///\param verb - verbose level
SFTimeConst::SFTimeConst(int seriesNo, double PE, bool verb){
  bool stat = SetDetails(seriesNo,PE,verb);
  if(stat==false){
    throw "##### Exception in SFTimeConst constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
SFTimeConst::~SFTimeConst(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------
///Sets values to private members of the calss. Loads TProfile histograms
///of average signals for requested PE vlaue. Creates vectors of SFFitResults
///objects.
///\param seriesNo - number of the series
///\param PE - PE value
///\param verb - verbose level
bool SFTimeConst::SetDetails(int seriesNo, double PE, bool verb){
  
  fSeriesNo = seriesNo;
  fPE       = PE;
  fVerb     = verb;
  
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char* message){
    cout << message << endl;
    throw "##### Exception in SFTimeConst::SetDetails()!";
    return false;
  }
  
  int     npoints = fData->GetNpoints();
  TString fiber   = fData->GetFiber();
  vector <double> positions = fData->GetPositions();
  TString selection = Form("ch_0.fPE>%.1f && ch_0.fPE<%.1f",fPE-0.5,fPE+0.5);
  TString results_name;
  
  int nsig = 0;
  if(fiber.Contains("LuAG")) 
    nsig = 50;
  else if(fiber.Contains("LYSO")) 
    nsig = 1;
  else{
    cout << "##### Error in SFTimeConst::SetDetails()! Unknown fiber material!" << endl;
    return false;
  }
  
  for(int i=0; i<npoints; i++){
   fSignalsCh0.push_back(fData->GetSignalAverage(0,positions[i],selection,nsig,true));
   results_name = fSignalsCh0[i]->GetName();
   fResultsCh0.push_back(new SFFitResults(results_name));
   
   fSignalsCh1.push_back(fData->GetSignalAverage(1,positions[i],selection,nsig,true));
   results_name = fSignalsCh1[i]->GetName();
   fResultsCh1.push_back(new SFFitResults(results_name));
  }
  
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  
  return true;
}
//------------------------------------------------------------------
///Resets private members of the calss to their default values.
void SFTimeConst::Reset(void){
  fSeriesNo = -1;
  fPE       = -1;
  fVerb     = false;
  fData     = NULL;
  if(!fSignalsCh0.empty()) fSignalsCh0.clear();
  if(!fSignalsCh1.empty()) fSignalsCh1.clear();
  if(!fResultsCh0.empty()) fResultsCh0.clear();
  if(!fResultsCh1.empty()) fResultsCh1.clear();
}
//------------------------------------------------------------------
///Private method to get an index of requested measurements based on source position.
int SFTimeConst::GetIndex(double position){
 
  int index = -1;
  vector <double> positions = fData->GetPositions();
  int npoints = fData->GetNpoints();
  TString desc = fData->GetDescription();
  
  if(!desc.Contains("Regular series")){
    index = position-1;
    return index;
  }
  
  for(int i=0; i<npoints; i++){
    if(fabs(positions[i]-position)<1){
      index = i;
      break;
    }
  }
  
  if(index==-1){
   cout << "##### Error in SFTimingRes::GetIndex()! Incorrecct position!" << endl;
   return index;
  }
  
  return index;
}
//------------------------------------------------------------------
///Double decay function describing falling slope of the signal.
double funDecayDouble(double *x, double *par){
 double fast_dec = par[0]*TMath::Exp(-(x[0]-par[1])/par[2]); 
 double slow_dec = par[3]*TMath::Exp(-(x[0]-par[1])/par[4]); 
 double constant = par[5];
 return fast_dec+slow_dec+constant;
}
//------------------------------------------------------------------
double funDecaySingle(double *x, double *par){
 double dec = par[0]*TMath::Exp(-(x[0]-par[1])/par[2]); 
 double constant = par[3];
 return dec+constant;
}
//------------------------------------------------------------------
bool SFTimeConst::FitDecayTimeSingle(TProfile *signal, double position){
  
  TString opt;
  if(fVerb) opt = "R0";
  else opt = "QR0";
  
  int ch = -1;
  TString hname = signal->GetName();
  if(hname.Contains("ch0"))      ch = 0;
  else if(hname.Contains("ch1")) ch = 1;
  else{
    cout << "Error in SFTimeConst::FitDecayTimeSingle()!" << endl;
    cout << "Could not interpret signal name!" << endl;
    return false;
  }
  
  int index = GetIndex(position);
  double xmin = signal->GetBinCenter(signal->GetMaximumBin())+20.;
  double xmax = signal->GetBinCenter(signal->GetNbinsX());
  
  TF1 *fun_BL = new TF1("fun_BL","pol0",0,50);
  signal->Fit(fun_BL,opt);
  
  TF1 *fun_dec = new TF1("fun_dec","[0]*exp(-(x-[1])/[2])",xmin,xmin+100);
  fun_dec->SetParameters(100.,100.,10.);
  fun_dec->SetParNames("A","t0","tau");
  signal->Fit(fun_dec,opt);
  
  TF1 *fun_all = new TF1("fall",funDecaySingle,xmin,xmax,4);
  
  fun_all->SetParNames("A","t0","tau","const");
  fun_all->SetParameter(0,fun_dec->GetParameter(0));
  fun_all->SetParameter(1,fun_dec->GetParameter(1));
  fun_all->SetParameter(2,fun_dec->GetParameter(2));
  fun_all->FixParameter(3,fun_BL->GetParameter(0));
  int fitStat = signal->Fit(fun_all,opt);
    
  if(fitStat!=0){
    cout << "##### Warning in SFTimeConst::FitDecayTimeSingle()" << endl;
    cout << "\t fit status: " << fitStat << endl;
  }
  
  if(ch==0){
    fResultsCh0[index]->SetFromFunction(fun_all);
    if(fitStat!=0) fResultsCh0[index]->SetStat(-1);
    fResultsCh0[index]->Print();
  }
  else if(ch==1){
    fResultsCh1[index]->SetFromFunction(fun_all);
    if(fitStat!=0) fResultsCh1[index]->SetStat(-1);
    fResultsCh1[index]->Print();
  }
  
  return true;
}
//------------------------------------------------------------------
///This function performs fitting to the given TProfile signal. 
///Double decay function if fitted to the falling slope of the signal.
///Results of the fit are subsequently written in the SFFitResults 
///class object. Function returns true if fitting was successful and
///fit results are valid.
bool SFTimeConst::FitDecayTimeDouble(TProfile *signal, double position){
  
  TString opt;
  if(fVerb) opt = "R0";
  else opt = "QR0";
   
  int ch = -1;
  TString hname = signal->GetName();
  if(hname.Contains("ch0"))      ch = 0;
  else if(hname.Contains("ch1")) ch = 1;
  else{
    cout << "##### Error in SFTimeConst::FitDecayTimeDouble()!" << endl;
    cout << "Could not interpret signal name!" << endl;
    return false;
  }
  
  int index = GetIndex(position);
  double xmin = signal->GetBinCenter(signal->GetMaximumBin())+20;
  double xmax = signal->GetBinCenter(signal->GetNbinsX());

  TF1 *fun_BL = new TF1("fun_BL","pol0",0,50);
  signal->Fit(fun_BL,opt);
  
  TF1 *fun_fast = new TF1("fun_fast","[0]*exp(-(x-[1])/[2])",xmin,xmin+130);
  fun_fast->SetParameters(100.,100.,10.);
  fun_fast->SetParNames("A","t0","tau");
  signal->Fit(fun_fast,opt);

  TF1 *fun_slow = new TF1("fun_slow","[0]*exp(-(x-[1])/[2])",400,xmax);
  fun_slow->SetParameters(100.,100.,400.);
  fun_slow->SetParNames("A","t0","tau");
  fun_slow->FixParameter(1,fun_fast->GetParameter(1));
  signal->Fit(fun_slow,opt);

  TF1* fun_all = new TF1("fall",funDecayDouble,xmin,xmax,6);

  fun_all->SetParNames("A_fast","t0","tau_fast","A_slow","tau_slow","const");
  fun_all->SetParameter(0,fun_fast->GetParameter(0));
  fun_all->FixParameter(1,fun_fast->GetParameter(1));
  fun_all->SetParameter(2,fun_fast->GetParameter(2));
  fun_all->SetParameter(3,fun_slow->GetParameter(0));
  fun_all->SetParameter(4,fun_slow->GetParameter(2));
  fun_all->FixParameter(5,fun_BL->GetParameter(0));
  int fitStat = signal->Fit(fun_all,opt);
  
  if(fitStat!=0){
    cout << "##### Warning in SFTimeConst::FitDecayTimeDouble()" << endl;
    cout << "\t fit status: " << fitStat << endl;
  }
  
  if(ch==0){
    fResultsCh0[index]->SetFromFunction(fun_all);
    if(fitStat!=0) fResultsCh0[index]->SetStat(-1);
    fResultsCh0[index]->Print();
  }
  else if(ch==1){
    fResultsCh1[index]->SetFromFunction(fun_all);
    if(fitStat!=0) fResultsCh1[index]->SetStat(-1);
    fResultsCh1[index]->Print();
  }
  
  return true;
}
//------------------------------------------------------------------
///Fits all signals of the analyzed series from both channels.
///This function also calculates average time constants with their
///uncertainties.
bool SFTimeConst::FitAllSignals(void){
 
  int n = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  TString fiber = fData->GetFiber();
  
  if(fiber.Contains("LuAG")){
    for(int i=0; i<n; i++){
     FitDecayTimeDouble(fSignalsCh0[i],positions[i]);
     FitDecayTimeDouble(fSignalsCh1[i],positions[i]);
    }
  }
  else if(fiber.Contains("LYSO")){
    for(int i=0; i<n; i++){
     FitDecayTimeSingle(fSignalsCh0[i],positions[i]);
     FitDecayTimeSingle(fSignalsCh1[i],positions[i]);
    }
  }
  else{
    cout << "##### Error in SFTimeConst::FitAllSignals()!" << endl;
    cout << "Unknown fiber material!" << endl;
    return false;
  }
  
  //----- Calculating average time constants
  //----- and intensities
  double statCh0, statCh1;
  int counter = 0;
  
  //----- For double decay 
  vector <double> fastDecAll;
  vector <double> slowDecAll;
  double fastDec, fastDecErr;
  double slowDec, slowDecErr;
  double fastDecSum = 0;
  double slowDecSum = 0;
  
  vector <double> fastAmpAll;
  vector <double> slowAmpAll;
  double fastAmp, fastAmpErr;
  double slowAmp, slowAmpErr;
  double fastAmpSum = 0;
  double slowAmpSum = 0;
  
  if(fiber.Contains("LuAG")){
  
    for(int i=0; i<n; i++){
      statCh0 = fResultsCh0[i]->GetStat();
      statCh1 = fResultsCh1[i]->GetStat();
    
      if(statCh0==0){
        fResultsCh0[i]->GetFastDecTime(fastDec,fastDecErr);
        fResultsCh0[i]->GetSlowDecTime(slowDec,slowDecErr);
        fastDecAll.push_back(fastDec);
        slowDecAll.push_back(slowDec);
        fastDecSum += fastDec;
        slowDecSum += slowDec;
      
        fResultsCh0[i]->GetAmpFast(fastAmp,fastAmpErr);
        fResultsCh0[i]->GetAmpSlow(slowAmp,slowAmpErr);
        fastAmpAll.push_back(fastAmp);
        slowAmpAll.push_back(slowAmp);
        fastAmpSum += fastAmp;
        slowAmpSum += slowAmp;
      
        counter++;
      }
    
      if(statCh1==0){
        fResultsCh1[i]->GetFastDecTime(fastDec,fastDecErr);
        fResultsCh1[i]->GetSlowDecTime(slowDec,slowDecErr);
        fastDecAll.push_back(fastDec);
        slowDecAll.push_back(slowDec);
        fastDecSum += fastDec;
        slowDecSum += slowDec;
      
        fResultsCh1[i]->GetAmpFast(fastAmp,fastAmpErr);
        fResultsCh1[i]->GetAmpSlow(slowAmp,slowAmpErr);
        fastAmpAll.push_back(fastAmp);
        slowAmpAll.push_back(slowAmp);
        fastAmpSum += fastAmp;
        slowAmpSum += slowAmp;
      
        counter++;
      }
    }
  
    double fastDecAverage = fastDecSum/counter;
    double slowDecAverage = slowDecSum/counter;
    double fastDecAvErr = TMath::StdDev(&fastDecAll[0],&fastDecAll[counter-1])/sqrt(counter);
    double slowDecAvErr = TMath::StdDev(&slowDecAll[0],&slowDecAll[counter-1])/sqrt(counter);
  
    double fastAmpAverage = fastAmpSum/counter;
    double slowAmpAverage = slowAmpSum/counter;
    double fastAmpAvErr = TMath::StdDev(&fastAmpAll[0],&fastAmpAll[counter-1])/sqrt(counter);
    double slowAmpAvErr = TMath::StdDev(&slowAmpAll[0],&slowAmpAll[counter-1])/sqrt(counter);
  
    double denom = fastAmpAverage*fastDecAverage + slowAmpAverage*slowDecAverage;
    double Ifast = ((fastAmpAverage*fastDecAverage)/denom) * 100;
    double Islow = ((slowAmpAverage*slowDecAverage)/denom) * 100;
  
    cout << "\n\n----------------------------------" << endl;
    cout << "Average decay constants for whole series:" << endl;
    cout << "Fast decay: " << fastDecAverage << " +/- " << fastDecAvErr << " ns" << endl;
    cout << "Fast component intensity: " << Ifast << " %" << endl; 
    cout << "Slow decay: " << slowDecAverage << " +/- " << slowDecAvErr << " ns" << endl;
    cout << "Slow component intensity: " << Islow << " %" << endl; 
    cout << "Counter: " << counter << endl;
    cout << "----------------------------------" << endl;
  }
  
  //----- For single decay time
  vector <double> decAll;
  double dec, decErr;
  double decSum = 0;
  
  if(fiber.Contains("LYSO")){
    
    for(int i=0; i<n; i++){
      statCh0 = fResultsCh0[i]->GetStat();
      statCh1 = fResultsCh1[i]->GetStat();
      
      if(statCh0==0){
	fResultsCh0[i]->GetDecTime(dec,decErr);
	decAll.push_back(dec);
	decSum+=dec;
	counter++;
      }
      
      if(statCh1==0){
	fResultsCh1[i]->GetDecTime(dec,decErr);
	decAll.push_back(dec);
	decSum+=dec;
	counter++;
      }
    }
    
    double decAverage = decSum/counter;
    double decAvErr   = TMath::StdDev(&decAll[0],&decAll[counter-1])/sqrt(counter);
    
    cout << "\n\n----------------------------------" << endl;
    cout << "Average decay constant for the whole series:" << endl;
    cout << "Decay constant: " << decAverage << " ns" << endl;
    cout << "Counter: " << counter << endl;
    cout << "----------------------------------" << endl;
    
  }
  
  return true;
}
//------------------------------------------------------------------
///Fits all signals of the analyzed series from the chosen channel.
///\param ch - channel number.
bool SFTimeConst::FitAllSignals(int ch){
 
  int n = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  TString fiber = fData->GetFiber();
  
  if(ch!=0 || ch!=1){
    cout << "##### Error in SFTimeConst::FitAllSignals()!" << endl;
    cout << "Incorrect channel number. Possible options: 0 or 1" << endl;
    return false;
  }
  
  if(fiber.Contains("LuAG")){
    for(int i=0; i<n; i++){
     if(ch==0)       FitDecayTimeDouble(fSignalsCh0[i],positions[i]);
     else if (ch==1) FitDecayTimeDouble(fSignalsCh1[i],positions[i]);
    }
  }
  else if(fiber.Contains("LYSO")){
    for(int i=0; i<n; i++){
     if(ch==0)       FitDecayTimeSingle(fSignalsCh0[i],positions[i]);
     else if (ch==1) FitDecayTimeSingle(fSignalsCh1[i],positions[i]);
    }
  }
  else{
    cout << "##### Error in SFTimeConst::FitAllSignal()!" << endl;
    cout << "Unknown fiber material!" << endl;
    return false;
  }
  
  return true;
}
//------------------------------------------------------------------
///Fits single signal chosen based on the channel number
///and source position.
///\param ch - channel number
///\param position - source position in mm.
bool SFTimeConst::FitSingleSignal(int ch, double position){

  int index = GetIndex(position);
  vector <double> positions = fData->GetPositions();
  TString fiber = fData->GetFiber();
  
  if(ch!=0 || ch!=1){
   cout << "##### Error in SFTimeConst::FitSingleSignal()" << endl;
   cout << "Incorrect channel number. Possible options are: 0 or 1" << endl;
   return false;
  }
  
  if(fiber.Contains("LuAG")){
    if(ch==0)      FitDecayTimeDouble(fSignalsCh0[index],positions[index]);
    else if(ch==1) FitDecayTimeDouble(fSignalsCh1[index],positions[index]);
  }
  else if(fiber.Contains("LYSO")){
    if(ch==0)      FitDecayTimeSingle(fSignalsCh0[index],positions[index]);
    else if(ch==1) FitDecayTimeSingle(fSignalsCh1[index],positions[index]);
  }
  else{
    cout << "##### Error in SFTimeConst::FitSingleSignal()!" << endl;
    cout << "Unknown fiber material!" << endl;
    return false;
  }

  return true;
}
//------------------------------------------------------------------
///Returns vector containing signals for whole series and   
///chosen channel
///\param ch - channel number
vector <TProfile*> SFTimeConst::GetAllSignals(int ch){
  if(!(ch==0 || ch==1)){
    cout << "##### Error in SFTimeConst::GetAllSignals()!" << endl;
    cout << "Incorrect channel number!" << endl;
  }
  if(ch==0)      return fSignalsCh0;
  else if(ch==1) return fSignalsCh1;
}
//------------------------------------------------------------------
///Returns signal for requested channel and source position 
///\param ch - channel number
///\param position - position of the source in mm.
TProfile* SFTimeConst::GetSingleSignal(int ch, double position){
  int index = GetIndex(position);
  if      (ch==0) return fSignalsCh0[index];
  else if (ch==1) return fSignalsCh1[index];
  else{   
    cout << "##### Error in SFTimeConst::GetSingleSignal()!" << endl;
    cout << "Incorrect channel number!" << endl;
    return NULL;
  } 
}
//------------------------------------------------------------------
///Returns vector containing results of fits for whole series and  
///chosen channel
///\param ch - channel number
vector <SFFitResults*> SFTimeConst::GetAllResults(int ch){
  if(!(ch==0 || ch==1)){
    cout << "##### Error in SFTimeConst::GetAllResults()!" << endl;
    cout << "Incorrect channel number!" << endl;
  }
  if(ch==0)      return fResultsCh0;
  else if(ch==1) return fResultsCh1;
}
//------------------------------------------------------------------
///Returns fitting results for requested signal.
///\param ch - channel
///\param position - position of the source in mm
SFFitResults* SFTimeConst::GetSingleResult(int ch, double position){
  int index = GetIndex(position);
  if(ch==0)      return fResultsCh0[index];
  else if(ch==1) return fResultsCh1[index];
  else{
    cout << "##### Error in SFTimeConst::GetSingleResult()!" << endl;
    cout << "Incorrect channel number!" << endl;
    return NULL;
  }
}
//------------------------------------------------------------------
///Prints details of the SFTimeConst class ojbect.
void SFTimeConst::Print(void){
  cout << "\n-------------------------------------------" << endl;
  cout << "This is print out of SFTimeConst class object" << endl;
  cout << "Exparimental series number: " << fSeriesNo << endl;
  cout << "Signals PE: " << fPE << endl;
  cout << "Verbose level: " << fVerb << endl;
  cout << "-------------------------------------------\n" << endl;
}
//------------------------------------------------------------------