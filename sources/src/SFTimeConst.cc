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
  bool flag = SetDetails(seriesNo,PE,verb);
  if(!flag){
    throw "##### Error in SFTimeConst constructor! Something went wrong in SetDetails()!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
SFTimeConst::~SFTimeConst(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------
///Sets values to private members of the calss. Loads TProfile histograms
///of average signals for requested PE vlaue. Creates vector of SFFitResults
///objects.
bool SFTimeConst::SetDetails(int seriesNo, double PE, bool verb){
  
  if(seriesNo<1 || seriesNo>9){
   cout << "##### Error in SFTimeConst::SetDetails()! Incorrect series number!" << endl;
   return false;
  }
  
  fSeriesNo = seriesNo;
  fPE       = PE;
  fVerb     = verb;
  fData     = new SFData(fSeriesNo);
  
  int npoints = fData->GetNpoints();
  double *positions = fData->GetPositions();
  TString selection = Form("ch_0.fPE>%.1f && ch_0.fPE<%.1f",fPE-0.5,fPE+0.5);
  TString results_name;
  
  for(int i=0; i<npoints; i++){
   fSignals.push_back(fData->GetSignalAverage(0,positions[i],selection, 100, true));
   results_name = fSignals[i]->GetName() + string("_single_decay");
   fResultsSingle.push_back(new SFFitResults(results_name));
   results_name = fSignals[i]->GetName() + string("_double_decay");
   fResultsDouble.push_back(new SFFitResults(results_name));
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
  if(!fSignals.empty()) fSignals.clear();
  if(!fResultsSingle.empty()) fResultsSingle.clear();
  if(!fResultsDouble.empty()) fResultsDouble.clear();
}
//------------------------------------------------------------------
///Private method to get an index of requested measurements based on source position.
int SFTimeConst::GetIndex(double position){
  int index = -1;
  double *positions = fData->GetPositions();
  int npoints = fData->GetNpoints();
  if(fSeriesNo>5){
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
bool SFTimeConst::FitDecTimeDouble(TProfile *signal, double position){
  
  TString opt;
  if(fVerb) opt = "R0";
  else opt = "QR0";
   
  int index = GetIndex(position);
  double xmin = signal->GetBinCenter(signal->GetMaximumBin())+20;
  double xmax = signal->GetBinCenter(signal->GetNbinsX());

  TF1 *fun_BL = new TF1("fun_BL","pol0",0,50);
  signal->Fit(fun_BL,opt);
  
  TF1 *fun_fast = new TF1("fun_fast","[0]*(exp(-(x-[1])/[2]))",xmin,xmin+50);
  fun_fast->SetParameters(100,100,50);
  fun_fast->SetParNames("A","t0","tau");
  signal->Fit(fun_fast,opt);

  TF1 *fun_slow = new TF1("fun_slow","[0]*(exp(-(x-[1])/[2]))",xmax-500,xmax);
  fun_slow->SetParameters(100,100,400);
  fun_slow->SetParNames("A","t0","tau");
  fun_slow->FixParameter(1,fun_fast->GetParameter(1));
  signal->Fit(fun_slow,opt);

  TF1* fun_all = new TF1("fall","[0]*(exp(-(x-[1])/[2])) + [3]*(exp(-(x-[1])/[4]))+[5]",xmin,xmax);

  fun_all->SetParNames("A_fast","t0","tau_fast","A_slow","tau_slow");
  fun_all->SetParameter(0,fun_fast->GetParameter(0));
  fun_all->FixParameter(1,fun_fast->GetParameter(1));
  fun_all->SetParameter(2,fun_fast->GetParameter(2));
  fun_all->SetParameter(3,fun_slow->GetParameter(0));
  fun_all->SetParameter(4,fun_slow->GetParameter(2));
  fun_all->FixParameter(5,fun_BL->GetParameter(0));
  
  fun_all->SetLineColor(kGreen+3);
  signal->Fit(fun_all,opt);
  fResultsDouble[index]->SetFromFunction(fun_all);
  fResultsDouble[index]->Print();
  
  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitDecTimeSingle(TProfile *signal, double position){

  TString opt;
  if(fVerb) opt = "R0";
  else opt = "QR0";
  
  int index = GetIndex(position);
  double xmin = signal->GetBinCenter(signal->GetMaximumBin())+20;
  double xmax = signal->GetBinCenter(signal->GetNbinsX());
  
  TF1 *fun_BL = new TF1("fun_BL","pol0",0,50);
  signal->Fit(fun_BL,opt);
  
  TF1 *fun_dec = new TF1("fun_dec","[0]*(exp(-(x-[1])/[2])) + [3]",xmin,xmax);
  fun_dec->SetParNames("A","t0","tau","const");
  fun_dec->SetParameter(0,100);
  fun_dec->SetParameter(1,100);
  fun_dec->SetParameter(2,100);
  fun_dec->FixParameter(3,fun_BL->GetParameter(0));
  
  signal->Fit(fun_dec,opt);
  fResultsSingle[index]->SetFromFunction(fun_dec);
  fResultsSingle[index]->Print();
  
  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitAllSignals(void){
 
  int n = fData->GetNpoints();
  double *positions = fData->GetPositions();
  
  for(int i=0; i<n; i++){
   FitDecTimeSingle(fSignals[i],positions[i]);
   FitDecTimeDouble(fSignals[i],positions[i]);
  }
  
  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitAllSignals(TString option){
  
  int n = fData->GetNpoints(); 
  double *positions = fData->GetPositions();
  
  for(int i=0; i<n; i++){
    if(option=="single")      FitDecTimeSingle(fSignals[i],positions[i]);
    else if(option=="double") FitDecTimeDouble(fSignals[i],positions[i]);
   }
  
  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitSingleSignal(double position){

  int index = GetIndex(position);
  double *positions = fData->GetPositions();
  
  FitDecTimeSingle(fSignals[index],positions[index]);
  FitDecTimeDouble(fSignals[index],positions[index]);
  
  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitSingleSignal(double position, TString option){

  int index = GetIndex(position);
  double *positions = fData->GetPositions();
  
  if(option=="single") FitDecTimeSingle(fSignals[index],positions[index]);
  if(option=="double") FitDecTimeDouble(fSignals[index],positions[index]);
  
  return true;
}
//------------------------------------------------------------------
///Returns signal for requested source position 
///\param position - position of the source in mm.
TProfile* SFTimeConst::GetSingleSignal(double position){
  int index = GetIndex(position);
  return fSignals[index];
}
//------------------------------------------------------------------
///Returns fitting results for requested signal.
///\param position - position of the source in mm
SFFitResults* SFTimeConst::GetSingleResult(double position, TString opt){
  int index = GetIndex(position);
  if(opt=="single")      return fResultsSingle[index];
  else if(opt=="double") return fResultsDouble[index];
  else{
    cout << "##### Error in SFTimeConst::GetSingleResult()! Incorrect option!" << endl;
    return NULL;
  }
}
//------------------------------------------------------------------
vector <SFFitResults*> SFTimeConst::GetAllResults(TString opt){
  if(opt!="single" || opt!="double"){
    cout << "##### Error in SFTimeConst::GetAllResults()! Incorrecct option!" << endl;
  }
  if(opt=="single")      return fResultsSingle;
  else if(opt=="double") return fResultsDouble;
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