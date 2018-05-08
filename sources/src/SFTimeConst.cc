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
   results_name = fSignals[i]->GetName();
   fResults.push_back(new SFFitResults(results_name));
  }

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  
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
  if(!fResults.empty()) fResults.clear();
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
bool SFTimeConst::FitDecTimeSingle(TProfile *signal, double position){

  TString opt;
  if(fVerb) opt = "0";
  else opt = "Q0";

  double tmax = signal->GetBinCenter(signal->GetMaximumBin());
  double xmin 	    = 0;
  double lastChiNDF = 0;
  double chiNDF     = 0;
  double results[12]; 
  
  TF1 *fun  = new TF1("exp_single","-[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]",0,1024);
  
  TF1 *pol0 = new TF1("pol0","pol0",0,50);
  signal->Fit(pol0,"Q0R");
  fun->FixParameter(4,pol0->GetParameter(0));
  
  for(int i=0; i<7; i++){
    fun->SetParameter(0,10);
    fun->SetParameter(1,10);
    fun->SetParameter(2,10);
    fun->SetParameter(3,100);
    xmin = tmax - (10*(i+3));
    signal->Fit(fun,opt,"",xmin,1024);
    chiNDF = fun->GetChisquare()/fun->GetNDF();
    if(lastChiNDF==0 || (fabs(1-lastChiNDF)>fabs(1-chiNDF))){
      lastChiNDF = chiNDF;
      results[0] = fun->GetParameter(0);	//parameter 0
      results[1] = fun->GetParameter(1);	//rise time [ns]
      results[2] = fun->GetParameter(2);	//parameter 2
      results[3] = fun->GetParameter(3);	//decay time [ns]
      results[4] = fun->GetParameter(4);	//parameter 4
      results[5] = fun->GetParError(0);		//parameter 0 error
      results[6] = fun->GetParError(1);		//rise time error [ns]
      results[7] = fun->GetParError(2);		//parameter 2 error
      results[8] = fun->GetParError(3);		//decay time error [ns]
      results[9] = fun->GetParError(4);		//para,eter 4 error
      results[10] = xmin;			//lower fitting range
      results[11] = chiNDF;			//Chi2/NDF
    }
  }
  
  int index = GetIndex(position);
  int npar = 5;
  vector <double> par(npar);
  vector <double> err(npar);
  
  for(int i=0; i<npar; i++){
    par[i] = results[i];
    err[i] = results[i+npar];
  }
  
  TF1 *decayFunction = new TF1("single_decay","-[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]",results[10],1024);
  
  for(int i=0; i<npar; i++){
    decayFunction->FixParameter(i,results[i]);
  }
 
  //fResults[index]->SetDecayTime(results[1],results[4]);  
  //fResults[index]->SetDecayChi2NDF(results[7]);
  //fResults[index]->SetRangesDecay(results[6],1024);
  fResults[index]->SetFormulaDecay(fun->GetExpFormula());
  fResults[index]->SetDecayPar(par);
  fResults[index]->SetDecayParErrors(err);
  fResults[index]->SetDecayFun(decayFunction);
    
  cout << "\n-----------------------" << endl;
  cout << signal->GetName() << endl;
  cout << "Rise time: " << results[1] << " +/- " << results[6] << " ns" << endl;
  cout << "Decay time: " << results[3] << " +/- " << results[8] << " ns" << endl;
  cout << "Chi2/NDF: " << results[11] << endl;
  cout << "Tmax = " << tmax << "\t Xmin = " << results[10] << endl;
  cout << "-----------------------\n" << endl; 

  return true;
}
//------------------------------------------------------------------
bool SFTimeConst::FitDecTimeDouble(TProfile *signal, double position){
 
  TString opt;
  if(fVerb) opt = "0";
  else opt = "Q0";

  double tmax = signal->GetBinCenter(signal->GetMaximumBin());
  double xmin  = 0;
  double split = 0;
  double lastChiNDF = 0;
  double chiNDF;
  double results[16];
  
  TString formula = "-[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]*TMath::Exp(-x/[5])+[6]";
  //TString formula = "[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])+[4]"
  TF1 *fun = new TF1("exp_double",formula,0,1024);
  
  TF1 *fastExp = new TF1("fastExp","[0]*TMath::Exp(-x/[1])+[2]",0,1024);
  TF1 *slowExp = new TF1("slowExp","[0]*TMath::Exp(-x/[1])+[2]",0,1024);
  
  TF1 *pol0 = new TF1("pol0","pol0",0,50); 
  signal->Fit(pol0,"R0Q");
  fun->FixParameter(6,pol0->GetParameter(0));
  
  //for(int i=0; i<7; i++){
  while(xmin<tmax-50){
    fun->SetParameter(0,10);
    fun->SetParameter(1,10);
    fun->SetParameter(2,10);
    fun->SetParameter(3,100);
    fun->SetParameter(4,10);
    fun->SetParameter(5,300);
    xmin = xmin + (10);
    signal->Fit(fun,opt,"",xmin,1024);
    chiNDF = fun->GetChisquare()/fun->GetNDF();
    if(lastChiNDF==0 || (fabs(1-lastChiNDF)>fabs(1-chiNDF))){
      lastChiNDF = chiNDF;
      results[0]  = fun->GetParameter(0);	//parameter 0
      results[1]  = fun->GetParameter(1); 	//rise time [ns]
      results[2]  = fun->GetParameter(2);	//parameter 2
      results[3]  = fun->GetParameter(3);	//fast decay time [ns]
      results[4]  = fun->GetParameter(4);	//parameter 4
      results[5]  = fun->GetParameter(5);	//slow decay time [ns]
      results[6]  = fun->GetParameter(6);	//parameter 4
      results[7]  = fun->GetParError(0);
      results[8]  = fun->GetParError(1);
      results[9]  = fun->GetParError(2);
      results[10] = fun->GetParError(3);
      results[11] = fun->GetParError(4);
      results[12] = fun->GetParError(5);
      results[13] = fun->GetParError(6);
      results[14] = xmin;
      results[15] = chiNDF;
    }
  }
  
  /*
  for(int i=0; i<5; i++){
    for(int ii=0; ii<5; ii++){
      for(int bin=signal->GetMaximumBin(); bin<signal->GetNbinsX(); bin++){
	if(signal->GetBinContent(bin) < (signal->GetMaximum()/(ii+2))){
	  split = signal->GetBinCenter(bin);
	  break;
	}
      }
      xmin = tmax + (i*10);
      fastExp->SetParameters(1,10,1);
      slowExp->SetParameters(1,10,1);
      signal->Fit(fastExp,opt,"",xmin,split);
      signal->Fit(slowExp,opt,"",split,1024);
      fun->SetParameter(0,fastExp->GetParameter(0));
      fun->SetParameter(1,fastExp->GetParameter(1));
      fun->SetParameter(2,slowExp->GetParameter(0));
      fun->SetParameter(3,slowExp->GetParameter(1));
      signal->Fit(fun,opt,"",xmin,1024);
      chiNDF = fun->GetChisquare()/fun->GetNDF();
      if(lastChiNDF==0 || (fabs(1-lastChiNDF) < fabs(1-chiNDF))){
	lastChiNDF = chiNDF;
	results[0]  = fun->GetParameter(0);	//parameter 0
	results[1]  = fun->GetParameter(1); 	//fast decay time [ns]
	results[2]  = fun->GetParameter(2);	//parameter 2
	results[3]  = fun->GetParameter(3);	//slow decay time [ns]
	results[4]  = fun->GetParameter(4);	//parameter 4
	results[5]  = fun->GetParError(0);	//parameter 0 error
	results[6]  = fun->GetParError(1);	//fast decay time error [ns]
	results[7]  = fun->GetParError(2);	//parameter 2 error 
	results[8]  = fun->GetParError(3);	//slow decay time error [ns]
	results[9]  = fun->GetParError(4);	//parameter 4 error 
	results[10] = xmin;			//lower fitting range
	results[11] = split;			//fast/slow split time
	results[12] = chiNDF;			//Chi2/NDF
      }
    }
  }
 */
  
  int index = GetIndex(position);
  int npar = 5; 
  vector <double> par(npar);
  vector <double> err(npar);
  
  for(int i=0; i<npar; i++){
    par[i] = results[i];
    err[i] = results[i+npar];
  }
  
  TF1* decayFunction = new TF1("double_decay",formula,results[14],1024);
  for(int i=0; i<npar; i++){
  decayFunction->FixParameter(i,results[i]);
  }
  
  double denom = results[2]*results[3] + results[4]*results[5];
  double Ifast = ((results[2]*results[3])/denom)*100;
  double Islow = ((results[4]*results[5])/denom)*100;


  //fResults[index]->SetFastDecTime(results[1],results[6]); 
  //fResults[index]->SetSlowDecTime(results[3],results[8]);
  fResults[index]->SetIntensities(Islow,Ifast);
  //fResults[index]->SetSplitTime(results[11]);
  //fResults[index]->SetRangesDecay(results[10],1024);
  //fResults[index]->SetDecayChi2NDF(results[12]);
  fResults[index]->SetDecayPar(par);
  fResults[index]->SetDecayParErrors(err);
  fResults[index]->SetDecayFun(decayFunction);
  fResults[index]->SetFormulaDecay(fun->GetExpFormula());
  
  cout << "\n-----------------------" << endl;
  cout << signal->GetName() << endl;
  cout << "Rise time: " << results[1] << " +/- " << results[7] << " ns" << endl;
  cout << "Fast decay time: " << results[3] << " +/- " << results[10] << " ns" << endl;
  cout << "Slow decay time: " << results[5] << " +/- " << results[12] << " ns" << endl;
  cout << "Fast component intensity: " << Form("%.2f",Ifast) << " %" << endl;
  cout << "Slow component intensity: " << Form("%.2f",Islow) << " %" << endl;
  cout << "Chi2/NDF: " << results[15] << endl;
  cout << "-----------------------\n" << endl; 
  
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
SFFitResults* SFTimeConst::GetSingleResult(double position){
  int index = GetIndex(position);
  return fResults[index];
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