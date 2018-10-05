// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFEnergyResolution.cc         *
// *             Jonas Kasper		   *
// *     kasper@physik.rwth-aachen.de      * 
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFEnergyResolution.hh"

ClassImp(SFEnergyResolution);

//------------------------------------------------------------------
///Default constructor.
SFEnergyResolution::SFEnergyResolution(){
  cout << "#### Warning in SFEnergyResolution constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  fSeriesNo = -1;
  fData = NULL;
  Reset();
}
//------------------------------------------------------------------
///Standard constructor (recommended)
///\param seriesNo is number of experimental series to be analyzed. 
SFEnergyResolution::SFEnergyResolution(int seriesNo){
  Reset();
  fSeriesNo = seriesNo;
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFEnergyResolution constructor!";
  }
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    cout << "##### Warning in SFEnergyResolution constructor!";
    cout << "Calculating lightouput length with non-regular series!" << endl;
  }
  
  TString cut = Form("ch_%i.fT0>0 && ch_%i.fT0<590 && ch_%i.fPE>0",0,0,0);
  fSpectraCh0 = fData->GetSpectra(0,"fPE",cut);
  cut = Form("ch_%i.fT0>0 && ch_%i.fT0<590 && ch_%i.fPE>0",1,1,1);
  fSpectraCh1 = fData->GetSpectra(1,"fPE",cut);
  
  for (int i=0;i<fData->GetNpoints();i++){
	fPFCh0.push_back(new SFPeakFinder(fSpectraCh0[i],"511",false));
	fPFCh1.push_back(new SFPeakFinder(fSpectraCh1[i],"511",false));
  }
  CalculateER(0);
  CalculateER(1);

  try{
    fAtt = new SFAttenuation(fSeriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFEnergyResolution constructor!";
  }

  fAtt->AttAveragedCh();
  fAttLen   = fAtt->GetAttenuation();

  CalculateERAve();
}
//------------------------------------------------------------------
///Default destructor.
SFEnergyResolution::~SFEnergyResolution(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------
///Calculates the Ligthoutput for the single channels with the averaged attenuation length 
void SFEnergyResolution::CalculateER(int ch){

  cout << "\n----- Inside SFEnergyResolution::SeparateCh() for series " << fSeriesNo << endl;
  cout << "----- Analyzing channel " << ch << endl;
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  vector <SFPeakFinder*> tempPF;
  if(ch==0){
	tempPF = fPFCh0;
  }
  else if(ch==1){
	tempPF = fPFCh1;
  }
  TString gname = Form("ER_s%i_ch%i",fSeriesNo,ch);
  
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  vector<double> peak_par;
  double ER_single=0;
  double ER_singleerr=0;
  double dist=0;

  double ER=0;
  double ERerr=0;

  for(int i=0; i<npoints; i++){
    peak_par=tempPF[i]->GetParameter();
    if(ch==0) dist=positions[i];
    else if(ch==1) dist=100-positions[i];
    //peak_par[1]=peak_par[1]*2*TMath::Sqrt(2+TMath::Log(2));
    ER_single= peak_par[1]/peak_par[0];
    ER_singleerr=ER_single*TMath::Sqrt((peak_par[2]*peak_par[2]/peak_par[0]/peak_par[0])+(peak_par[3]*peak_par[3]/peak_par[1]/peak_par[1]));
    ER_single=ER_single*100;
    ER_singleerr=ER_singleerr*100;
    graph->SetPoint(i,positions[i],ER_single);
    graph->SetPointError(i,0,ER_singleerr);
    ER+=ER_single*(1/ER_singleerr/ER_singleerr);
    ERerr+=(1/ER_singleerr/ER_singleerr);
  }
   
  ER=ER/ERerr;
  ERerr=TMath::Sqrt(1/ERerr);
  if(ch==0){
		fEnergyResGraphCh0=graph;
		fEnergyResCh0=ER;
		fEnergyResErrCh0=ERerr;
  }
  else if(ch==1){
		fEnergyResGraphCh1=graph;
		fEnergyResCh1=ER;
		fEnergyResErrCh1=ERerr;
	}
  cout << "\n      For channel " << ch << " is the EnergyRes: " << ER << " +/- " << ERerr << "%" << endl;
}

void SFEnergyResolution::CalculateERAve(){

  cout << "\n----- Inside SFEnergyResolution::Averaged() for series " << fSeriesNo << endl;

  vector <SFPeakFinder*> tempPF;
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  double dist0,dist1;
 
  for (int i=0;i<fData->GetNpoints();i++){
    	dist0=positions[i];
    	dist1=100-positions[i];
	fSpectraCorCh0.push_back(fData->GetCustomHistogram(Form("AttLength ch_0.fPE/exp(%f/%f)",-dist0,fAttLen[0]),"",positions[i]));
	fSpectraCorCh1.push_back(fData->GetCustomHistogram(Form("AttLength ch_1.fPE/exp(%f/%f)",-dist1,fAttLen[0]),"",positions[i]));
 	fSpectraAve.push_back(fData->GetCustomHistogram(Form("AttLength ch_1.fPE/exp(%f/%f)+ch_1.fPE/exp(%f/%f)",-dist1,fAttLen[0],-dist0,fAttLen[0]),"",positions[i]));

  }
 
  for (int i=0;i<fData->GetNpoints();i++){
	tempPF.push_back(new SFPeakFinder(fSpectraAve[i],"511Sum",true));
  }
  TString gname = Form("ER_s%i_ave",fSeriesNo);
  
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  vector<double> peak_par;
  double ER_single=0;
  double ER_singleerr=0;

  double ER=0;
  double ERerr=0;

  for(int i=0; i<npoints; i++){
    peak_par=tempPF[i]->GetParameter();
    ER_single= peak_par[1]/peak_par[0];
    ER_singleerr=ER_single*TMath::Sqrt((peak_par[2]*peak_par[2]/peak_par[0]/peak_par[0])+(peak_par[3]*peak_par[3]/peak_par[1]/peak_par[1]));
    ER_single=ER_single*100;
    ER_singleerr=ER_singleerr*100;
    graph->SetPoint(i,positions[i],ER_single);
    graph->SetPointError(i,0,ER_singleerr);
    ER+=ER_single*(1/ER_singleerr/ER_singleerr);
    ERerr+=(1/ER_singleerr/ER_singleerr);
  }
   
  ER=ER/ERerr;
  ERerr=TMath::Sqrt(1/ERerr);
  fEnergyResGraphAve=graph;
  fEnergyResAve=ER;
  fEnergyResErrAve=ERerr;
  cout << "\n      Averaged is the EnergyRes: " << ER << " +/- " << ERerr << "%" << endl;
}
//------------------------------------------------------------------
///Resets values of private members of the class to their default values.
void SFEnergyResolution::Reset(void){
 if(!fSpectraCh0.empty()) fSpectraCh0.clear();
 if(!fSpectraCh1.empty()) fSpectraCh1.clear();
 if(!fPFCh0.empty())fPFCh0.clear();
 if(!fPFCh1.empty())fPFCh1.clear();
 fEnergyResCh0=0;	  
 fEnergyResErrCh0=0;
 fEnergyResCh1=0;	  
 fEnergyResErrCh1=0;
}
//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFEnergyResolution::Print(void){
 cout << "\n-------------------------------------------" << endl;
 cout << "This is print out of SFEnergyResolution class object" << endl;
 cout << "Experimental series number " << fSeriesNo << endl;
 cout << "-------------------------------------------\n" << endl;
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyResolution::GetEnergyResolutionGraph(){
  return fEnergyResGraphAve;
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyResolution::GetEnergyResolutionGraph(int ch){
	if (ch==0)return fEnergyResGraphCh0;
  	else if (ch==1)return fEnergyResGraphCh1;
}
//------------------------------------------------------------------
vector<TH1D*> SFEnergyResolution::GetSpectra(int ch){
  vector<TH1D*> temp_spec;
  if (ch==0) temp_spec=fSpectraCh0;
  else if (ch==1) temp_spec=fSpectraCh1;
  return temp_spec;
}
//------------------------------------------------------------------
vector<TH1D*> SFEnergyResolution::GetAveSpectra(){
  return fSpectraAve;
}
//------------------------------------------------------------------
vector<TH1D*> SFEnergyResolution::GetAveSpectra(int ch){
  vector<TH1D*> temp_spec;
  if (ch==0) temp_spec=fSpectraCorCh0;
  else if (ch==1) temp_spec=fSpectraCorCh1;
  return temp_spec;
}
//------------------------------------------------------------------
vector<double> SFEnergyResolution::GetEnergyResolution(int ch){
  vector<double> er_return;
  if(ch==0){
	  er_return.push_back(fEnergyResCh0);
	  er_return.push_back(fEnergyResErrCh0);
  }
  else if(ch==1){
	  er_return.push_back(fEnergyResCh1);
	  er_return.push_back(fEnergyResErrCh1);
  }
  return er_return;
}
//------------------------------------------------------------------
vector<double> SFEnergyResolution::GetEnergyResolution(){
  vector<double> er_return;
  er_return.push_back(fEnergyResAve);
  er_return.push_back(fEnergyResErrAve);
  return er_return;
}
