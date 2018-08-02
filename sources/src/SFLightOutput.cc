// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFLightOutput.cc            *
// *             Jonas Kasper		   *
// *     kasper@physik.rwth-aachen.de      * 
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFLightOutput.hh"

ClassImp(SFLightOutput);

//------------------------------------------------------------------
///Default constructor.
SFLightOutput::SFLightOutput(){
  cout << "#### Warning in SFLightOutput constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  fSeriesNo = -1;
  fData = NULL;
  fAtt = NULL;
  Reset();
}
//------------------------------------------------------------------
///Standard constructor (recommended)
///\param seriesNo is number of experimental series to be analyzed. 
SFLightOutput::SFLightOutput(int seriesNo){
  Reset();
  fSeriesNo = seriesNo;
  fPDE=0.31; 
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFLightOutput constructor!";
  }
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    cout << "##### Warning in SFLightOutput constructor!";
    cout << "Calculating lightouput length with non-regular series!" << endl;
  }

  try{
    fAtt = new SFAttenuation(fSeriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFLightOutput constructor!";
  }
  
  //----- averaged channels method
  fAtt->AttAveragedCh();
  fAttLen   = fAtt->GetAttenuation();
  //----- separate channels method
  fAtt->AttSeparateCh(0);
  fSpectraCh0  = fAtt->GetSpectra(0);
  fPeaksCh0    = fAtt->GetPeaks(0);
  fAttLenCh0   = fAtt->GetAttenuation(0);
  
  fAtt->AttSeparateCh(1);
  fSpectraCh1 = fAtt->GetSpectra(1);
  fPeaksCh1   = fAtt->GetPeaks(1);
  fAttLenCh1   = fAtt->GetAttenuation(1);
  CalculateLO(0,0);
  CalculateLO(0,1);
  CalculateSLO(0);
  CalculateLO(1,0);
  CalculateLO(1,1);
  CalculateSLO(1);
}
//------------------------------------------------------------------
///Default destructor.
SFLightOutput::~SFLightOutput(){
  if(fData!=NULL) delete fData;
  if(fAtt!=NULL) delete fAtt;
}

//------------------------------------------------------------------
///Calculates the Ligthoutput for the single channels with the averaged attenuation length 
void SFLightOutput::CalculateLO(bool mode,int ch){

  cout << "\n----- Inside SFLightOutput::SeparateCh() for series " << fSeriesNo << endl;
  cout << "----- Analyzing channel " << ch << endl;
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  vector <TH1D*> spectra;
  vector <TH1D*> peaks;
  if(ch==0){
	spectra =fSpectraCh0;
  	peaks = fPeaksCh0;
  }
  else if(ch==1){
	spectra =fSpectraCh1;
  	peaks = fPeaksCh1;
  }
  TString gname;
  if(mode==0)gname = Form("LOAve_s%i_ch%i",fSeriesNo,ch);
  if(mode==1)gname = Form("LOSep_s%i_ch%i",fSeriesNo,ch);
  
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("LightOutput [P.E./MeV]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  vector <TF1*> fun;
  double mean=0;
  double sigma=0;
  double correctedLO=0;
  double correctedLOE=0;
  double dist=0;

  double lightout=0;
  double lightouterr=0;

  vector <double> temp_fat; 

  if(mode==0) temp_fat=fAttLen;
  if(mode==1){
	if(ch==0) temp_fat=fAttLenCh0;
	else if(ch==1) temp_fat=fAttLenCh1;
  }
  for(int i=0; i<npoints; i++){
    mean = peaks[i]->GetMean();
    sigma = peaks[i]->GetRMS();
    fun.push_back(new TF1("fun","gaus",mean-3*sigma,mean+3*sigma));
    peaks[i]->Fit(fun[i],"QR");
    if(ch==0) dist=positions[i];
    else if(ch==1) dist=100-positions[i];
    correctedLO= fun[i]->GetParameter(1)/fPDE/0.511/TMath::Exp(-dist/temp_fat[0]);
    correctedLOE= TMath::Sqrt((correctedLO*correctedLO/fun[i]->GetParameter(1)/fun[i]->GetParameter(1)*fun[i]->GetParameter(2)*fun[i]->GetParameter(2))+(correctedLO*correctedLO/temp_fat[0]/temp_fat[0]/temp_fat[0]/temp_fat[0]*temp_fat[1]*temp_fat[1]));
    graph->SetPoint(i,positions[i],correctedLO);
    graph->SetPointError(i,0,correctedLOE);
    lightout+=correctedLO*(1/correctedLOE/correctedLOE);
    lightouterr+=(1/correctedLOE/correctedLOE);
  }
   
  lightout=lightout/lightouterr;
  lightouterr=TMath::Sqrt(1/lightouterr);
  if(mode==0){
	if(ch==0){
		fLightOutAveGraphCh0=graph;
		fLightOutAveCh0=lightout;
		fLightOutAveErrCh0=lightouterr;
	}
	else if(ch==1){
		fLightOutAveGraphCh1=graph;
		fLightOutAveCh1=lightout;
		fLightOutAveErrCh1=lightouterr;
	}
  }
  else if(mode==1){
	if(ch==0){
		fLightOutSepGraphCh0=graph;
		fLightOutSepCh0=lightout;
		fLightOutSepErrCh0=lightouterr;
	}
	else if(ch==1){
		fLightOutSepGraphCh1=graph;
		fLightOutSepCh1=lightout;
		fLightOutSepErrCh1=lightouterr;
	}
  }
  cout << "\n      For channel " << ch << " is the Lightoutput: " << lightout << " +/- " << lightouterr << "P.E./MeV" << endl;
}

//------------------------------------------------------------------
///Calculates the Summed Lightouput
void SFLightOutput::CalculateSLO(bool mode){
  cout << "\n----- Inside SFLightOutput::Summed() for series " << fSeriesNo << endl;
  
   int npoints = fData->GetNpoints();
   vector <double> positions = fData->GetPositions();
   
   TString gname; 
   if(mode==0) gname = Form("LOAve_s%i",fSeriesNo);
   else if(mode==1) gname = Form("LOAve_s%i",fSeriesNo);
   TGraphErrors* graph = new TGraphErrors(npoints);
   graph->GetXaxis()->SetTitle("source position [mm]");
   graph->GetYaxis()->SetTitle("LightOutput [P.E./MeV]");
   graph->SetTitle(gname);
   graph->SetName(gname);
   graph->SetMarkerStyle(4);
   double x=0;
   double y=0;
   double SumLO=0;
   double SumLOE=0;
   double lightout=0;
   double lightouterr=0;

   for(int i=0; i<npoints; i++){
	   if(mode==0)fLightOutAveGraphCh0->GetPoint(i,x,y);
	   if(mode==1)fLightOutSepGraphCh0->GetPoint(i,x,y);
	   SumLO=y;
	   if(mode==0)fLightOutAveGraphCh1->GetPoint(i,x,y);
	   if(mode==1)fLightOutSepGraphCh1->GetPoint(i,x,y);
	   SumLO+=y;
	   if(mode==0)SumLOE=TMath::Sqrt((fLightOutAveGraphCh0->GetErrorY(i)*fLightOutAveGraphCh0->GetErrorY(i))+(fLightOutAveGraphCh1->GetErrorY(i)*fLightOutAveGraphCh1->GetErrorY(i)));
	   else if(mode==1)SumLOE=TMath::Sqrt((fLightOutSepGraphCh0->GetErrorY(i)*fLightOutSepGraphCh0->GetErrorY(i))+(fLightOutSepGraphCh1->GetErrorY(i)*fLightOutSepGraphCh1->GetErrorY(i)));
	   graph->SetPoint(i,positions[i],SumLO);
	   graph->SetPointError(i,0,SumLOE);

	   lightout+=SumLO*(1/SumLOE/SumLOE);
	   lightouterr+=(1/SumLOE/SumLOE);
   }
   lightout=lightout/lightouterr;
   lightouterr=TMath::Sqrt(1/lightouterr);

   if(mode==0){
	fLightOutAve=lightout;
	fLightOutAveErr=lightouterr;
	fLightOutAveGraph=graph;
   }
   if(mode==1){
	fLightOutSep=lightout;
	fLightOutSepErr=lightouterr;
	fLightOutSepGraph=graph;
   }
  cout << " \n      The summed Lightoutput is:" << lightout << " +/- " << lightouterr << "P.E./MeV" << endl;
}
//------------------------------------------------------------------
///Resets values of private members of the class to their default values.
void SFLightOutput::Reset(void){
 fLightOutAve= 0;
 fLightOutAveErr      = 0;
 fLightOutAveCh0= 0;
 fLightOutAveErrCh0      = 0;
 fLightOutAveCh1= 0;
 fLightOutAveErrCh1      = 0;
 fLightOutAveGraph    = NULL;
 fLightOutAveGraphCh0 = NULL;
 fLightOutAveGraphCh1 = NULL;
 fLightOutSep= 0;
 fLightOutSepErr      = 0;
 fLightOutSepCh0= 0;
 fLightOutSepErrCh0      = 0;
 fLightOutSepCh1= 0;
 fLightOutSepErrCh1      = 0;
 fLightOutSepGraph    = NULL;
 fLightOutSepGraphCh0 = NULL;
 fLightOutSepGraphCh1 = NULL;
 if(!fSpectraCh0.empty()) fSpectraCh0.clear();
 if(!fSpectraCh1.empty()) fSpectraCh1.clear();
 if(!fPeaksCh0.empty()) fPeaksCh0.clear();
 if(!fPeaksCh1.empty()) fPeaksCh1.clear();
 if(!fAttLen.empty()) fAttLen.clear();
 if(!fAttLenCh0.empty()) fAttLenCh1.clear();
 if(!fAttLenCh1.empty()) fAttLenCh0.clear();
}
//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFLightOutput::Print(void){
 cout << "\n-------------------------------------------" << endl;
 cout << "This is print out of SFLightOutput class object" << endl;
 cout << "Experimental series number " << fSeriesNo << endl;
 cout << "-------------------------------------------\n" << endl;
}
//------------------------------------------------------------------
TGraphErrors* SFLightOutput::GetLightOutputGraph(TString mode){
  if(mode.Contains("Averaged") | mode.Contains("Ave"))  return fLightOutAveGraph;
  else if(mode.Contains("Seperate") | mode.Contains("Sep"))  return fLightOutSepGraph;
}
//------------------------------------------------------------------
TGraphErrors* SFLightOutput::GetLightOutputGraph(TString mode, int ch){
  if(mode.Contains("Averaged") | mode.Contains("Ave")){ 
	if (ch==0)return fLightOutAveGraphCh0;
  	else if (ch==1)return fLightOutAveGraphCh1;
  }
  else if(mode.Contains("Seperate") | mode.Contains("Sep")){ 
	if (ch==0)return fLightOutSepGraphCh0;
  	else if (ch==1)return fLightOutSepGraphCh1;
  }
}
//------------------------------------------------------------------
vector<double> SFLightOutput::GetLightOutput(TString mode){
  vector<double> flo_return;
  if(mode.Contains("Averaged") | mode.Contains("Ave")){
	flo_return.push_back(fLightOutAve);
  	flo_return.push_back(fLightOutAveErr);
  }
  else if(mode.Contains("Seperate") | mode.Contains("Sep")){
	flo_return.push_back(fLightOutSep);
  	flo_return.push_back(fLightOutSepErr);
  }
  return flo_return;
}
//------------------------------------------------------------------
vector<double> SFLightOutput::GetLightOutput(TString mode,int ch){
  vector<double> flo_return;
  if(mode.Contains("Averaged") | mode.Contains("Ave")){
 	 if(ch==0){
 	       flo_return.push_back(fLightOutAveCh0);
 	 	flo_return.push_back(fLightOutAveErrCh0);
 	 }
 	 else if(ch==1){
 	       flo_return.push_back(fLightOutAveCh1);
 	 	flo_return.push_back(fLightOutAveErrCh1);
 	 }
  }
  else if(mode.Contains("Seperate") | mode.Contains("Sep")){
 	 if(ch==0){
 	       flo_return.push_back(fLightOutSepCh0);
 	 	flo_return.push_back(fLightOutSepErrCh0);
 	 }
 	 else if(ch==1){
 	       flo_return.push_back(fLightOutSepCh1);
 	 	flo_return.push_back(fLightOutSepErrCh1);
 	 }
  }
  return flo_return;
}
