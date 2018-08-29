// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFAttenuation.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFAttenuation.hh"

ClassImp(SFAttenuation);

//------------------------------------------------------------------
///Default constructor.
SFAttenuation::SFAttenuation(){
  cout << "#### Warning in SFAttenuation constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  fSeriesNo = -1;
  fData = NULL;
  Reset();
}
//------------------------------------------------------------------
///Standard constructor (recommended)
///\param seriesNo is number of experimental series to be analyzed. 
SFAttenuation::SFAttenuation(int seriesNo){
  fSeriesNo = seriesNo;
  
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFAttenuation constructor!";
  }
  
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    cout << "##### Warning in SFAttenuation constructor!";
    cout << "Calculating attenuation length with non-regular series!" << endl;
  }
  
  Reset();
}
//------------------------------------------------------------------
///Default destructor.
SFAttenuation::~SFAttenuation(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------
///Method to determine attenuation length used in Pauwels et al., JINST 8 (2013) P09019.
///For both ends of the fiber one value is calculated, since averaged signal from both channels
///is taken into account.
bool SFAttenuation::AttAveragedCh(void){
 
  cout << "\n----- Inside SFAttenuation::AttAveragedCh() for series " << fSeriesNo << endl;
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  TString selection = "log(sqrt(ch_1.fPE/ch_0.fPE))";
  TString cut = "ch_0.fT0>0 && ch_0.fT0<590 && ch_1.fT0>0 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0";
  fRatios = fData->GetCustomHistograms(selection,cut);
  
  double mean, sigma;
  double fit_min, fit_max;
  vector <TF1*> fun;
  TString gname = Form("att_s%i",fSeriesNo);
  fAttnGraph = new TGraphErrors(npoints);
  fAttnGraph->GetXaxis()->SetTitle("source position [mm]");
  fAttnGraph->GetYaxis()->SetTitle("ln(M_{FB})");
  fAttnGraph->SetTitle(gname);
  fAttnGraph->SetName(gname);
  fAttnGraph->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    mean = fRatios[i]->GetMean();
    sigma = fRatios[i]->GetRMS();
    if(i<npoints/2){
      fit_min = mean-(1.5*sigma);
      fit_max = mean+(0.5*sigma);
    }
    else if(i==npoints/2 && npoints%2==1){
      fit_min = mean-sigma;
      fit_max = mean+sigma;
    }
    else{
      fit_min = mean-(0.5*sigma);
      fit_max = mean+(1.5*sigma);
    }
    fun.push_back(new TF1("fun","gaus",fit_min,fit_max));
    fRatios[i]->Fit(fun[i],"QR");
    fAttnGraph->SetPoint(i,positions[i],fun[i]->GetParameter(1));
    fAttnGraph->SetPointError(i,0,fun[i]->GetParError(1));
  }
  
  TF1 *fpol1 = new TF1("fpol1","pol1",positions[0],positions[npoints-1]);
  fAttnGraph->Fit(fpol1,"QR");
  fAttnLen = fabs(1./fpol1->GetParameter(1));
  fAttnErr = fpol1->GetParError(1)/pow(fpol1->GetParameter(1),2);
  
  cout << "\n\tAttenuation lenght is: " << fAttnLen << " +/- " << fAttnErr << " mm\n" << endl;
  
  return true;
}
//------------------------------------------------------------------
///Method to determine attenuation length for both channels independently. 
///\param ch - channel number
bool SFAttenuation::AttSeparateCh(int ch){
 
  cout << "\n----- Inside SFAttenuation::AttSeparateCh() for series " << fSeriesNo << endl;
  cout << "----- Analyzing channel " << ch << endl;
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  TString cut = Form("ch_%i.fT0>0 && ch_%i.fT0<590 && ch_%i.fPE>0",ch,ch,ch);
  vector <TH1D*> spectra = fData->GetSpectra(ch,"fPE",cut);
  
  TString gname = Form("att_s%i_ch%i",fSeriesNo,ch);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  vector <SFPeakFinder*> peakfin;
  vector <TH1D*> peaks;
  vector <double> parameter;
  
  for(int i=0; i<npoints; i++){
    peakfin.push_back(new SFPeakFinder(spectra[i],"511",false));
    peaks.push_back(peakfin[i]->GetPeak());
    parameter=peakfin[i]->GetParameter();
    graph->SetPoint(i,positions[i],parameter[0]);
    graph->SetPointError(i,0,parameter[2]);
  }
  
  TF1 *fexp = new TF1("fexp","expo",positions[0],positions[npoints-1]);
  graph->Fit(fexp,"QR");
  double attenuation = fabs(1./fexp->GetParameter(1));
  double att_error = fexp->GetParError(1)/pow(fexp->GetParameter(1),2);
 
  cout << "\n\tAttenuation for channel "<< ch << ": " << attenuation << " +/- " << att_error << " mm\n" << endl;
  
  if(ch==0){
    fAttnLenCh0 = attenuation;
    fAttnErrCh0 = att_error;
    fSpectraCh0 = spectra;
    fPeaksCh0 = peaks;
    fAttnGraphCh0 = graph;
  }
  else if(ch==1){
    fAttnLenCh1 = attenuation;
    fAttnErrCh1 = att_error;
    fSpectraCh1 = spectra;
    fPeaksCh1 = peaks;
    fAttnGraphCh1 = graph;
  }
  
  return true;
}
//------------------------------------------------------------------
///Returns vector containing histograms with signal ratios from both channels. 
///Histograms are used during averaged channels analysis i.e. in AttAveragedCh().
vector <TH1D*> SFAttenuation::GetRatios(void){
  if(fRatios.empty()){
    cout << "#### Error in SFAttenuation::GetRatios(). Empty vector!" << endl;
   }
   return fRatios;
}
//------------------------------------------------------------------
///Returns attenuation graph created in averaged channels method i.e. AttAveragedCh().
TGraphErrors* SFAttenuation::GetAttGraph(void){
  if(fAttnGraph==NULL){
    cout << "##### Error in SFAttenuation::GetAttGraph(). Empty pointer!" << endl;
    return NULL;
  }
  return fAttnGraph;
}
//------------------------------------------------------------------
///Returns attenuation length determined by averaged channels method i.e. AttAveragedCh().
double SFAttenuation::GetAttLength(void){
  if(fAttnLen==-1){
    cout << "##### Error in SFAttenuation::GetAttLength(). Incorrect value!" << endl;
    return 0;
  }
  return fAttnLen;
}
//------------------------------------------------------------------
///Returns error on attenuation length determined by averaged channels 
///method i.e. AttAveragedCh().
double SFAttenuation::GetAttError(void){
  if(fAttnErr==-1){
    cout << "##### Error in SFAttenuation::GetAttError(). Incorrect value!" << endl;
    return 0;
  }
  return fAttnErr;
}
//------------------------------------------------------------------
///Returns attenuation lenght and its error in form of a vector. Order in the vector:
///attenuation length, error. Both are in mm.
vector <double> SFAttenuation::GetAttenuation(void){
  vector <double> temp;
  if(fAttnLen==-1 || fAttnErr==-1){
   cout << "##### Error in SFAttenuation::GetAttenuation(). Incorrect attenuation length or error" << endl;
   cout << "\t" << fAttnLen << "\t" << fAttnErr << endl;
  }
  temp.push_back(fAttnLen);
  temp.push_back(fAttnErr);
  return temp;
}
//------------------------------------------------------------------
///Returns vector containing PE spectra used in determination of attenuation
///length with separate channels method i.e. AttSeparateCh().
///\param ch - channel number
vector <TH1D*> SFAttenuation::GetSpectra(int ch){
  if((ch==0 && fSpectraCh0.empty()) || (ch==1 && fSpectraCh1.empty())){
    cout << "##### Error in SFAttenuation::GetSpectra(). Empty vector!" << endl;
  }
  if(ch==0) return fSpectraCh0;
  else if(ch==1) return fSpectraCh1;
}
//------------------------------------------------------------------
///Returns vector containing peaks (spectra after background subtraction with 
///SFPeakFinder class) used in determination of attenuation length with separate
///channels method i.e. AttSeparateCh(). 
///\param ch - channel number.
vector <TH1D*> SFAttenuation::GetPeaks(int ch){
  if((ch==0 && fPeaksCh0.empty()) || (ch==1 && fPeaksCh1.empty())){
    cout << "##### Error in SFAttenuation::GetPeaks(). Empty vector!" << endl;
  }
  if(ch==0) return fPeaksCh0;
  else if(ch==1) return fPeaksCh1;
}
//------------------------------------------------------------------
///Returns attenuation graph for requested channel ch. Graph produced 
///with separate channels method - AttSeparateCh().
TGraphErrors* SFAttenuation::GetAttGraph(int ch){
   if((ch==0 && fAttnGraphCh0==NULL) || (ch==1 && fAttnGraphCh1==NULL)){
     cout << "##### Error in SFAttenuation::GetAttnGraph(int). Empty pointer!" << endl;
     return NULL;
   }
   if(ch==0) return fAttnGraphCh0;
   else if(ch==1) return fAttnGraphCh1;
}
//------------------------------------------------------------------
///Returns determined value of attenuation length and its error in a vector. 
///Order in a vector: attenuation length, error, both in mm. Values produced in
///AttSeparateCh().
///\param ch - channel number.
vector <double> SFAttenuation::GetAttenuation(int ch){
  vector <double> temp; 
  if(ch==0){
    if(fAttnLenCh0==-1 || fAttnErrCh0==-1){
      cout << "##### Error in SFAttenuation::GetAttenuation(int). Incorrect value!" << endl;
      cout << fAttnLenCh0 << "\t" << fAttnErrCh0 << endl;
    }
    temp.push_back(fAttnLenCh0);
    temp.push_back(fAttnErrCh0);
  }
  else if(ch==1){
    if(fAttnLenCh1==-1 || fAttnErrCh1==-1){
      cout << "##### Error in SFAttenuation::GetAttenuation(int). Incorrect value!" << endl;
      cout << fAttnLenCh1 << "\t" << fAttnErrCh1 << endl;
    }
    temp.push_back(fAttnLenCh1);
    temp.push_back(fAttnErrCh1);
  }
  return temp;
}
//------------------------------------------------------------------
///Returns attenuation length for requested channel determined by separate channels method
double SFAttenuation::GetAttLength(int ch){
  if((ch==0 && fAttnLenCh0==-1) || (ch==1 && fAttnLenCh1==-1)){
    cout << "##### Error in SFAttenuation:GetAttLength(int). Incorrect value!" << endl;
    return 0;
  }
  if(ch==0) return fAttnLenCh0;
  else if(ch==1) return fAttnLenCh1;
}
//------------------------------------------------------------------
///Returns error on attenuation length for requested channel determined by 
///separate channels method
double SFAttenuation::GetAttError(int ch){
  if((ch==0 && fAttnErrCh0==-1) || (ch==1 && fAttnErrCh1==-1)){
    cout << "##### Error in SFAttenuation:GetAttError(int). Incorrect value!" << endl;
    return 0;
  }
  if(ch==0) return fAttnErrCh0;
  else if(ch==1) return fAttnErrCh1;
}
//------------------------------------------------------------------
///Resets values of private members of the class to their default values.
void SFAttenuation::Reset(void){
 fAttnLen      = -1;
 fAttnErr      = -1;
 fAttnLenCh0   = -1;
 fAttnLenCh1   = -1;
 fAttnErrCh0   = -1;
 fAttnErrCh1   = -1;
 fAttnGraph    = NULL;
 fAttnGraphCh0 = NULL;
 fAttnGraphCh1 = NULL;
 if(!fRatios.empty()) fRatios.clear();
 if(!fSpectraCh0.empty()) fSpectraCh0.clear();
 if(!fSpectraCh1.empty()) fSpectraCh1.clear();
 if(!fPeaksCh0.empty()) fPeaksCh0.clear();
 if(!fPeaksCh1.empty()) fPeaksCh1.clear();
}
//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFAttenuation::Print(void){
 cout << "\n-------------------------------------------" << endl;
 cout << "This is print out of SFAttenuation class object" << endl;
 cout << "Experimental series number " << fSeriesNo << endl;
 cout << "-------------------------------------------\n" << endl;
}
//------------------------------------------------------------------
