// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimingRes.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFTimingRes.hh"

ClassImp(SFTimingRes);

//------------------------------------------------------------------
/// Default constructor.
SFTimingRes::SFTimingRes(){
  cout << "##### Warning in SFTimingRes constructor!" << endl;
  cout << "You are using default constructor!" <<endl;
  Reset();
}
//------------------------------------------------------------------
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed.
/// \param threshold identifies which tree should be opened. Possible
/// options are: ft - fixed threshold and cf - constant fraction.
/// \param method is a flag to choose type of analysis. Possible options are:
/// with cut - for analysis with energy cut on 511 keV, no cut - cutting only
/// scattered events.
SFTimingRes::SFTimingRes(int seriesNo, TString threshold, TString method){
  
  Reset();
  
  if(!(threshold=="ft" || threshold=="cf")){
    throw "##### Error in SFTimingRes constructor! Incorrect threshold type!\nPossible options are: ft, cf";
  }
  
  if(!(method=="with cut" || method=="no cut")){
    throw "##### Error in SFTimingRes constructor! Incorrect method!\nPossible options are: with cut, no cut";
  }
  
  fSeriesNo = seriesNo;
  fThreshold = threshold;
  fMethod = method;
  
  try{
    fData = new SFData(fSeriesNo,threshold);
  }
  catch(const char *message){
    cout << message << endl;
    throw "##### Exception in SFTimingRes constructor!";
  }
  fType = fData->GetMeasureType();

  
  if(fMethod=="no cut")        AnalyzeNoECut();
  else if(fMethod=="with cut") AnalyzeWithECut();
}
//------------------------------------------------------------------
/// Default destructor.
SFTimingRes::~SFTimingRes(){
 if(fData!=NULL) delete fData; 
}
//------------------------------------------------------------------
///Private method to get index of requested measurement based on source position.
int SFTimingRes::GetIndex(double position){
  
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
///Private method to get ratio histograms necesarry to impose cuts. 
///Depending on measurement type single or double gaussian function
///is fitted to the histograms.
bool SFTimingRes::LoadRatios(void){
  
  int npoints = fData->GetNpoints();
  
  TString selection = "log(sqrt(ch_1.fPE/ch_0.fPE))";
  TString cut = "ch_0.fT0<590 && ch_0.fPE>0 && ch_0.fT0>0 && ch_1.fT0<590 && ch_1.fT0>0 && ch_1.fPE>0"; 
  fRatios = fData->GetCustomHistograms(selection,cut);
  
  vector <TF1*> fun;
  double min, max;
  
  for(int i=0; i<npoints; i++){
    if(fType.Contains("Lead")){
      fun.push_back(new TF1("fun","gaus(0)+gaus(3)",-1,1));
      fun[i]->SetParameter(0,fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));		//thin gauss
      fun[i]->SetParameter(1,fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin()));
      fun[i]->SetParameter(2,6E-2);
      fun[i]->SetParameter(3,0.5*fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));	//thick gauss
      if(i<npoints/2)
        fun[i]->SetParameter(4,fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())+0.2);
      else
	fun[i]->SetParameter(4,fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())-0.2);
      fun[i]->SetParameter(5,2E-1);
      fRatios[i]->Fit(fun[i],"QR");
    }
    else if(fType.Contains("Electric")){
      min = fRatios[i]->GetMean()-2*fRatios[i]->GetRMS();
      max = fRatios[i]->GetMean()+2*fRatios[i]->GetRMS();
      fun.push_back(new TF1("fun","gaus",min,max));		//single gauss
      fRatios[i]->Fit(fun[i],"QR");
    }
  }

  return true;
}
//------------------------------------------------------------------
///Formula for Lorentzian function (if needed).
double LorentzianFun(double *x, double *par){
 //par0 - integral
 //par1 - width (FWHM)
 //par2 - mean  
 return (0.5*par[0]*par[1]/TMath::Pi())/
        ((x[0]-par[2])*(x[0]-par[2])+0.25*par[1]*par[1]);
}
//------------------------------------------------------------------
///Method for simple timing resolution determination. Timing resolution spectrum
///is build as difference between channel 0 and 1 signals with additional cut
///on scattered events. Timing resolution and its uncertainty is determined 
///based on Gaussian function as its sigma.
///If measurement was taken with lead collimator - double Gauss fitted.
///If measurement was taken with electronic collimatir - single Gaussian.
bool SFTimingRes::AnalyzeNoECut(void){
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  
  TString selection = "ch_0.fT0-ch_1.fT0";
  TString cut;
  double mean, sigma;
  
  if(fRatios.empty()) LoadRatios();
  vector <TF1*> fun;

  fT0Graph = new TGraphErrors(npoints);
  fT0Graph->GetXaxis()->SetTitle("source position [mm]");
  fT0Graph->GetYaxis()->SetTitle("mean of time difference distribution [ns]");
  fT0Graph->SetTitle("");
  fT0Graph->SetMarkerStyle(8);
  
  for(int i=0; i<npoints; i++){
    
    mean = fRatios[i]->GetFunction("fun")->GetParameter(1);
    sigma = fRatios[i]->GetFunction("fun")->GetParameter(2);
    cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",mean-0.5*sigma, mean+0.5*sigma);
    fT0Diff.push_back(fData->GetCustomHistogram(selection,cut,positions[i]));
    
    if(fType.Contains("Lead")){
      fun.push_back(new TF1("fun","gaus(0)+gaus(3)",-100,100));
      fun[i]->SetParameter(0,fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1,fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2,fT0Diff[i]->GetRMS());
      fun[i]->SetParameter(3,fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin())/4);
      if(i<npoints/2)
	fun[i]->SetParameter(4,fT0Diff[i]->GetMean()-5);
      else
	fun[i]->SetParameter(4,fT0Diff[i]->GetMean()+5);
      fun[i]->SetParameter(5,fT0Diff[i]->GetRMS()*2);
      fT0Diff[i]->Fit(fun[i],"QR");
    }
    else if(fType.Contains("Electric")){
      fT0Diff.back()->Rebin(2);
      fun.push_back(new TF1("fun","gaus",-50,50));
      fun[i]->SetParameter(0,fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1,fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2,fT0Diff[i]->GetRMS());
      fT0Diff[i]->Fit(fun[i],"QR");
    }

    fTimeRes.push_back(fabs(fun[i]->GetParameter(2)));
    fTimeResErr.push_back(fun[i]->GetParError(2));
    fT0Graph->SetPoint(i,positions[i],fun[i]->GetParameter(1));
    fT0Graph->SetPointError(i,2,fabs(fun[i]->GetParameter(2)));	//position uncertainty 2mm
  }

  return true;
}
//------------------------------------------------------------------
///Method for timing resolution determination with an energy cut on 511 keV peak.
///Timing resolution based on the Gaussian fit to the histogram - sigma of the fitted
///function. Regardless the measurement type - always sigle Gauss fitted.
bool SFTimingRes::AnalyzeWithECut(void){
  
  int npoints = fData->GetNpoints();
  vector <double> positions = fData->GetPositions();
  
  if(fRatios.empty()) LoadRatios();
  
  fPEch0 = fData->GetSpectra(0,"fPE","ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  fPEch1 = fData->GetSpectra(1,"fPE","ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0");
  vector <SFPeakFinder*> peakFin_ch0;
  vector <SFPeakFinder*> peakFin_ch1;
  TF1* fun = new TF1("fun","gaus",-200,200);
  double xmin_ch0, xmax_ch0;
  double xmin_ch1, xmax_ch1;
  double mean_ratio, sigma_ratio;
  double mean, sigma;
  //double f = 2*sqrt(2*log(2));		//to recalculate sigma into FWHM
  TString selection = "ch_0.fT0-ch_1.fT0";
  TString cut;
  
  double center_ch0, delta_ch0;		//changed here for smaller cut
  double center_ch1, delta_ch1;		//
  
  fT0Graph = new TGraphErrors(npoints);
  fT0Graph->GetXaxis()->SetTitle("source position [mm]");
  fT0Graph->GetYaxis()->SetTitle("mean of time difference distribution [ns]");
  fT0Graph->SetTitle("");
  fT0Graph->SetMarkerStyle(8);
  
  for(int i=0; i<npoints; i++){
    peakFin_ch0.push_back(new SFPeakFinder(fPEch0[i],"511",false));
    peakFin_ch1.push_back(new SFPeakFinder(fPEch1[i],"511",false));
    peakFin_ch0[i]->FindPeakRange(xmin_ch0,xmax_ch0);
    peakFin_ch1[i]->FindPeakRange(xmin_ch1,xmax_ch1);
    center_ch0 = xmin_ch0+(xmax_ch0-xmin_ch0)/2.;	//changed here for smaller cut
    delta_ch0  = (xmax_ch0-xmin_ch0)/6.;		//
    center_ch1 = xmin_ch1+(xmax_ch1-xmin_ch1)/2.;	//
    delta_ch1  = (xmax_ch1-xmin_ch1)/6.;		//
    mean_ratio = fRatios[i]->GetFunction("fun")->GetParameter(1);
    sigma_ratio = fRatios[i]->GetFunction("fun")->GetParameter(2);
    cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", center_ch0-delta_ch0,center_ch0+delta_ch0,center_ch1-delta_ch1,center_ch1+delta_ch1,mean_ratio-0.5*sigma_ratio, mean_ratio+0.5*sigma_ratio);	//changed here for smaller cut
    fT0Diff.push_back(fData->GetCustomHistogram(selection,cut,positions[i]));
    mean = fT0Diff[i]->GetMean();
    sigma = fT0Diff[i]->GetRMS();
    fT0Diff[i]->Fit(fun,"Q","",mean-5*sigma,mean+5*sigma);
    fTimeRes.push_back(fun->GetParameter(2));	//Timing resolution as sigma, if FWHM needed multiply by f
    fTimeResErr.push_back(fun->GetParError(2));		//FWHM - multiply by f
    fT0Graph->SetPoint(i,positions[i],fun->GetParameter(1));
    fT0Graph->SetPointError(i,2,fun->GetParameter(2));	//position uncertainy 2 mm
							//FWHM - multiply by f
  }
  
  return true;
}
//------------------------------------------------------------------
///Returns vector of histograms, which contains timing resolution spectra
///for whole series.
vector <TH1D*> SFTimingRes::GetT0Diff(void){
  if(fT0Diff.empty()){
   cout << "##### Error in SFTimingRes::GetT0Diff()! Empty vector!" << endl; 
  }
  return fT0Diff;
}
//------------------------------------------------------------------
///Returns single timing resolution spectrum for requested measurement.
///Measurement is identified by source position. If source position is not unique
///number of measurement in the series should be given.
TH1D* SFTimingRes::GetT0Diff(double position){
  int index = GetIndex(position);
  if(index>fT0Diff.size()){
    cout << "##### Error in SFTimingRes::GetT0Diff(position)! Incorrect position!" << endl;
    return NULL;
  }
  return fT0Diff[index];
}
//------------------------------------------------------------------
///Returns graph ch_0.T0-ch_1.T0 vs. source position.
TGraphErrors* SFTimingRes::GetT0Graph(void){
  if(fT0Graph==NULL){
    cout << "##### Error in SFTimingRes::GetT0Graph()! Graph was not created!" << endl;
    return NULL;
  }
  return fT0Graph;
}
//------------------------------------------------------------------
///Returns timing resolution value and its uncertainty for requested measurement.
///Measurement is identified like in GetT0Diff() function. Order in the returned vector:
///timing resolution, uncertainty.
vector <double> SFTimingRes::GetTimingResolution(double position){
  
  vector <double> temp;
  int index = GetIndex(position);
  
  if(index<fTimeRes.size()){
    cout << "##### Error in SFTimingRes::GetTimingResolution(position)! Incorrect position!" << endl;
  }
  
  temp.push_back(fTimeRes[index]);
  temp.push_back(fTimeResErr[index]);
  
  return temp;
}
//------------------------------------------------------------------
///Returns vector containing timing resolutions for all measurements 
///in the series
vector <double> SFTimingRes::GetTimingResolutions(void){
  if(fTimeRes.empty()){
    cout << "##### Error in SFTimingRes::GetTimingResolutions()! Empty vector!" << endl;
  }
  return fTimeRes;
}
//------------------------------------------------------------------
///Returns vector containing uncertainties of timing resolutions for 
///all measurements in the series.
vector <double> SFTimingRes::GetTimingResErrors(void){
  if(fTimeResErr.empty()){
    cout << "##### Error in SFTimingRes::GetTimingResErrors()! Empty vector!" << endl;
  }
  return fTimeResErr;
}
//------------------------------------------------------------------
/// Returns vector of ratio histograms used for cut on scattered events.
vector <TH1D*> SFTimingRes::GetRatios(void){
 if(fRatios.empty()){
   cout << "##### Error in SFTimingRes::GetRatios()! Empty vector!" << endl;
 }
 return fRatios; 
}
//------------------------------------------------------------------
/// Returns vector of PE spectra used for 511 keV energy cut. 
/// \param ch - channel number (0 or 1).
vector <TH1D*> SFTimingRes::GetSpectra(int ch){
 if(fMethod=="no cut"){
   cout << "##### Error in SFTimingRes::GetSpectra()! Incorrect method!" << endl;
 }
 if((ch==0 && fPEch0.empty()) || (ch==1 && fPEch1.empty())){
   cout << "##### Error in SFTimingRes::GetSpectra()! Empty vector!" << endl;
 }
 if(ch==0)       return fPEch0; 
 else if (ch==1) return fPEch1;
}
//------------------------------------------------------------------
///Resets private members of the class to their default values.
void SFTimingRes::Reset(void){
  fSeriesNo  = -1;
  fThreshold = "dummy";
  fMethod    = "dummy";
  fData      = NULL;
  fT0Graph   = NULL;
  if(!fRatios.empty())      fRatios.clear();
  if(!fPEch0.empty())       fPEch0.clear();
  if(!fPEch1.empty())       fPEch1.clear();
  if(!fT0Diff.empty())      fT0Diff.clear();
  if(!fTimeRes.empty())     fTimeRes.clear();
  if(!fTimeResErr.empty())  fTimeResErr.clear();
}
//------------------------------------------------------------------
/// Prints details of SFTimingRes class object.
void SFTimingRes::Print(void){
  cout << "\n--------------------------------------------" << endl;
  cout << "This is print out of SFTimingRes class object" << endl;
  cout << "Experimental series number " << fSeriesNo << endl;
  cout << "Threshold type: " << fThreshold << endl;
  cout << "Analysis method: " << fMethod << endl;
  cout << "--------------------------------------------\n" << endl;
}
//------------------------------------------------------------------
