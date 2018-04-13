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
  fSeriesNo  = -1;
  fThreshold = "dummy";
  fMethod    = "dummy";
  fData      = NULL;
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
  if(!(threshold=="ft" || threshold=="cf")){
    throw "##### Error in SFTimingRes constructor! Incorrect threshold type!\nPossible options are: ft, cf";
  }
  if(!(method=="with cut" || method=="no cut")){
    throw "##### Error in SFTimingRes constructor! Incorrect method!\nPossible options are: with cut, no cut";
  }
  fSeriesNo = seriesNo;
  fThreshold = threshold;
  fMethod = method;
  fData = new SFData(fSeriesNo,threshold);
  if(fMethod=="no cut") AnalyzeNoECut();
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
///Method for simple timing resolution determination. Timing resolution spectrum
///is build as difference between channel 0 and 1 signals with additional cut
///on scattered events. Timing resolution and its uncertainty is determined as 
///histogram's mean and RMS. 
bool SFTimingRes::AnalyzeNoECut(void){
  
  int npoints = fData->GetNpoints();
  double *positions = fData->GetPositions();
  
  TString selection_ratio = "log(ch_0.fPE/ch_1.fPE)";
  TString cut_ratio = "ch_0.fT0<590 && ch_0.fPE>0 && ch_0.fT0>0 && ch_1.fT0<590 && ch_1.fT0>0 && ch_1.fPE>0"; 
  TString selection_T0 = "ch_0.fT0-ch_1.fT0";
  TString cut_T0;
  
  vector <TH1D*> ratios = fData->GetRatios(selection_ratio,cut_ratio);
  double mean, sigma;
  
  for(int i=0; i<npoints; i++){
    mean = ratios[i]->GetMean();
    sigma = ratios[i]->GetRMS();
    cut_T0 = Form("ch_0.fT0>0 && ch_1.fT0>0 && log(ch_0.fPE/ch_1.fPE)>%f && log(ch_0.fPE/ch_1.fPE)<%f",mean-sigma, mean+sigma);
    fT0Diff.push_back(fData->GetRatio(selection_T0,cut_T0,positions[i]));
    fTimeRes.push_back(fT0Diff[i]->GetRMS());
    fTimeResErr.push_back(0);
  }
  
  return true;
}
//------------------------------------------------------------------
///Method for timing resolution determination with an energy cut on 511 keV peak.
bool SFTimingRes::AnalyzeWithECut(void){
  
  int npoints = fData->GetNpoints();
  double *positions = fData->GetPositions();
  
  vector <TH1D*> hPEch0 = fData->GetSpectra(0,"fPE","ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  vector <TH1D*> hPEch1 = fData->GetSpectra(1,"fPE","ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0");
  vector <SFPeakFinder*> peakFin_ch0;
  vector <SFPeakFinder*> peakFin_ch1;
  vector <TF1*> fun;
  double xmin_ch0, xmax_ch0;
  double xmin_ch1, xmax_ch1;
  double mean, sigma;
  TString selection = "ch_0.fT0-ch_1.fT0";
  TString cut;
  
  for(int i=0; i<npoints; i++){
    peakFin_ch0.push_back(new SFPeakFinder(hPEch0[i],"511"));
    peakFin_ch1.push_back(new SFPeakFinder(hPEch1[i],"511"));
    peakFin_ch0[i]->FindPeakRange(xmin_ch0,xmax_ch0);
    peakFin_ch1[i]->FindPeakRange(xmin_ch1,xmax_ch1);
    cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f", xmin_ch0,xmax_ch0,xmin_ch1,xmax_ch1);
    fT0Diff.push_back(fData->GetRatio(selection,cut,positions[i]));
    mean = fT0Diff[i]->GetMean();
    sigma = fT0Diff[i]->GetRMS();
    fun.push_back(new TF1("fun","gaus",mean-3*sigma,mean+3*sigma));
    fT0Diff[i]->Fit(fun[i],"QR");
    fTimeRes.push_back(fun[i]->GetParameter(2));
    fTimeResErr.push_back(fun[i]->GetParError(2));
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
/// Prints details of SFTimingRes class object.
void SFTimingRes::Print(void){
  cout << "This is print out of SFTimingRes class object" << endl;
  cout << "Experimental series number " << fSeriesNo << endl;
  cout << "Threshold type: " << fThreshold << endl;
  cout << "Analysis method: " << fMethod << endl;
}
//------------------------------------------------------------------