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
/// Standard constructor
/// \param seriesNo is number of experimental series to be analyzed.
SFTimingRes::SFTimingRes(int seriesNo): fSeriesNo(seriesNo),
                                        fData(nullptr), 
                                        fTResGraph(nullptr),
                                        fTResECutGraph(nullptr),
                                        fTimeRes(-1),
                                        fTimeResErr(-1),
                                        fTimeResECut(-1),
                                        fTimeResECutErr(-1) {
  
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFTimingRes constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
SFTimingRes::~SFTimingRes(){
 
  if(fData!=nullptr) 
    delete fData; 
}
//------------------------------------------------------------------
///Private method to get ratio histograms necesarry to impose cuts. 
///Depending on measurement type single or double gaussian function
///is fitted to the histograms.
bool SFTimingRes::LoadRatios(void){
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  
  TString cut = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0 && ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0"; 
  fRatios = fData->GetCustomHistograms(SFSelectionType::LogSqrtPERatio, cut);
  
  std::vector <TF1*> fun;
  double min, max;
  
  for(int i=0; i<npoints; i++){
    if(collimator.Contains("Lead")){
      fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -1, 1));
      fun[i]->SetParameter(0, fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));   //thin gauss
      fun[i]->SetParameter(1, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin()));
      fun[i]->SetParameter(2, 6E-2);
      fun[i]->SetParameter(3, 0.5*fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));   //thick gauss
      if(i<npoints/2)
        fun[i]->SetParameter(4, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())+0.2);
      else
	fun[i]->SetParameter(4, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())-0.2);
      fun[i]->SetParameter(5, 2E-1);
      fRatios[i]->Fit(fun[i], "QR");
    }
    else if(collimator.Contains("Electronic")){
      min = fRatios[i]->GetMean()-2*fRatios[i]->GetRMS();
      max = fRatios[i]->GetMean()+2*fRatios[i]->GetRMS();
      fun.push_back(new TF1("fun", "gaus", min, max));   //single gauss
      fRatios[i]->Fit(fun[i], "QR");
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
///If measurement was taken with electronic collimator - single Gaussian.
bool SFTimingRes::AnalyzeNoECut(void){
  
  std::cout << "----- Inside SFTimingRes::AnalyzeNoECut()" << std::endl;
  std::cout << "----- Series number: " << fSeriesNo << std::endl;  
    
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measIDs = fData->GetMeasurementsIDs();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  TString cut;
  double mean, sigma;
  double tResAv = 0;
  double tResAvErr = 0;
  
  if(fRatios.empty()) LoadRatios();
  std::vector <TF1*> fun;

  TString gname = Form("timeDiff_S%i", fSeriesNo);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("T_{D} [ns]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    
    mean = fRatios[i]->GetFunction("fun")->GetParameter(1);
    sigma = fRatios[i]->GetFunction("fun")->GetParameter(2);
    
    if(collimator.Contains("Lead")){
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-0.5*sigma, mean+0.5*sigma);
      fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
      fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -100, 100));
      fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
      fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin())/4);
      if(i<npoints/2)
        fun[i]->SetParameter(4, fT0Diff[i]->GetMean()-5);
      else
        fun[i]->SetParameter(4,fT0Diff[i]->GetMean()+5);
      fun[i]->SetParameter(5,fT0Diff[i]->GetRMS()*2);
      fT0Diff[i]->Fit(fun[i],"QR");
    }
    else if(collimator.Contains("Electronic")){
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-3*sigma,  mean+3*sigma);
      fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
      fT0Diff.back()->Rebin(2);
      fun.push_back(new TF1("fun", "gaus", -50, 50));
      fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
      fT0Diff[i]->Fit(fun[i], "QR");
    }

    int parNum=0;
    
    if(collimator=="Lead"){
      if(fun[i]->GetParameter(2)<fun[i]->GetParameter(5))
        parNum=2;
      else 
        parNum=5;
    }
    else 
      parNum=2;
    
    fTimeResAll.push_back(fun[i]->GetParameter(parNum));
    fTimeResAllErr.push_back(fun[i]->GetParError(parNum));
    
    graph->SetPoint(i, positions[i], fun[i]->GetParameter(parNum-1));
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun[i]->GetParameter(parNum));
    
    tResAv += fTimeResAll[i]*(1./pow(fTimeResAllErr[i], 2));
    tResAvErr += 1./pow(fTimeResAllErr[i], 2);
  }

  
  fTResGraph  = graph;
  fTimeRes    = tResAv/tResAvErr;
  fTimeResErr = sqrt(1./tResAvErr);
  
  std::cout << "Average timing resolution: " << fTimeRes << " +/- " << fTimeResErr << " ns" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
///Method for timing resolution determination with an energy cut on 511 keV peak.
///Timing resolution based on the Gaussian fit to the histogram - sigma of the fitted
///function. Regardless the measurement type - always sigle Gauss fitted.
bool SFTimingRes::AnalyzeWithECut(void){
  
  std::cout << "----- Inside SFTimingRes::AnalyzeWithECut()" << std::endl;  
  std::cout << "----- Series number " << fSeriesNo << std::endl;
    
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measIDs = fData->GetMeasurementsIDs();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  if(fRatios.empty()) LoadRatios();
  
  fSpecCh0 = fData->GetSpectra(0, SFSelectionType::PE, "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  fSpecCh1 = fData->GetSpectra(1, SFSelectionType::PE, "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0");
  std::vector <SFPeakFinder*> peakFin_ch0;
  std::vector <SFPeakFinder*> peakFin_ch1;
  TF1* fun = new TF1("fun", "gaus", -200, 200);
  double xmin_ch0, xmax_ch0;
  double xmin_ch1, xmax_ch1;
  double mean_ratio, sigma_ratio;
  double mean, sigma;
  //double f = 2*sqrt(2*log(2));   //to recalculate sigma into FWHM
  TString cut;
  
  double center_ch0, delta_ch0;  //changed here for smaller cut
  double center_ch1, delta_ch1;  //
  
  double tResAv = 0;
  double tResAvErr = 0;
  
  TString gname = Form("timeDiffECut_S%i", fSeriesNo);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("mean of time difference distribution [ns]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    peakFin_ch0.push_back(new SFPeakFinder(fSpecCh0[i],false));
    peakFin_ch1.push_back(new SFPeakFinder(fSpecCh1[i],false));
    peakFin_ch0[i]->FindPeakRange(xmin_ch0, xmax_ch0);
    peakFin_ch1[i]->FindPeakRange(xmin_ch1, xmax_ch1);
    
    center_ch0 = xmin_ch0+(xmax_ch0-xmin_ch0)/2.;   //changed here for smaller cut
    delta_ch0  = (xmax_ch0-xmin_ch0)/6.;            //
    center_ch1 = xmin_ch1+(xmax_ch1-xmin_ch1)/2.;   //
    delta_ch1  = (xmax_ch1-xmin_ch1)/6.;            //
    mean_ratio = fRatios[i]->GetFunction("fun")->GetParameter(1);
    sigma_ratio = fRatios[i]->GetFunction("fun")->GetParameter(2);
    
    if(collimator.Contains("Lead"))
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", center_ch0-delta_ch0, center_ch0+delta_ch0, center_ch1-delta_ch1, center_ch1+delta_ch1, mean_ratio-0.5*sigma_ratio, mean_ratio+0.5*sigma_ratio);   //changed here for smaller cut
    else if(collimator.Contains("Electronic"))
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0, center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-3*sigma_ratio, mean_ratio+3*sigma_ratio);   //changed here for smaller cut
      
    fT0DiffECut.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
    mean = fT0DiffECut[i]->GetMean();
    sigma = fT0DiffECut[i]->GetRMS();
    fT0DiffECut[i]->Fit(fun, "Q", "", mean-5*sigma, mean+5*sigma);
    fTimeResECutAll.push_back(fun->GetParameter(2));	//Timing resolution as sigma, if FWHM needed multiply by f
    fTimeResECutAllErr.push_back(fun->GetParError(2));		//FWHM - multiply by f
    
    graph->SetPoint(i, positions[i], fun->GetParameter(1));
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun->GetParameter(2));   //FWHM - multiply by f
    
    tResAv += fTimeResECutAll[i]*(1./pow(fTimeResECutAllErr[i], 2));
    tResAvErr += 1./pow(fTimeResECutAllErr[i], 2);
  }
  
  fTResECutGraph  = graph;
  fTimeResECut    = tResAv/tResAvErr;
  fTimeResECutErr = sqrt(1./tResAvErr);
  
  std::cout << "Average timing resolution: " << fTimeResECut << " +/- " << fTimeResECutErr << " ns" << std::endl; 
  
  return true;
}
//------------------------------------------------------------------
/// Returns vector of histograms, which contains timing resolution 
/// spectra for whole series.
/// \param type - type of analysis (1 - with cut, 0 - no cut)
std::vector <TH1D*> SFTimingRes::GetT0Diff(bool type){
  
  if((type==0 && fT0Diff.empty()) || 
     (type==1 && fT0DiffECut.empty())){
    std::cerr << "##### Error in SFTimingRes::GetT0Diff()!" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
  }
  
  if(type==0)
    return fT0Diff;
  else if(type==1)
    return fT0DiffECut;
}
//------------------------------------------------------------------
/// Returns graph ch_0.T0-ch_1.T0 vs. source position.
/// \param type - type of analysis (0 - no cut, 1 - with cut)
TGraphErrors* SFTimingRes::GetTimingResGraph(bool type){
  
  if((type==0 && fTResGraph==nullptr) ||
     (type==1 && fTResECutGraph==nullptr)){
    std::cerr << "##### Error in SFTimingRes::GetTimingResGraph()!" << std::endl;
  std::cerr << "No graph available!" << std::endl;
  std::abort();
  }
  
  if(type==0)
    return fTResGraph;
  else if(type==1)
    return fTResECutGraph;
}
//------------------------------------------------------------------
/// Returns vector containing average timing resolution for the 
/// series along with the uncertainty.
/// \param type - type of analysis (0 - no cut, 1 - with cut)
std::vector <double> SFTimingRes::GetTimingResolution(bool type){
    
  std::vector <double> tmp;
  
  if(type==0){
    tmp.push_back(fTimeRes);
    tmp.push_back(fTimeResErr);
  }
  else if(type==1){
    tmp.push_back(fTimeResECut);
    tmp.push_back(fTimeResECutErr);
  }

  if(tmp[0]==-1 || tmp[1]==-1){
    std::cerr << "##### Error in SFTimingRes::GetTimingResolution()!" << std::endl;
    std::cerr << "Incorrect values: " << tmp[0] << " +/- " << tmp[1] << " ns" << std::endl;
    std::abort();
  }

  return tmp;
}
//------------------------------------------------------------------
/// Returns vector containing timing resolution values for 
/// all measurements in the series.
/// \param type - type of analysis (0 - no cut, 1 - with cut)
std::vector <double> SFTimingRes::GetTimingResolutionAll(bool type){

  if((type==0 && fTimeResAll.empty()) || 
     (type==1 && fTimeResECutAll.empty())){
    std::cerr << "##### Error in SFTimingRes::GetTimingResolutionAll()!" << std::endl;
    std::cerr << "No results available!" << std::endl;
    std::abort();
  }
  
  if(type==0)
    return fTimeResAll;
  else if(type==1)
    return fTimeResECutAll;
}
//------------------------------------------------------------------
/// Returns vector containing uncertainties of timing resolutions for 
/// all measurements in the series.
/// \param type - type of analysis (0 - no cut, 1 - with cut)
std::vector <double> SFTimingRes::GetTimingResolutionAllErr(bool type){
  
  if((type==0 && fTimeResAllErr.empty()) || 
     (type==1 && fTimeResECutAllErr.empty())){
    std::cerr << "##### Error in SFTimingRes::GetTimingResolutionAllErr()!" << std::endl;
    std::cerr << "No results available!" << std::endl;
    std::abort();
  }
  
  if(type==0)
    return fTimeResAllErr;
  else if(type==1)
    return fTimeResECutAllErr;   
}
//------------------------------------------------------------------
/// Returns vector of ratio histograms used for cut on scattered events.
std::vector <TH1D*> SFTimingRes::GetRatios(void){
    
  if(fRatios.empty()){
    std::cout << "##### Error in SFTimingRes::GetRatios()! Empty vector!" << std::endl;
    std::abort();
  }
  return fRatios; 
}
//------------------------------------------------------------------
/// Returns vector of PE spectra used for 511 keV energy cut. 
/// \param ch - channel number (0 or 1).
std::vector <TH1D*> SFTimingRes::GetSpectra(int ch){
 
  if((ch==0 && fSpecCh0.empty()) || 
     (ch==1 && fSpecCh1.empty())){
    std::cerr << "##### Error in SFTimingRes::GetSpectra()!" << std::endl; 
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }

  if(ch==0)       
    return fSpecCh0; 
  else if(ch==1) 
    return fSpecCh1;
  else{
    std::cerr << "##### Error in SFTimingRes::GetSpectra()" << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
/// Prints details of SFTimingRes class object.
void SFTimingRes::Print(void){
  std::cout << "\n--------------------------------------------" << std::endl;
  std::cout << "This is print out of SFTimingRes class object" << std::endl;
  std::cout << "Experimental series number " << fSeriesNo << std::endl;
  std::cout << "--------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
