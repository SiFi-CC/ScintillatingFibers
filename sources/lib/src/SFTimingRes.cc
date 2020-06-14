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
                                        fTResECutGraph(nullptr) {
  
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
  TString sipm = fData->GetSiPM();
  
  TString cut = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0 && ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0"; 
  fRatios = fData->GetCustomHistograms(SFSelectionType::LogSqrtPERatio, cut);
  
  std::vector <TF1*> fun;
  
  for(int i=0; i<npoints; i++){
    if(collimator.Contains("Lead")){
      SFTools::RatiosFitDoubleGauss(fRatios, 5);
    }
    else if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
      SFTools::RatiosFitDoubleGauss(fRatios, 2);
    }
    else if(collimator.Contains("Electronic") && sipm.Contains("Hamamatsu")){
      SFTools::RatiosFitGauss(fRatios, 1);  
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
  TString sipm = fData->GetSiPM();
  TString testBench = fData->GetTestBench();
  TString fiber = fData->GetFiber();
  
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
    
    int parNum = 0;  
    
    if(collimator.Contains("Lead")){
      
    if(fRatios[i]->GetFunction("fDGauss")->GetParameter(0) > 
       fRatios[i]->GetFunction("fDGauss")->GetParameter(3))
        parNum = 1;
    else 
        parNum = 4;
    
      mean = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum);
      sigma = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum+1);
    
      //cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-0.5*sigma, mean+0.5*sigma);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>15 && ch_1.fPE>15 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-0.5*sigma, mean+0.5*sigma);
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
    else if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
      
      if(fRatios[i]->GetFunction("fDGauss")->GetParameter(0) > 
         fRatios[i]->GetFunction("fDGauss")->GetParameter(3))
          parNum = 1;
      else 
          parNum = 4;
    
      mean = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum);
      sigma = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum+1);
      
      //cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-2*sigma,  mean+2*sigma);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>20 && ch_1.fPE>20 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-2*sigma,  mean+2*sigma);
      fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
      fT0Diff.back()->Rebin(2);
      fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -30, 30));
      
      fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
      fun[i]->SetParLimits(2, 0, 10);
      fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin())/10.);
      fun[i]->SetParameter(4, fT0Diff[i]->GetMean());
      fun[i]->SetParameter(5, fT0Diff[i]->GetRMS()*10);
      fun[i]->SetParLimits(5, 0, 20);
      fT0Diff[i]->Fit(fun[i], "R");
    }
    else if(collimator.Contains("Electronic") && sipm.Contains("Hamamatsu")){
    
      mean = fRatios[i]->GetFunction("fGauss")->GetParameter(1);
      sigma = fRatios[i]->GetFunction("fGauss")->GetParameter(2);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 && ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-3*sigma,  mean+3*sigma);
      fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
      fT0Diff.back()->Rebin(2);
      
      fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -50, 50));
      fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
      fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
      fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
      fun[i]->SetParLimits(2, 0, 50);
      if(fiber.Contains("LuAG")){
        fun[i]->FixParameter(3, 0);
        fun[i]->FixParameter(4, 0);
        fun[i]->FixParameter(5, 0);
      }
      else if(fiber.Contains("LYSO")){
        fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin())/10);
        if(i<npoints/2)
          fun[i]->SetParameter(4, fT0Diff[i]->GetMean()-2);
        else
          fun[i]->SetParameter(4, fT0Diff[i]->GetMean()+2);
        fun[i]->SetParameter(5, fT0Diff[i]->GetRMS()*5);
        fun[i]->SetParLimits(5, 0, 50);
      }
      else if(fiber.Contains("GAGG")){
        fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin())/10);
        if(i<npoints/2)
          fun[i]->SetParameter(4, fT0Diff[i]->GetMean()-1);
        else
          fun[i]->SetParameter(4, fT0Diff[i]->GetMean()+1);
        fun[i]->SetParameter(5, fT0Diff[i]->GetRMS()*2);
        fun[i]->SetParLimits(5, 0, 50);
      }
      fT0Diff[i]->Fit(fun[i], "QR");
    }

    parNum=0;
    
    if(fun[i]->GetNpar()>3){
      if(fun[i]->GetParameter(0)>fun[i]->GetParameter(3))
        parNum=2;
      else 
        parNum=5;
    }
    else 
      parNum=2;
    
    fResults.fTimeResAll.push_back(fun[i]->GetParameter(parNum));
    fResults.fTimeResAllErr.push_back(fun[i]->GetParError(parNum));
    
    graph->SetPoint(i, positions[i], fun[i]->GetParameter(parNum-1));
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun[i]->GetParameter(parNum));
    
    tResAv += fResults.fTimeResAll[i]*(1./pow(fResults.fTimeResAllErr[i], 2));
    tResAvErr += 1./pow(fResults.fTimeResAllErr[i], 2);
  }

  
  fTResGraph           = graph;
  fResults.fTimeRes    = tResAv/tResAvErr;
  fResults.fTimeResErr = sqrt(1./tResAvErr);
  
  std::cout << "Average timing resolution: " << fResults.fTimeRes 
            << " +/- " << fResults.fTimeResErr << " ns" << std::endl;
  
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
  TString sipm = fData->GetSiPM();
  
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
  double const_1, const_2;
  //double f = 2*sqrt(2*log(2));   //to recalculate sigma into FWHM
  TString cut;
  
  double center_ch0, delta_ch0;  //change here for smaller cut
  double center_ch1, delta_ch1;  //
  
  double tResAv = 0;
  double tResAvErr = 0;
  
  int parNo = 0;
  
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
    
    center_ch0 = xmin_ch0+(xmax_ch0-xmin_ch0)/2.;   //change here for smaller cut
    delta_ch0  = (xmax_ch0-xmin_ch0)/6.;            //
    center_ch1 = xmin_ch1+(xmax_ch1-xmin_ch1)/2.;   //
    delta_ch1  = (xmax_ch1-xmin_ch1)/6.;            //
    
    if(collimator.Contains("Lead")){
      const_1 = fRatios[i]->GetFunction("fDGauss")->GetParameter(0); 
      const_2 = fRatios[i]->GetFunction("fDGauss")->GetParameter(3);
      if(const_1 > const_2)
          parNo = 1;
      else
          parNo = 4;
      mean_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo);
      sigma_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo+1);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f", center_ch0-delta_ch0, center_ch0+delta_ch0, center_ch1-delta_ch1, center_ch1+delta_ch1, mean_ratio-0.5*sigma_ratio, mean_ratio+0.5*sigma_ratio);   //changed here for smaller cut
    }
    else if(collimator.Contains("Electronic") && sipm.Contains("Hamamatsu")){
      parNo = 1;
      mean_ratio = fRatios[i]->GetFunction("fGauss")->GetParameter(parNo);
      sigma_ratio = fRatios[i]->GetFunction("fGauss")->GetParameter(parNo+1);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0, center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-3*sigma_ratio, mean_ratio+3*sigma_ratio);   //changed here for smaller cut
    }
    else if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
      const_1 = fRatios[i]->GetFunction("fDGauss")->GetParameter(0);
      const_2 = fRatios[i]->GetFunction("fDGauss")->GetParameter(3);
      if(const_1 > const_2)
          parNo = 1;
      else
          parNo = 4;
      mean_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo);
      sigma_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo+1);
      cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0, center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-2*sigma_ratio, mean_ratio+2*sigma_ratio);
    }
    
    fT0DiffECut.push_back(fData->GetCustomHistogram(SFSelectionType::T0Difference, cut, measIDs[i]));
    mean = fT0DiffECut[i]->GetMean();
    sigma = fT0DiffECut[i]->GetRMS();
    fT0DiffECut[i]->Fit(fun, "Q", "", mean-5*sigma, mean+5*sigma);
    fResults.fTimeResECutAll.push_back(fun->GetParameter(2));   //Timing resolution as sigma, if FWHM needed multiply by f
    fResults.fTimeResECutAllErr.push_back(fun->GetParError(2)); //FWHM - multiply by f
    
    graph->SetPoint(i, positions[i], fun->GetParameter(1));
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun->GetParameter(2));   //FWHM - multiply by f
    
    tResAv += fResults.fTimeResECutAll[i]*(1./pow(fResults.fTimeResECutAllErr[i], 2));
    tResAvErr += 1./pow(fResults.fTimeResECutAllErr[i], 2);
  }
  
  fTResECutGraph           = graph;
  fResults.fTimeResECut    = tResAv/tResAvErr;
  fResults.fTimeResECutErr = sqrt(1./tResAvErr);
  
  std::cout << "Average timing resolution: " << fResults.fTimeResECut 
            << " +/- " << fResults.fTimeResECutErr << " ns" << std::endl; 
  
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
