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
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed. 
SFAttenuation::SFAttenuation(int seriesNo): fSeriesNo(seriesNo),
                                            fData(nullptr),
                                            fAttnGraph(nullptr),
                                            fSigmaGraph(nullptr),
                                            fAttnGraphCh0(nullptr),
                                            fAttnGraphCh1(nullptr) {
  
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFAttenuation constructor!";
  }
  
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cout << "##### Error in SFAttenuation constructor! Non-regular series!" << std::endl;
    throw "##### Exception in SFAttenuation constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
SFAttenuation::~SFAttenuation(){
    
  if(fData!=nullptr) 
     delete fData;
}
//------------------------------------------------------------------
/// Method to determine attenuation length used in Pauwels et al., JINST 8 (2013) P09019.
/// For both ends of the fiber one value is calculated, since combined signal from both channels
/// is taken into account.
bool SFAttenuation::AttAveragedCh(void){
 
  std::cout << "\n----- Inside SFAttenuation::AttAveragedCh() for series " << fSeriesNo << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  TString cut =  "ch_0.fPE>0 && ch_1.fPE>0 && ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590";
  fRatios = fData->GetCustomHistograms(SFSelectionType::LogSqrtPERatio, cut);
  
  double mean, sigma;
  double fit_min, fit_max;
  std::vector <TF1*> fun;
  
  TString gname = Form("att_s%i", fSeriesNo);
  fAttnGraph = new TGraphErrors(npoints);
  fAttnGraph->GetXaxis()->SetTitle("source position [mm]");
  fAttnGraph->GetYaxis()->SetTitle("ln(M_{LR})");
  fAttnGraph->SetTitle(gname);
  fAttnGraph->SetName(gname);
  fAttnGraph->SetMarkerStyle(4);
  
  gname = Form("sigmas_s%i", fSeriesNo);
  fSigmaGraph = new TGraphErrors(npoints);
  fSigmaGraph->GetXaxis()->SetTitle("source position [mm]");
  fSigmaGraph->GetYaxis()->SetTitle("#sigma M_{LR}");
  fSigmaGraph->SetTitle(gname);
  fSigmaGraph->SetName(gname);
  fSigmaGraph->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    mean = fRatios[i]->GetMean();
    sigma = fRatios[i]->GetRMS();
    if(collimator.Contains("Lead")){
      fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -1, 1));
      fun[i]->SetParameter(0, fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));		//thin gauss
      fun[i]->SetParameter(1, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin()));
      fun[i]->SetParameter(2, 6E-2);
      fun[i]->SetParameter(3, 0.5*fRatios[i]->GetBinContent(fRatios[i]->GetMaximumBin()));	//thick gauss
      if(i<npoints/2)
        fun[i]->SetParameter(4, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())+0.2);
      else
      fun[i]->SetParameter(4, fRatios[i]->GetBinCenter(fRatios[i]->GetMaximumBin())-0.2);
      fun[i]->SetParameter(5, 2E-1);
    }
    else if(collimator.Contains("Electronic")){
      fit_min = mean - sigma;
      fit_max = mean + sigma;
      fun.push_back(new TF1("fun", "gaus", fit_min, fit_max));	//single gauss
    }
    fRatios[i]->Fit(fun[i],"QR");
    fAttnGraph->SetPoint(i, positions[i], fun[i]->GetParameter(1));
    fAttnGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun[i]->GetParError(1));
    fSigmaGraph->SetPoint(i, positions[i], fun[i]->GetParameter(2));
    fSigmaGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), fun[i]->GetParError(2));
  }
  
  TF1 *fpol1 = new TF1("fpol1", "pol1", positions[0], positions[npoints-1]);
  fpol1->SetParameters(-0.15, 0.005);
  fAttnGraph->Fit(fpol1, "QR");
  
  fResults.fAttCombPol1    = fabs(1./fpol1->GetParameter(1));
  fResults.fAttCombPol1Err = fpol1->GetParError(1)/pow(fpol1->GetParameter(1), 2);

  std::cout << "Attenuation lenght from pol1 fit is: " << fResults.fAttCombPol1  
            << " +/- " << fResults.fAttCombPol1Err << " mm\n" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
double myPol3(double *x, double *par){
    
  double xx = (x[0]-par[3]);  
  double f = par[0] + par[1]*(xx + par[2]*pow(xx,3));
  return f;  
}
//------------------------------------------------------------------
bool SFAttenuation::Fit3rdOrder(void){
    
  if(fAttnGraph==nullptr)
      AttAveragedCh();
  
  TF1 *fpol3 = new TF1("fpol3", myPol3, -50, 150, 4); 
  double fiberLengthHalf = fData->GetFiberLength()/2;
  fpol3->SetParameter(0, fAttnGraph->GetFunction("fpol1")->GetParameter(0));
  fpol3->SetParameter(1, fAttnGraph->GetFunction("fpol1")->GetParameter(1));
  fpol3->SetParameter(3, fiberLengthHalf);

  fpol3->SetLineColor(kBlue-7);
  
  fAttnGraph->Fit(fpol3, "QR+");
  
  fResults.fAttCombPol3    = fabs(1./fpol3->GetParameter(1));
  fResults.fAttCombPol3Err = fpol3->GetParError(1)/pow(fpol3->GetParameter(1), 2);
  
  std::cout << "Attenuation lenght from pol3 fit is: " << fResults.fAttCombPol3  
            << " +/- " << fResults.fAttCombPol3Err << " mm\n" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
/// Method to determine attenuation length for both channels independently. 
/// If series was measured with lead collimator peak position is determied 
/// with the FindPeakNoBackground() method of the SFPeakFinder class. If
/// series was measured with electronic collimator - FindPeakFit() method
/// of the SFPeakFinder class is used.
/// \param ch - channel number
bool SFAttenuation::AttSeparateCh(int ch){
 
  std::cout << "\n----- Inside SFAttenuation::AttSeparateCh() for series " << fSeriesNo << std::endl;
  std::cout << "----- Analyzing channel " << ch << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  TString cut = Form("ch_%i.fT0>0 && ch_%i.fT0<590 && ch_%i.fPE>0", ch, ch, ch);
  std::vector <TH1D*> spectra = fData->GetSpectra(ch, SFSelectionType::PE, cut);
  
  TString gname = Form("att_s%i_ch%i", fSeriesNo, ch);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  std::vector <SFPeakFinder*> peakfin;
  PeakParams peakParams;
  
  for(int i=0; i<npoints; i++){
    peakfin.push_back(new SFPeakFinder(spectra[i], false));
    peakfin[i]->FindPeakFit();
    peakParams = peakfin[i]->GetParameters();
    graph->SetPoint(i, positions[i], peakParams.fPosition);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), peakParams.fPositionErr);
  }
  
  //----- fitting 
  TF1 *fexp = new TF1("fexp", "expo", positions[0], positions[npoints-1]);
  graph->Fit(fexp, "QR");
  
  //----- calculating attenuation length
  
  double attenuation = fabs(1./fexp->GetParameter(1));
  double att_error = fexp->GetParError(1)/pow(fexp->GetParameter(1),2);
  
  std::cout << "\n\tAttenuation for channel "<< ch << ": " << attenuation << " +/- " << att_error << " mm\n" << std::endl;
  
  if(ch==0){
    fResults.fAttCh0    = attenuation;
    fResults.fAttCh0Err = att_error;
    fSpectraCh0   = spectra;
    fAttnGraphCh0 = graph;
  }
  else if(ch==1){
    fResults.fAttCh1    = attenuation;
    fResults.fAttCh1Err = att_error;
    fSpectraCh1   = spectra;
    fAttnGraphCh1 = graph;
  }
  
  return true;
}
//------------------------------------------------------------------
///Returns attenuation graph created in averaged channels method i.e. AttAveragedCh().
TGraphErrors* SFAttenuation::GetAttGraph(void){
    
  if(fAttnGraph==nullptr){
    std::cerr << "##### Error in SFAttenuation::GetAttGraph(). Empty pointer!" << std::endl;
    std::abort();
  }
  return fAttnGraph;
}
//------------------------------------------------------------------
///Returns attenuation graph for requested channel ch. Graph produced 
///with separate channels method - AttSeparateCh().
TGraphErrors* SFAttenuation::GetAttGraph(int ch){
    
   if((ch==0 && fAttnGraphCh0==nullptr) || (ch==1 && fAttnGraphCh1==nullptr)){
     std::cerr << "##### Error in SFAttenuation::GetAttnGraph(int). Empty pointer!" << std::endl;
     std::abort();
   }
   if(ch==0) return fAttnGraphCh0;
   else if(ch==1) return fAttnGraphCh1;
}
//------------------------------------------------------------------
TGraphErrors* SFAttenuation::GetSigmaGraph(void){
    
  if(fSigmaGraph==nullptr){
    std::cerr << "##### Error in SFAttenuation::GetSigmaGraph(). Empty pointer!" << std::endl;
    std::abort();
  }
  return fSigmaGraph;
}
//------------------------------------------------------------------
///Returns vector containing PE spectra used in determination of attenuation
///length with separate channels method i.e. AttSeparateCh().
///\param ch - channel number
std::vector <TH1D*> SFAttenuation::GetSpectra(int ch){
    
  if((ch==0 && fSpectraCh0.empty()) || (ch==1 && fSpectraCh1.empty())){
    std::cerr << "##### Error in SFAttenuation::GetSpectra(). Empty vector!" << std::endl;
    std::abort();
  }
  if(ch==0) return fSpectraCh0;
  else if(ch==1) return fSpectraCh1;
}
//------------------------------------------------------------------
///Returns vector containing peaks (spectra after background subtraction with 
///SFPeakFinder class) used in determination of attenuation length with separate
///channels method i.e. AttSeparateCh(). 
///\param ch - channel number.
std::vector <TH1D*> SFAttenuation::GetPeaks(int ch){
    
  if((ch==0 && fPeaksCh0.empty()) || (ch==1 && fPeaksCh1.empty())){
    std::cerr << "##### Error in SFAttenuation::GetPeaks(). Empty vector!" << std::endl;
    std::abort();
  }
  if(ch==0) return fPeaksCh0;
  else if(ch==1) return fPeaksCh1;
}
//------------------------------------------------------------------
///Returns vector containing histograms with signal ratios from both channels. 
///Histograms are used during averaged channels analysis i.e. in AttAveragedCh().
std::vector <TH1D*> SFAttenuation::GetRatios(void){
    
  if(fRatios.empty()){
    std::cerr << "#### Error in SFAttenuation::GetRatios(). Empty vector!" << std::endl;
    std::abort();
   }
   return fRatios;
}
//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFAttenuation::Print(void){
 std::cout << "\n-------------------------------------------" << std::endl;
 std::cout << "This is print out of SFAttenuation class object" << std::endl;
 std::cout << "Experimental series number " << fSeriesNo << std::endl;
 std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------