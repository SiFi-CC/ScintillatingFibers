// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFEnergyRes.cc              *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyRes.hh"
#include <cmath>

ClassImp(SFEnergyRes);

//------------------------------------------------------------------
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed. 
SFEnergyRes::SFEnergyRes(int seriesNo): fSeriesNo(seriesNo),
                                        fData(nullptr),
                                        fEnergyResGraphCh0(nullptr),
                                        fEnergyResGraphCh1(nullptr),
                                        fEnergyResGraphSum(nullptr) {

  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFEnergyRes constructor!";
  }
  
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cout << "##### Warning in SFEnergyRes constructor!" << std::endl;
    std::cout << "Calculating energy resolution with non-regular series!" << std::endl;
  }
  
  //--- getting necessary PE spectra
  TString cut_ch0 = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0";
  TString cut_ch1 = "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0";
  fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::PE, cut_ch0);
  fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::PE, cut_ch1);
   
}
//------------------------------------------------------------------
/// Default destructor.
SFEnergyRes::~SFEnergyRes(){
  
  if(fData!=nullptr)
    delete fData;
}
//------------------------------------------------------------------
bool SFEnergyRes::CalculateEnergyRes(int ch){
    
  std::cout << "\n----- Inside SFEnergyRes::CalculateEnergyRes()" << std::endl;
  std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
  std::cout << "----- Analyzing channel: " << ch << std::endl; 
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <SFPeakFinder*> peakFin;
  
  for(int i=0; i<npoints; i++){
    if(ch==0)
      peakFin.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
    else if(ch==1)
      peakFin.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
    else{
      std::cerr << "##### Error in SFEnergyRes::CalculateEnergyRes() for ch"
                <<  ch << std::endl;
      std::cerr << "Incorrect channel number!" << std::endl;
      return false;
    }
  }
  
  TString gname = Form("ER_s%i_ch%i", fSeriesNo, ch);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  PeakParams parameters;
  double enRes = 0;
  double enResErr = 0;
  double enResAve = 0;
  double enResAveErr = 0;
  
  for(int i=0; i<npoints; i++){
    peakFin[i]->FindPeakFit();
    parameters = peakFin[i]->GetParameters();
    enRes = parameters.fSigma/parameters.fPosition;
    enResErr = enRes * sqrt(pow(parameters.fPositionErr, 2)/pow(parameters.fPosition, 2) +
                            pow(parameters.fSigmaErr, 2)/pow(parameters.fSigma, 2));
    enRes = enRes*100;
    enResErr = enResErr*100;
    graph->SetPoint(i, positions[i], enRes);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), enResErr);
    enResAve += enRes * (1./pow(enResErr,2));
    enResAveErr += (1./pow(enResErr,2));
  }
  
  enResAve = enResAve/enResAveErr;
  enResAveErr = sqrt(1./enResAveErr);
  
  std::cout << "Average energy resolution for channel " << ch 
            << ": " << enResAve << " +/- " << enResAveErr << " % \n" << std::endl; 
 
  if(ch==0){
      fEnergyResGraphCh0 = graph;
      fResults.fEnergyResCh0 = enResAve;
      fResults.fEnergyResCh0Err = enResAveErr;
  }
  else if(ch==1){
      fEnergyResGraphCh1 = graph;
      fResults.fEnergyResCh1 = enResAve;
      fResults.fEnergyResCh1Err = enResAveErr;
  }
            
  return true;
}
//------------------------------------------------------------------
bool SFEnergyRes::CalculateEnergyRes(void){
    
  std::cout << "\n----- Inside SFEnergyRes::CalculateEnergyRes()" << std::endl;
  std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measIDs = fData->GetMeasurementsIDs();
  std::vector <SFPeakFinder*> peakFin;
  
  SFAttenuation *att;
  
  try{
    att = new SFAttenuation(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in SFEnergyRes::CalculateEergyRes()" << std::endl;
    std::cerr << "Couldn't access attenuation length data!" << std::endl;
    std::abort();
  }
  
  att->AttAveragedCh();
  AttenuationResults results = att->GetResults();
  
  TString cut_ch0 = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0";
  TString cut_ch1 = "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0";
  TString cut = cut_ch0 + std::string(" && ") + cut_ch1;
  
  std::vector <double> customNumCh0(2);
  customNumCh0[1] = results.fAttCombPol1; 
  
  std::vector <double> customNumCh1(2);
  customNumCh1[1] = results.fAttCombPol1; 
  
  std::vector <double> customNumSum(4);
  customNumSum[1] = results.fAttCombPol1; 
  customNumSum[3] = results.fAttCombPol1;
  
  TString gname = Form("ER_s%i_ave", fSeriesNo);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  PeakParams parameters;
  double distCh0, distCh1;
  double enRes, enResErr;
  double enResAve, enResAveErr;
  
  for(int i=0; i<npoints; i++){
    distCh0 = positions[i];
    distCh1 = 100. - positions[i]; 
    customNumCh0[0] = -distCh0;
    customNumCh1[0] = -distCh1;
    customNumSum[0] = -distCh0;
    customNumSum[2] = -distCh1;
    
    fSpectraCorrCh0.push_back(fData->GetCustomHistogram(0, SFSelectionType::PEAttCorrected, 
                                                        cut_ch0, measIDs[i], customNumCh0));
    fSpectraCorrCh1.push_back(fData->GetCustomHistogram(1, SFSelectionType::PEAttCorrected, 
                                                        cut_ch1, measIDs[i], customNumCh1));
    fSpectraSum.push_back(fData->GetCustomHistogram(SFSelectionType::PEAttCorrectedSum, 
                                                    cut, measIDs[i], customNumSum));
    
    peakFin.push_back(new SFPeakFinder(fSpectraSum[i], 0));
    peakFin[i]->FindPeakFit();
    parameters = peakFin[i]->GetParameters();
    
    enRes = parameters.fSigma/parameters.fPosition;
    enResErr = enRes * sqrt(pow(parameters.fPositionErr, 2)/pow(parameters.fPosition, 2) + 
                            pow(parameters.fSigmaErr, 2)/pow(parameters.fSigma, 2));
    enRes = enRes*100;
    enResErr = enResErr*100;
    graph->SetPoint(i, positions[i], enRes);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), enResErr);
    enResAve += enRes * (1./pow(enResErr, 2));
    enResAveErr += (1./pow(enResErr, 2));
  }
  
  fResults.fEnergyResSum = enResAve/enResAveErr;
  fResults.fEnergyResSumErr = sqrt(1./enResAveErr);
  fEnergyResGraphSum = graph;
  
  if(std::isnan(fResults.fEnergyResSum) || fResults.fEnergyResSum<0 || fResults.fEnergyResSum>100){
    std::cout << "##### Warning in SFEnergyRes class!" << std::endl;
    std::cout << "Incorrect energy resolution value: " << fResults.fEnergyResSum  
              << " +/- " << fResults.fEnergyResSumErr << std::endl;
    fResults.fEnergyResSum = 0;
    fResults.fEnergyResSumErr = 0;
  }
  
  std::cout << "Average energy resolution calculated from summed and attenuation length corrected histograms is: "  << fResults.fEnergyResSum << " +/- " << fResults.fEnergyResSumErr << " % \n" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyRes::GetEnergyResolutionGraph(int ch){
    
  if((ch==0 && fEnergyResGraphCh0==nullptr) || 
     (ch==1 && fEnergyResGraphCh1==nullptr)) {
      std::cerr << "##### Error in SFEnergyRes::GetEnergyResolutionGraph() for ch" 
                << ch << std::endl;
      std::cerr << "Requested graph doesnt exist!" << std::endl;
      std::abort();
  }
  
  if(ch==0)
    return fEnergyResGraphCh0;
  else if(ch==1) 
    return fEnergyResGraphCh1;
  else{
    std::cerr << "##### Error in SFEnergyRes::GetEnergyResolutionGraph() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyRes::GetEnergyResolutionGraph(void){
 
  if(fEnergyResGraphSum==nullptr){
    std::cerr << "##### Error in SFEnergyRes::GetEnergyResolution()" << std::endl;  
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();
  }
  
  return fEnergyResGraphSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyRes::GetSpectra(int ch){
    
  if((ch==0 && fSpectraCh0.empty()) ||
     (ch==1 && fSpectraCh1.empty())) {
      std::cerr << "##### Error in SFEnergyRes::GetSpectra() fo ch" << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fSpectraCh0;
  else if(ch==1)
    return fSpectraCh1;
  else{
    std::cerr << "##### Error in SFEnergyRes::GetSpectra for ch" << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyRes::GetSpectraCorrected(int ch){
    
  if((ch==0 && fSpectraCorrCh0.empty()) ||
     (ch==1 && fSpectraCorrCh1.empty())) {
      std::cerr << "##### Error in SFEnergyRes::GetSpectraCorrected() for ch" 
                << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fSpectraCorrCh0;
  else if(ch==1)
    return fSpectraCorrCh1;
  else{
    std::cerr << "##### Error in SFEnergyRes::GetSpectraCorrected() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyRes::GetSpectraSum(void){
  
  if(fSpectraSum.empty()){
    std::cerr << "##### Error in SFEnergyRes::GetSpectraSum()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fSpectraSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyRes::GetPeaks(void){
  
  if(fPeaksSum.empty()){
    std::cerr << "##### Error in SFEnergyRes::GetPeaks()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fPeaksSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyRes::GetPeaks(int ch){
    
  if((ch==0 && fPeaksCh0.empty()) ||
     (ch==1 && fPeaksCh1.empty())) {
      std::cerr << "##### Error in SFEnergyRes::GetPeaks() for ch" 
                << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fPeaksCh0;
  else if(ch==1)
    return fPeaksCh1;
  else{
    std::cerr << "##### Error in SFEnergyRes::GetPeaks() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
/// Prints details of the SFEnergyResolution class object.
void SFEnergyRes::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl;
  std::cout << "This is print out of SFEnergyRes class object" << std::endl;
  std::cout << "Experimental series number " << fSeriesNo << std::endl;
  std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------