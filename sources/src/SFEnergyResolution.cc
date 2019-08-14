// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFEnergyResolution.cc         *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyResolution.hh"
#include <cmath>

ClassImp(SFEnergyResolution);

//-----------------------------------------------------------------
/// Default constructor.
SFEnergyResolution::SFEnergyResolution(): fSeriesNo(-1),
                                          fData(nullptr),
                                          fEnergyResCh0(-1),
                                          fEnergyResCh0Err(-1),
                                          fEnergyResCh1(-1),
                                          fEnergyResCh1Err(-1),
                                          fEnergyResSum(-1),
                                          fEnergyResSumErr(-1),
                                          fEnergyResGraphCh0(nullptr),
                                          fEnergyResGraphCh1(nullptr),
                                          fEnergyResGraphSum(nullptr) {
                                              
  std::cout << "#### Warning in SFEnergyResolution constructor!" << std::endl;
  std::cout << "You are using default constructor!" << std::endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed. 
SFEnergyResolution::SFEnergyResolution(int seriesNo): fSeriesNo(seriesNo),
                                                      fData(nullptr),
                                                      fEnergyResCh0(-1),
                                                      fEnergyResCh0Err(-1),
                                                      fEnergyResCh1(-1),
                                                      fEnergyResCh1Err(-1),
                                                      fEnergyResSum(-1),
                                                      fEnergyResSumErr(-1),
                                                      fEnergyResGraphCh0(nullptr),
                                                      fEnergyResGraphCh1(nullptr),
                                                      fEnergyResGraphSum(nullptr) {

  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFEnergyResolution constructor!";
  }
  
  TString desc = fData->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cout << "##### Warning in SFEnergyResolution constructor!" << std::endl;
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
SFEnergyResolution::~SFEnergyResolution(){
  
  if(fData!=nullptr)
    delete fData;
}
//------------------------------------------------------------------
bool SFEnergyResolution::CalculateEnergyRes(int ch){
    
  std::cout << "\n----- Inside SFEnergyResolution::CalculateEnergyRes()" << std::endl;
  std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
  std::cout << "----- Analyzing channel: " << ch << std::endl; 
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <SFPeakFinder*> peakFin;
  std::vector <TH1D*> peaks;
  
  for(int i=0; i<npoints; i++){
    if(ch==0)
      peakFin.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
    else if(ch==1)
      peakFin.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
    else{
      std::cerr << "##### Error in SFEnergyResolution::CalculateEnergyRes() for ch"
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
  
  std::vector <double> parameters;
  double enRes = 0;
  double enResErr = 0;
  double enResAve = 0;
  double enResAveErr = 0;
  
  for(int i=0; i<npoints; i++){
    if(collimator=="Lead"){
      peakFin[i]->FindPeakNoBackground();
      peaks.push_back(peakFin[i]->GetPeak());
    }
    else if(collimator=="Electronic")
      peakFin[i]->FindPeakFit();
    parameters = peakFin[i]->GetParameters();
    enRes = parameters[1]/parameters[0];
    enResErr = enRes * sqrt(pow(parameters[2], 2)/pow(parameters[0], 2) +
                            pow(parameters[3], 2)/pow(parameters[1], 2));
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
      fEnergyResCh0 = enResAve;
      fEnergyResCh0Err = enResAveErr;
      fPeaksCh0 = peaks;
  }
  else if(ch==1){
      fEnergyResGraphCh1 = graph;
      fEnergyResCh1 = enResAve;
      fEnergyResCh1Err = enResAveErr;
      fPeaksCh1 = peaks;
  }
            
  return true;
}
//------------------------------------------------------------------
bool SFEnergyResolution::CalculateEnergyRes(void){
    
  std::cout << "\n----- Inside SFEnergyResolution::CalculateEnergyRes()" << std::endl;
  std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <SFPeakFinder*> peakFin;
  
  SFAttenuation *att;
  
  try{
    att = new SFAttenuation(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in SFEnergyResolution::CalculateEergyRes()" << std::endl;
    std::cerr << "Couldn't access attenuation length data!" << std::endl;
    std::abort();
  }
  
  att->AttAveragedCh();
  
  TString cut_ch0 = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0";
  TString cut_ch1 = "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0";
  TString cut = cut_ch0 + std::string(" && ") + cut_ch1;
  
  std::vector <double> customNumCh0;
  customNumCh0.resize(2);
  customNumCh0[1] = att->GetAttLength();
  
  std::vector <double> customNumCh1;
  customNumCh1.resize(2);
  customNumCh1[1] = att->GetAttLength();
  
  std::vector <double> customNumSum;
  customNumSum.resize(4);
  customNumSum[1] = att->GetAttLength();
  customNumSum[3] = att->GetAttLength();
  
  TString gname = Form("ER_s%i_ave", fSeriesNo);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  std::vector <double> parameters;
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
                                                        cut_ch0, positions[i], customNumCh0));
    fSpectraCorrCh1.push_back(fData->GetCustomHistogram(1, SFSelectionType::PEAttCorrected, 
                                                        cut_ch1, positions[i], customNumCh1));
    fSpectraSum.push_back(fData->GetCustomHistogram(SFSelectionType::PEAttCorrectedSum, 
                                                    cut, positions[i], customNumSum));
    peakFin.push_back(new SFPeakFinder(fSpectraSum[i], 0));
    
    if(collimator=="Lead"){
      peakFin[i]->FindPeakNoBackground();
      fPeaksSum.push_back(peakFin[i]->GetPeak());
    }
    else if(collimator=="Electronic")
      peakFin[i]->FindPeakFit();
    
    parameters = peakFin[i]->GetParameters();
    
    enRes = parameters[1]/parameters[0];
    enResErr = enRes * sqrt(pow(parameters[2], 2)/pow(parameters[0], 2) + 
                            pow(parameters[3], 2)/pow(parameters[1], 2));
    enRes = enRes*100;
    enResErr = enResErr*100;
    graph->SetPoint(i, positions[i], enRes);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), enResErr);
    enResAve += enRes * (1./pow(enResErr, 2));
    enResAveErr += (1./pow(enResErr, 2));
  }
  
  fEnergyResSum = enResAve/enResAveErr;
  fEnergyResSumErr = sqrt(1./enResAveErr);
  fEnergyResGraphSum = graph;
  
  if(std::isnan(fEnergyResSum) || fEnergyResSum<0 || fEnergyResSum>100){
    std::cout << "##### Warning in SFEnergyResolution class!" << std::endl;
    std::cout << "Incorrect energy resolution value: " << fEnergyResSum  
              << " +/- " << fEnergyResSumErr << std::endl;
    fEnergyResSum=0;
    fEnergyResSumErr=0;
  }
  
  std::cout << "Average energy resolution calculated from summed and attenuation length corrected histograms is: "  << fEnergyResSum << " +/- " << fEnergyResSumErr << " % \n" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
std::vector <double> SFEnergyResolution::GetEnergyResolution(int ch){
    
  std::vector <double> temp;
  if(ch==0){
    temp.push_back(fEnergyResCh0);
    temp.push_back(fEnergyResCh0Err);
  }
  else if(ch==1){
    temp.push_back(fEnergyResCh1);
    temp.push_back(fEnergyResCh1Err);
  }
  else{
    std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolution() for ch" << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolution() for ch" << ch << std::endl;
    std::cerr << "Incorrect values: " << temp[0] << " +/- " << temp[1] << " %" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFEnergyResolution::GetEnergyResolution(void){
  
  std::vector <double> temp;
  temp.push_back(fEnergyResSum);
  temp.push_back(fEnergyResSumErr);
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolution()" << std::endl;
    std::cerr << "Incorrect values: " << temp[0] << " +/- " << temp[1] << " %" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyResolution::GetEnergyResolutionGraph(int ch){
    
  if((ch==0 && fEnergyResGraphCh0==nullptr) || 
     (ch==1 && fEnergyResGraphCh1==nullptr)) {
      std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolutionGraph() for ch" 
                << ch << std::endl;
      std::cerr << "Requested graph doesnt exist!" << std::endl;
      std::abort();
  }
  
  if(ch==0)
    return fEnergyResGraphCh0;
  else if(ch==1) 
    return fEnergyResGraphCh1;
  else{
    std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolutionGraph() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
TGraphErrors* SFEnergyResolution::GetEnergyResolutionGraph(void){
 
  if(fEnergyResGraphSum==nullptr){
    std::cerr << "##### Error in SFEnergyResolution::GetEnergyResolution()" << std::endl;  
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();
  }
  
  return fEnergyResGraphSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyResolution::GetSpectra(int ch){
    
  if((ch==0 && fSpectraCh0.empty()) ||
     (ch==1 && fSpectraCh1.empty())) {
      std::cerr << "##### Error in SFEnergyResolution::GetSpectra() fo ch" << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fSpectraCh0;
  else if(ch==1)
    return fSpectraCh1;
  else{
    std::cerr << "##### Error in SFEnergyResolution::GetSpectra for ch" << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyResolution::GetSpectraCorrected(int ch){
    
  if((ch==0 && fSpectraCorrCh0.empty()) ||
     (ch==1 && fSpectraCorrCh1.empty())) {
      std::cerr << "##### Error in SFEnergyResolution::GetSpectraCorrected() for ch" 
                << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fSpectraCorrCh0;
  else if(ch==1)
    return fSpectraCorrCh1;
  else{
    std::cerr << "##### Error in SFEnergyResolution::GetSpectraCorrected() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyResolution::GetSpectraSum(void){
  
  if(fSpectraSum.empty()){
    std::cerr << "##### Error in SFEnergyResolution::GetSpectraSum()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fSpectraSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyResolution::GetPeaks(void){
  
  if(fPeaksSum.empty()){
    std::cerr << "##### Error in SFEnergyResolution::GetPeaks()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fPeaksSum;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFEnergyResolution::GetPeaks(int ch){
    
  if((ch==0 && fPeaksCh0.empty()) ||
     (ch==1 && fPeaksCh1.empty())) {
      std::cerr << "##### Error in SFEnergyResolution::GetPeaks() for ch" 
                << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fPeaksCh0;
  else if(ch==1)
    return fPeaksCh1;
  else{
    std::cerr << "##### Error in SFEnergyResolution::GetPeaks() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
/// Prints details of the SFEnergyResolution class object.
void SFEnergyResolution::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl;
  std::cout << "This is print out of SFEnergyResolution class object" << std::endl;
  std::cout << "Experimental series number " << fSeriesNo << std::endl;
  std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
