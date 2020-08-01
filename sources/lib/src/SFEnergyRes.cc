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
                                        fEnergyResGraphAve(nullptr) {

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
  double s = SFTools::GetSigmaBL(fData->GetSiPM());
  std::vector <double> sigmas = {s, s};
  TString cut_ch0    = SFDrawCommands::GetCut(SFCutType::SpecCh0, sigmas);
  TString cut_ch1    = SFDrawCommands::GetCut(SFCutType::SpecCh1, sigmas);
  TString cut_ch0ch1 = SFDrawCommands::GetCut(SFCutType::CombCh0Ch1, sigmas);
  fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::PE, cut_ch0);
  fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::PE, cut_ch1);
  fSpectraAve = fData->GetCustomHistograms(SFSelectionType::PEAverage, cut_ch0ch1);
   
}
//------------------------------------------------------------------
/// Default destructor.
SFEnergyRes::~SFEnergyRes(){
  
  if(fData!=nullptr)
    delete fData;
}
//------------------------------------------------------------------
/// Calculates energy resolution based on charge spectra of requested channel 
/// for all measurements in analyzed measurement series. In this function
/// fEnergyResGrapCh0 and fEnergyResGraphCh1 graphs are filled and values of
/// fResults.fEnergyResCh0, fEnergyResCh0Err, fResults.fEnergyResCh1 and
/// fResults.fEnergyResAveCh1Err are assigned. 
/// \param ch - channel number
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
  
  SFPeakParams parameters;
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
/// Calculates energy resolution based on averaged charge spectra for 
/// all measurements in analyzed measurement series. In this function
/// fEnergyResGraphAve graph is filled and values of fResults.fEnergyResAve
/// and fResults.fEnergyResAveErr are assigned.  
bool SFEnergyRes::CalculateEnergyRes(void){
    
  std::cout << "\n----- Inside SFEnergyRes::CalculateEnergyRes()" << std::endl;
  std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measIDs = fData->GetMeasurementsIDs();
  std::vector <SFPeakFinder*> peakFin;
  
  TString gname = Form("ER_s%i_ave", fSeriesNo);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("energy resolution [%]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  SFPeakParams parameters;
  double enRes, enResErr;
  double enResAve, enResAveErr;
  
  for(int i=0; i<npoints; i++){    
    peakFin.push_back(new SFPeakFinder(fSpectraAve[i], 0));
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
  
  fResults.fEnergyResAve = enResAve/enResAveErr;
  fResults.fEnergyResAveErr = sqrt(1./enResAveErr);
  fEnergyResGraphAve = graph;
  
  if(std::isnan(fResults.fEnergyResAve) || fResults.fEnergyResAve<0 || fResults.fEnergyResAve>100){
    std::cout << "##### Warning in SFEnergyRes class!" << std::endl;
    std::cout << "Incorrect energy resolution value: " << fResults.fEnergyResAve  
              << " +/- " << fResults.fEnergyResAveErr << std::endl;
    fResults.fEnergyResAve = 0;
    fResults.fEnergyResAveErr = 0;
  }
  
  std::cout << "Average energy resolution calculated from summed and attenuation length corrected histograms is: "  << fResults.fEnergyResAve << " +/- " << fResults.fEnergyResAveErr << " % \n" << std::endl;
  
  return true;
}
//------------------------------------------------------------------
/// Returns energy resolution graph for requested channel. Energy resolution
/// graph shows dependency of energy resolution [%] determined for given source 
/// position vs. the position [mm]
/// \param ch - channel number
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
/// Returns energy resolution graph for averaged channels. Energy resolution 
/// graph shows dependency of energy resolution [%] determined for given source
/// position vs. the position [mm].
TGraphErrors* SFEnergyRes::GetEnergyResolutionGraph(void){
 
  if(fEnergyResGraphAve==nullptr){
    std::cerr << "##### Error in SFEnergyRes::GetEnergyResolution()" << std::endl;  
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();
  }
  
  return fEnergyResGraphAve;
}
//------------------------------------------------------------------
/// Returns vector containing charge spectra of requested channel.
/// \par ch - channel number 
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
/// Returns vector containing averaged charge spectra.
std::vector <TH1D*> SFEnergyRes::GetSpectra(void){
  
  if(fSpectraAve.empty()){
    std::cerr << "##### Error in SFEnergyRes::GetSpectraSum()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fSpectraAve;
}
//------------------------------------------------------------------
/// Returns vector containing background subtracted averaged charge spectra.
std::vector <TH1D*> SFEnergyRes::GetPeaks(void){
  
  if(fPeaksAve.empty()){
    std::cerr << "##### Error in SFEnergyRes::GetPeaks()" << std::endl;
    std::cerr << "No spectra available!" << std::endl;
    std::abort();
  }
  
  return fPeaksAve;
}
//------------------------------------------------------------------
/// Returns vector containing background subtracted charge spectra of requested 
/// channel.
/// \param ch - channel number
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
