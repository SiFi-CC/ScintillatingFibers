// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFLightOutput.cc           *
// *             Jonas Kasper              *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFLightOutput.hh"

ClassImp(SFLightOutput);

//------------------------------------------------------------------
/// Default constructor.
SFLightOutput::SFLightOutput(): fSeriesNo(-1),
                                fPDE(-1), 
                                fCrossTalk(-1),
                                fLightOut(-1),
                                fLightOutErr(-1),
                                fLightOutCh0(-1),
                                fLightOutCh0Err(-1),
                                fLightOutCh1(-1),
                                fLightOutCh1Err(-1),
                                fLightOutGraph(nullptr),
                                fLightOutCh0Graph(nullptr),
                                fLightOutCh1Graph(nullptr),
                                fData(nullptr),
                                fAtt(nullptr) {
                                    
  std::cout << "##### Warning in SFLightOutput constructor!" << std::endl;
  std::cout << "You are using default constructor!" << std::endl;
}
//------------------------------------------------------------------
SFLightOutput::SFLightOutput(int seriesNo): fSeriesNo(seriesNo),
                                            fPDE(-1), 
                                            fCrossTalk(-1),
                                            fLightOut(-1),
                                            fLightOutErr(-1),
                                            fLightOutCh0(-1),
                                            fLightOutCh0Err(-1),
                                            fLightOutCh1(-1),
                                            fLightOutCh1Err(-1),
                                            fLightOutGraph(nullptr),
                                            fLightOutCh0Graph(nullptr),
                                            fLightOutCh1Graph(nullptr),
                                            fData(nullptr),
                                            fAtt(nullptr) {
                                                
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cout << message << std::endl;
    throw "##### Exception in SFLightOutput constructor!";
  }
  
  TString description = fData->GetDescription();
  TString SiPM = fData->GetSiPM();
  
  if(!description.Contains("Regular series")){
    std::cout << "##### Warning in SFLightOutput constructor!" << std::endl;
    std::cout << "Calculating light output for non-regular series!" << std::endl;
  }
  
  TString cutCh0 = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0";
  TString cutCh1 = "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0";
  fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::PE, cutCh0);
  fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::PE, cutCh1);
  
  if(SiPM=="Hamamatsu"){
    fPDE       = 0.31; 
    fCrossTalk = 0.07; 
  }
  else if(SiPM=="SensL"){
    fPDE       = 0.4;
    fCrossTalk = 0.03;
  }
  
  try{
    fAtt = new SFAttenuation(fSeriesNo);
  }
  catch(const char *message){
    std::cout << message << std::endl;
    throw "##### Exception in SFLightOutput constructor!";
  }
  fAtt->AttAveragedCh();
}
//------------------------------------------------------------------
/// Default destructor.
SFLightOutput::~SFLightOutput(){
 
  if(fData != nullptr)
    delete fData;
  
  if(fAtt != nullptr)
    delete fAtt;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightOut(void){
    
  std::cout << "\n----- Inside SFLightOutput::CalculateLightOut()" << std::endl;
  std::cout << "----- Series: " << fSeriesNo << std::endl;
  
  if(fLightOutCh0Graph==nullptr || 
     fLightOutCh1Graph==nullptr){
    std::cerr << "##### Error in SFLightOutput::CalculateLightOut()" << std::endl;
    std::cerr << "fLightOutCh0Graph and fLightOutCh1Graph don't exist!" << std::endl;
    std::cerr << "Run analysis for channels 0 and 1 first!" << std::endl;
    return false;
  }
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  
  TString gname = Form("LightOutSum_S%i", fSeriesNo);
  
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("light output [ph/MeV]");
  graph->SetName(gname);
  graph->SetTitle(gname);
  graph->SetMarkerStyle(4);
  
  double lightOut, lightOutErr;
  double lightOutAv, lightOutAvErr;
  double x, yCh0, yCh1;
  
  for(int i=0; i<npoints; i++){
    fLightOutCh0Graph->GetPoint(i, x, yCh0);
    fLightOutCh1Graph->GetPoint(i, x, yCh1);
    lightOut = yCh0+yCh1;
    lightOutErr = sqrt((pow(fLightOutCh0Graph->GetErrorY(i), 2)) + 
                       (pow(fLightOutCh1Graph->GetErrorY(i), 2)));
    graph->SetPoint(i, positions[i], lightOut);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightOutErr);
    lightOutAv += lightOut*(1./pow(lightOutErr, 2));
    lightOutAvErr += 1./pow(lightOutErr, 2);
  }
  
  fLightOut = lightOutAv/lightOutAvErr;
  fLightOutErr = sqrt(1./lightOutAvErr);
  fLightOutGraph = graph;
  
  std::cout << "Averaged and summed light output: " << fLightOut 
            <<" +/- " << fLightOutErr << " ph/MeV" << std::endl;

  return true;  
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightOut(int ch){
    
  std::cout << "\n----- Inside SFLightOutput::CalculateLightOut()" << std::endl;
  std::cout << "----- Series: " << fSeriesNo << "\t Channel: " << ch << std::endl; 
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <double> attenuation = fAtt->GetAttenuation();
  std::vector <SFPeakFinder*> peakFin;
  std::vector <TH1D*> peaks;
  
  for(int i=0; i<npoints; i++){
    if(ch==0)
      peakFin.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
    else if(ch==1)
      peakFin.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
    else{
      std::cerr << "##### Error in SFLightOutput::CalculateLightOut!" << std::endl;
      std::cerr << "Incorrect channel number" << std::endl;
      return false;
    }
  }
  
  TString gname = Form("LightOut_S%i_Ch%i", fSeriesNo, ch);

  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("light output [ph/MeV]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  std::vector <double> parameters;
  double lightOut = 0;
  double lightOutErr = 0;
  double lightOutAv = 0;
  double lightOutAvErr = 0;
  double distance = 0;
  
  for(int i=0; i<npoints; i++){
    
    if(collimator=="Lead"){
      peakFin[i]->FindPeakNoBackground();
      peaks.push_back(peakFin[i]->GetPeak());
    }
    else if(collimator=="Electronic")
      peakFin[i]->FindPeakFit();
    
    if(ch==0) distance = positions[i];
    if(ch==1) distance = 100.-positions[i];
    parameters = peakFin[i]->GetParameters();
    
    lightOut = parameters[0]*(1-fCrossTalk)/fPDE/0.511/TMath::Exp(-distance/attenuation[0]);
    lightOutErr = sqrt((pow(lightOut, 2)*pow(parameters[1], 2)/pow(parameters[0], 2)) + 
                       (pow(lightOut, 2)*pow(attenuation[1],2)/pow(attenuation[0], 4)));
    graph->SetPoint(i, positions[i], lightOut);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightOutErr);
    lightOutAv += lightOut * (1./pow(lightOutErr, 2));
    lightOutAvErr += (1./pow(lightOutErr, 2));
  }
  
  lightOutAv = lightOutAv/lightOutAvErr;
  lightOutAvErr = sqrt(1./lightOutAvErr);
  
  if(ch==0){
    fLightOutCh0 = lightOutAv;
    fLightOutCh0Err = lightOutAvErr;
    fLightOutCh0Graph = graph;
    fPeaksCh0 = peaks;
  }
  else if(ch==1){
    fLightOutCh1 = lightOutAv;
    fLightOutCh1Err = lightOutAvErr;
    fLightOutCh1Graph = graph;
    fPeaksCh1 = peaks;
  }
  
  std::cout << "Average light output for channel " << ch << ": " << lightOutAv 
            << " +/- " << lightOutAvErr << " ph/MeV \n" << std::endl;
  
  return true;  
}
//------------------------------------------------------------------
std::vector <double> SFLightOutput::GetLightOutput(void){
    
  std::vector <double> temp;
  temp.push_back(fLightOut);
  temp.push_back(fLightOutErr);
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFLightOutout::GetLightOutput()" << std::endl;
    std::cerr << "Incorrect values: " << temp[0] << " +/- " << temp[1] << " ph/MeV" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFLightOutput::GetLightOutput(int ch){

  std::vector <double> temp;
  
  if(ch==0){
    temp.push_back(fLightOutCh0);
    temp.push_back(fLightOutCh0Err);
  }
  else if(ch==1){
    temp.push_back(fLightOutCh1);
    temp.push_back(fLightOutCh1Err);
  }
  else{
    std::cerr << "##### Error in SFLightOutput::GetLightOutput() for ch " << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
    
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFLightOutput::GetLightOutput() for ch " << ch << std::endl;
    std::cerr << "Incorrect values: " << temp[0] << " +/- " << temp[1] << " ph/MeV" << std::endl;
    std::abort();
  }
  
  return temp;  
}
//------------------------------------------------------------------
std::vector <TH1D*> SFLightOutput::GetSpectra(int ch){
    
  if((ch==0 && fSpectraCh0.empty()) ||
     (ch==1 && fSpectraCh1.empty())) {
      std::cerr << "##### Error in SFLightOutput::GetSpectra() fo ch" << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fSpectraCh0;
  else if(ch==1)
    return fSpectraCh1;
  else{
    std::cerr << "##### Error in SFLightOutput::GetSpectra for ch" << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
std::vector <TH1D*> SFLightOutput::GetPeaks(int ch){
    
  if((ch==0 && fPeaksCh0.empty()) ||
     (ch==1 && fPeaksCh1.empty())) {
      std::cerr << "##### Error in SFLightOutput::GetPeaks() for ch" 
                << ch << std::endl; 
      std::cerr << "No spectra available!" << std::endl; 
      std::abort();
  }

  if(ch==0)
    return fPeaksCh0;
  else if(ch==1)
    return fPeaksCh1;
  else{
    std::cerr << "##### Error in SFLightOutput::GetPeaks() for ch" 
              << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
TGraphErrors* SFLightOutput::GetLightOutputGraph(void){
    
  if(fLightOutGraph==nullptr){
    std::cerr << "##### Error in SFLightOut::GetLightOutputGraph()" << std::endl;
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();
  }
  
  return fLightOutGraph;
}
//------------------------------------------------------------------
TGraphErrors* SFLightOutput::GetLightOutputGraph(int ch){
  
  if((ch==0 && fLightOutCh0Graph==nullptr) ||
     (ch==1 && fLightOutCh1Graph==nullptr)){
    std::cerr << "##### Error in SFLightOutput::GetLightOutputGraph() for ch " << ch << std::endl;
    std::cerr << "Requested graph doesn't exist!" << std::endl;
  }
    
  if(ch==0)
    return fLightOutCh0Graph;
  else if(ch==1)
    return fLightOutCh1Graph;
  else{
    std::cerr << "##### Error in SFLightOutput::GetLightOutputGraph() for ch " << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
}
//------------------------------------------------------------------
/// Prints details of the SFLightOutput class object.
void SFLightOutput::Print(void){
 std::cout << "\n-------------------------------------------" << std::endl;
 std::cout << "This is print out of SFLightOutput class object" << std::endl;
 std::cout << "Experimental series number " << fSeriesNo << std::endl;
 std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
