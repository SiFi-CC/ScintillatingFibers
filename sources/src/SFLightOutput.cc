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
SFLightOutput::SFLightOutput(int seriesNo): fSeriesNo(seriesNo),
                                            fPDE(-1), 
                                            fCrossTalk(-1),
                                            fLightOut(-1),
                                            fLightOutErr(-1),
                                            fLightOutCh0(-1),
                                            fLightOutCh0Err(-1),
                                            fLightOutCh1(-1),
                                            fLightOutCh1Err(-1),
                                            fLightCol(-1),
                                            fLightColErr(-1),
                                            fLightColCh0(-1),
                                            fLightColCh0Err(-1),
                                            fLightColCh1(-1),
                                            fLightColCh1Err(-1),
                                            fLightOutGraph(nullptr),
                                            fLightOutCh0Graph(nullptr),
                                            fLightOutCh1Graph(nullptr),
                                            fLightColGraph(nullptr),
                                            fLightColCh0Graph(nullptr),
                                            fLightColCh1Graph(nullptr),
                                            fData(nullptr),
                                            fAtt(nullptr) {
                                                
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cout << message << std::endl;
    throw "##### Exception in SFLightOutput constructor!";
  }
  
  int npoints = fData->GetNpoints();
  TString description = fData->GetDescription();
  
  if(!description.Contains("Regular series")){
    std::cout << "##### Warning in SFLightOutput constructor!" << std::endl;
    std::cout << "Calculating light output for non-regular series!" << std::endl;
  }
  
  TString cutCh0 = "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0";
  TString cutCh1 = "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0";
  fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::PE, cutCh0);
  fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::PE, cutCh1);

  for(int i=0; i<npoints; i++){
      fPFCh0.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
      fPFCh1.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
  }
  
  fCrossTalk = GetCrossTalk();
  fPDE = GetPDE();
  
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
double SFLightOutput::GetCrossTalk(void){
    
  double crossTalk = -1;
  TString SiPM = fData->GetSiPM();
  double overvol = fData->GetOvervoltage();
  
  TString fname = std::string(getenv("SFDATA")) + 
                  "/DB/SiPMs_Data.root";
  TFile *file = new TFile(fname, "READ");
  TGraph *gCrossTalk;
  
  if(SiPM == "Hamamatsu"){
    gCrossTalk = (TGraph*)file->Get("gHamamatsu_CT");
  }
  else if(SiPM == "SensL"){
    crossTalk = 0.03;
    return crossTalk;
  }
  else{
    std::cerr << "##### Error in SFLightOutput::GetCrossTalk()!" << std::endl;
    std::cerr << "Unknown SiPM type! Please check!" << std::endl; 
    std::abort();
  }
  
  if(gCrossTalk == nullptr){
    std::cerr << "##### Error in SFLightOutput::GetCrossTalk()!" << std::endl;
    std::cerr << "Couldn't get requested cross talk graph!" << std::endl; 
    std::abort();
  }
  
  crossTalk = gCrossTalk->Eval(overvol)/100.;
  file->Close();
  
  std::cout << "\t----- Cross talk is: " << crossTalk << std::endl;
  
  return crossTalk;
}
//------------------------------------------------------------------
double SFLightOutput::GetPDE(void){
    
  double PDE = 0;
  double overvol = fData->GetOvervoltage();
  TString fiber = fData->GetFiber();
  TString SiPM = fData->GetSiPM();
    
  TString fname = std::string(getenv("SFDATA")) + 
                  "/DB/SiPMs_Data.root";
  TFile *file = new TFile(fname, "READ");

  TGraph *gPDEvsVol;
  TH1F *hPDEvsWavelen;
  TH1F *hLightOutvsWavelen;
  
  if(fiber.Contains("LuAG:Ce")){
     hLightOutvsWavelen = (TH1F*)file->Get("hLuAG_LOvsWave");
  }
  else if(fiber.Contains("LYSO:Ce")){
    hLightOutvsWavelen = (TH1F*)file->Get("hLYSO_LOvsWave");
  }
  else if(fiber.Contains("GAGG:Ce")){
    hLightOutvsWavelen = (TH1F*)file->Get("hGAGG_LOvsWave");
  }
  else{
    std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
    std::cerr << "Unknown fiber type! Please check!" << std::endl;
  }
  
  if(SiPM == "Hamamatsu"){
    gPDEvsVol = (TGraph*)file->Get("gHamamatsu_PDEvsVol");
    hPDEvsWavelen = (TH1F*)file->Get("hHamamatsu_PDEvsWave");
  }
  else if(SiPM == "SensL"){
    gPDEvsVol = (TGraph*)file->Get("gSensL_PDEvsVol");
    hPDEvsWavelen = (TH1F*)file->Get("hSensL_PDEvsWave");
  }
  else{
    std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
    std::cerr << "Unknown SiPM type! Please check!" << std::endl;
  }
  
  if(gPDEvsVol == nullptr ||
     hPDEvsWavelen == nullptr ||
     hLightOutvsWavelen == nullptr){
    std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
    std::cerr << "Couldn't get one of requested graphs/histograms!" << std::endl;
    std::abort();
  }
  
  int nbins = hPDEvsWavelen->GetXaxis()->GetNbins();
  double f = gPDEvsVol->Eval(overvol)/hPDEvsWavelen->GetBinContent(hPDEvsWavelen->GetMaximumBin());

  for(int i=1; i<nbins+1; i++){
    PDE += hLightOutvsWavelen->GetBinContent(i) * hPDEvsWavelen->GetBinContent(i)/100.;
  }
  
  double PDE_final = PDE*f;
  TLatex text;

  TCanvas *can = new TCanvas("can","can",1800,1800);
  can->Divide(2,2);
  can->cd(1);
  gPDEvsVol->Draw("AP");
  can->cd(2);
  hLightOutvsWavelen->Draw("hist");
  can->cd(3);
  text.DrawLatex(0.2,0.9,Form("f = %.2f",f));
  text.DrawLatex(0.2,0.8,Form("fPDE = %.2f",gPDEvsVol->Eval(overvol)));
  text.DrawLatex(0.2,0.7,Form("fSiPM = %.2f",hPDEvsWavelen->GetBinContent(hPDEvsWavelen->GetMaximumBin())));
  text.DrawLatex(0.2,0.6,Form("Nbins emission = %i",hLightOutvsWavelen->GetNbinsX()));
  text.DrawLatex(0.2,0.5,Form("Nbins SiPM = %i",hPDEvsWavelen->GetNbinsX()));
  text.DrawLatex(0.2,0.4,Form("Conv = %.2f",PDE));
  text.DrawLatex(0.2,0.3,Form("PDE = %.2f",PDE_final));
  can->cd(4);
  hPDEvsWavelen->Draw();
  can->SaveAs(Form("lightout_S%i.root",fSeriesNo));
  
  file->Close();
  std::cout << "\t----- PDE is: " << PDE << std::endl;
  
  return PDE_final;
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
    lightOutErr = sqrt(pow(fLightOutCh0Graph->GetErrorY(i),2) + pow(fLightOutCh1Graph->GetErrorY(i),2));
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
  std::vector <double> attenuation = fAtt->GetAttLenPol1();
  std::vector <SFPeakFinder*> peakFin;
  
  for(int i=0; i<npoints; i++){
    if(ch==0)
      peakFin = fPFCh0; 
    else if(ch==1)
      peakFin = fPFCh1; 
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
    
    peakFin[i]->FindPeakFit();
    if(ch==0) distance = positions[i];
    if(ch==1) distance = 100.-positions[i];
    parameters = peakFin[i]->GetParameters();
    
    lightOut = parameters[0]*(1-fCrossTalk)/fPDE/0.511/TMath::Exp(-distance/attenuation[0]);
    lightOutErr = sqrt((pow(lightOut, 2)*pow(parameters[1], 2)/pow(parameters[0], 2)) + 
                       (pow(lightOut, 2)*pow(attenuation[1], 2)/pow(attenuation[0], 4)));
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
  }
  else if(ch==1){
    fLightOutCh1 = lightOutAv;
    fLightOutCh1Err = lightOutAvErr;
    fLightOutCh1Graph = graph;
  }
  
  std::cout << "Average light output for channel " << ch << ": " << lightOutAv 
            << " +/- " << lightOutAvErr << " ph/MeV \n" << std::endl;
  
  return true;  
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightCol(void){
  
  std::cout << "\n----- Inside SFLightOutput::CalculateLCol()" << std::endl;
  std::cout << "----- Series: " << fSeriesNo << std::endl;
  
  if(fLightColCh0Graph==nullptr || 
     fLightColCh1Graph==nullptr){
    std::cerr << "##### Error in SFLightOutput::CalculateLightCol()" << std::endl;
    std::cerr << "fLightColCh0Graph and fLightColCh1Graph don't exist!" << std::endl;
    std::cerr << "Run analysis for channels 0 and 1 first!" << std::endl;
    return false;
  }
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  
  TString gname = Form("LCol_s%i_sum",fSeriesNo);
  TGraphErrors* graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("collected light [P.E./MeV]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  double lightCol, lightColErr;
  double lightColAv, lightColAvErr;
  double x, yCh0, yCh1;
  
  for(int i=0; i<npoints; i++){
    fLightColCh0Graph->GetPoint(i, x, yCh0);
    fLightColCh1Graph->GetPoint(i, x, yCh1);
    lightCol = yCh0+yCh1;
    lightColErr = fLightColCh0Graph->GetErrorY(i) + fLightColCh1Graph->GetErrorY(i);
    graph->SetPoint(i, positions[i], lightCol);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightColErr);
    lightColAv += lightCol*(1./pow(lightColErr, 2));
    lightColAvErr += 1./pow(lightColErr, 2);
  }
  
  fLightCol = lightColAv/lightColAvErr;
  fLightColErr = sqrt(1./lightColAvErr);
  fLightColGraph = graph;
  
  std::cout << "Average light collection for this series: "
            << fLightCol << " +/- " << fLightColErr << " ph/MeV" << std::endl;
 
  return true;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightCol(int ch){
  
  std::cout << "\n----- Inside SFLightOutput::CalculateLightCol() for series " << fSeriesNo << std::endl;
  std::cout << "----- Analyzing channel " << ch << std::endl;
  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <double> positions = fData->GetPositions();
  std::vector <SFPeakFinder*> tempPF;
  
  TString gname = Form("LCol_s%i_ch%i",fSeriesNo,ch);
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("collected light [P.E./MeV]");
  graph->SetTitle(gname);
  graph->SetName(gname);
  graph->SetMarkerStyle(4);
  
  std::vector <double> peak_par;
  double lightCol = 0;
  double lightColErr = 0;
  double lightColAv = 0;
  double lightColAvErr = 0;
  
  for(int i=0; i<npoints; i++){
    if(ch==0)  
      tempPF = fPFCh0; 
    else if(ch==1)
      tempPF = fPFCh1;
    else{
      std::cerr << "##### Error in SFLightOutput::CalculateLightCol()!" << std::endl;
      std::cerr << "Incorrect channel number! PLease check!" << std::endl;
      return false;
    }
    
    peak_par = tempPF[i]->GetParameters();
    lightCol = peak_par[0]/0.511;
    lightColErr = peak_par[1]/0.511;
    graph->SetPoint(i, positions[i], lightCol);
    graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightColErr);
    
    lightColAv += lightCol*(1./pow(lightColErr, 2));
    lightColAvErr += (1./pow(lightColErr, 2));
  }
 
  lightColAv = lightColAv/lightColAvErr;
  lightColAvErr = sqrt(1./lightColAvErr);
  
  if(ch==0){
    fLightColCh0Graph = graph;
    fLightColCh0 = lightColAv;
    fLightColCh0Err = lightColAvErr;
  }
  else if(ch==1){
    fLightColCh1Graph = graph;
    fLightColCh1 = lightColAv;
    fLightColCh1Err = lightColAvErr;
  }
 
  std::cout << "Average light collection for this series and channel " << ch 
            << ": " << lightColAv << " +/- " << lightColAvErr << " ph/MeV" << std::endl;
 
  return true;
}
//------------------------------------------------------------------
std::vector <double> SFLightOutput::GetLightOutput(void){
    
  std::vector <double> temp;
  temp.push_back(fLightOut);
  temp.push_back(fLightOutErr);
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFLightOutput::GetLightOutput()" << std::endl;
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
std::vector <double> SFLightOutput::GetLightCol(void){
    
  std::vector <double> temp;
  temp.push_back(fLightCol);
  temp.push_back(fLightColErr);  
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFLightOutput::GetLightCol()" << std::endl;
    std::cerr << "Incorrect values: " << temp[0] << " +/- " << temp[1] << " ph/MeV" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFLightOutput::GetLightCol(int ch){
  
  std::vector <double> temp;
  
  if(ch==0){
    temp.push_back(fLightColCh0);
    temp.push_back(fLightColCh0Err);
  }
  else if(ch==1){
    temp.push_back(fLightColCh1);
    temp.push_back(fLightColCh1Err);
  }
  else{
    std::cerr << "##### Error in SFLightOutput::GetLightCol() for ch " << ch << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
    
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFLightOutput::GetLightCol() for ch " << ch << std::endl;
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
TGraphErrors* SFLightOutput::GetLightColGraph(void){
    
  if(fLightColGraph==nullptr){
    std::cerr << "##### Error in SFLightOUtput::GetLightColGraph()" << std::endl;
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();
  }
  
  return fLightColGraph;
}
//------------------------------------------------------------------
TGraphErrors* SFLightOutput::GetLightColGraph(int ch){
    
  if((ch==0 && fLightColCh0Graph==nullptr) ||
     (ch==1 && fLightColCh1Graph==nullptr)){
    std::cerr << "##### Error in SFLightOutput::GetLightColGraph() for ch " << ch << std::endl;
    std::cerr << "Requested graph doesn't exist!" << std::endl;
  }
    
  if(ch==0)
    return fLightColCh0Graph;
  else if(ch==1)
    return fLightColCh1Graph;
  else{
    std::cerr << "##### Error in SFLightOutput::GetLightColGraph() for ch " << ch << std::endl;
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
