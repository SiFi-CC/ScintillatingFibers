// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionRes.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// ***************************************** 

#include "SFPositionRes.hh" 

ClassImp(SFPositionRes);

//------------------------------------------------------------------
SFPositionRes::SFPositionRes(int seriesNo): fSeriesNo(seriesNo),
                                            fPosResSeries(-1),
                                            fPosResSeriesErr(-1),
                                            fData(nullptr),
                                            fPosVsMean(nullptr),
                                            fPosVsSigmaAtt(nullptr),
                                            fRecoPosition(nullptr) {
                                                
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFPositionRes constructor!";
  }
  
  TString desc = fData->GetDescription();
  
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in SFPositionRes constructor! Non-regular series!" << std::endl;  
    throw "##### Exception in SFPositionRes constructor!";
  }
  
  
  try{
    fAtt = new SFAttenuation(fSeriesNo);  
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Error in SFPositionRes constructor! Problem with SFAttenuation!" << std::endl;
    throw "##### Exception in SFPositionRes constructor!";
  }
  
  fAtt->AttAveragedCh();
}
//------------------------------------------------------------------
SFPositionRes::~SFPositionRes(){
}
//------------------------------------------------------------------
bool SFPositionRes::LoadRatios(void){
 
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  
  double mean, sigma;
  double xminCh0, xmaxCh0;
  double xminCh1, xmaxCh1;
  TString cut = "";
  
  TF1 *fun = new TF1("fun", "gaus", -100, 100);
  std::vector <SFPeakFinder*> peakFinCh0;
  std::vector <SFPeakFinder*> peakFinCh1;
  fSpecCh0 = fData->GetSpectra(0, SFSelectionType::PE, 
                               "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  fSpecCh1 = fData->GetSpectra(1, SFSelectionType::PE, 
                               "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0");
  
  for(int i=0; i<npoints; i++){
    //fSpecCh0[i]->Print();
    //fSpecCh1[i]->Print();
    peakFinCh0.push_back(new SFPeakFinder(fSpecCh0[i], false));
    peakFinCh1.push_back(new SFPeakFinder(fSpecCh1[i], false));
    peakFinCh0[i]->FindPeakRange(xminCh0, xmaxCh0);
    peakFinCh1[i]->FindPeakRange(xminCh1, xmaxCh1);
    //std::cout << xminCh0 << " " << xmaxCh0 << std::endl;
    //std::cout << xminCh1 << " " << xmaxCh1 << std::endl;
    cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f", xminCh0, xmaxCh0, xminCh1, xmaxCh1);
    fQRatios.push_back(fData->GetCustomHistogram(SFSelectionType::LogSqrtPERatio, cut, positions[i]));
    mean = fQRatios[i]->GetMean();
    sigma = fQRatios[i]->GetRMS();
    //std::cout << mean << " " << sigma << std::endl;
    //fQRatios[i]->Print();
    fQRatios[i]->Fit(fun, "Q", "", mean-5*sigma, mean+5*sigma);
  }
  
  return true;  
}
//------------------------------------------------------------------
bool SFPositionRes::AnalyzePositionRes(void){
 
  std::cout << "\n\n----- Position resolution analysis for series " << fSeriesNo << std::endl;  
  
  if(fQRatios.empty())
    LoadRatios();  
    
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  fPosVsMean = new TGraphErrors(npoints);
  fPosVsMean->SetMarkerStyle(4);
  fPosVsMean->GetXaxis()->SetTitle("source position [mm]");
  fPosVsMean->GetYaxis()->SetTitle("mean of #sqrt{Q1/Q0} distribution");
  fPosVsMean->SetName(Form("PosVsMean_S%i", fSeriesNo)); 
  fPosVsMean->SetTitle(Form("Mean of #sqrt{Q1/Q0} S%i", fSeriesNo));
  
  fPosVsSigmaAtt = new TGraphErrors(npoints);
  fPosVsSigmaAtt->SetMarkerStyle(4);
  fPosVsSigmaAtt->GetXaxis()->SetTitle("source position [mm]");
  fPosVsSigmaAtt->GetYaxis()->SetTitle("position resolution [mm]");
  fPosVsSigmaAtt->SetName(Form("PosVsSigmaAtt_S%i", fSeriesNo));
  fPosVsSigmaAtt->SetTitle(Form("Position resolution S%i", fSeriesNo));
  
  double mean, sigma;
  double meanErr, sigmaErr;
  double posResAv = 0;
  double posResAvErr = 0;
  std::vector <double> attlen = fAtt->GetAttenuation();
  
  for(int i=0; i<npoints; i++){
    mean = fQRatios[i]->GetFunction("fun")->GetParameter(1);
    meanErr = fQRatios[i]->GetFunction("fun")->GetParError(1);
    sigma = fQRatios[i]->GetFunction("fun")->GetParameter(2);
    sigmaErr = fQRatios[i]->GetFunction("fun")->GetParError(2);
    fPosVsMean->SetPoint(i, positions[i], mean);
    fPosVsMean->SetPointError(i, SFTools::GetPosError(collimator, testBench), meanErr);
    fPosRes.push_back(sigma*attlen[0]);
    fPosResErr.push_back(attlen[0]*sigmaErr + sigma*attlen[1]);
    fPosVsSigmaAtt->SetPoint(i, positions[i], fPosRes[i]);
    fPosVsSigmaAtt->SetPointError(i, SFTools::GetPosError(collimator, testBench), fPosResErr[i]);
    posResAv += fPosRes[i] * (1./pow(fPosResErr[i], 2));
    posResAvErr += (1./pow(fPosResErr[i], 2));
  }
  
  fPosResSeries = posResAv/posResAvErr;
  fPosResSeriesErr = sqrt(1./posResAvErr);
  
  std::cout << "Average position resolution for this series is: ";
  std::cout << fPosResSeries << " +/- " << fPosResSeriesErr << " mm"  << std::endl;
    
  return true;  
}
//------------------------------------------------------------------
bool SFPositionRes::ReconstructPos(void){
    
  std::cout << "\n\n ----- Position reconstruction for series " << fSeriesNo << std::endl;
  
  if(fQRatios.empty())
    LoadRatios();
  
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  fRecoPosition = new TGraphErrors(npoints);
  fRecoPosition->SetMarkerStyle(4);
  fRecoPosition->GetXaxis()->SetTitle("true source position [mm]");
  fRecoPosition->GetYaxis()->SetTitle("reconstructed source position [mm]");
  fRecoPosition->SetName(Form("RecoPos_S%i", fSeriesNo)); 
  fRecoPosition->SetTitle(Form("Reconstructed source position S%i", fSeriesNo));
  
  double mean, meanErr;
  double recoPos, recoPosErr;
  std::vector <double> attlen = fAtt->GetAttenuation();
  std::vector <double> a0 = fAtt->GetA0();
  
  for(int i=0; i<npoints; i++){
    mean = fQRatios[i]->GetFunction("fun")->GetParameter(1);
    meanErr = fQRatios[i]->GetFunction("fun")->GetParError(1);
    recoPos = attlen[0] * (mean-a0[0]);
    recoPosErr = attlen[1]*(mean-a0[0]) + attlen[0]*(meanErr+a0[1]);
    fRecoPosition->SetPoint(i, positions[i], recoPos);
    fRecoPosition->SetPointError(i, SFTools::GetPosError(collimator, testBench), recoPosErr);
  }
  
  return true;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetMeanGraph(void){
    
  if(fPosVsMean == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetMeanGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  return fPosVsMean;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetPositionResGraph(void){
    
  if(fPosVsSigmaAtt == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetSigmaGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  return fPosVsSigmaAtt;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetRecoPosGraph(void){

  if(fRecoPosition==nullptr){
    std::cerr << "##### Error in SFPositionRes::GetRecoPosGraph()!" << std::endl;
    std::cerr << "Requested graph doesn't exist!" << std::endl;
    std::abort();    
  }
  return fRecoPosition;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFPositionRes::GetRatios(void){

  if(fQRatios.empty()){
    std::cerr << "##### Error in SFPositionRes::GetRatios()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return fQRatios;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFPositionRes::GetSpectra(int ch){
  
  std::vector <TH1D*> temp;

  if(ch==0)
    temp = fSpecCh0;
  else if(ch==1)
    temp = fSpecCh1;
  else{
    std::cerr << "##### Error in SFPositionRes::GetSpectra()!" << std::endl;
    std::cerr << "Incorrect channel number!" << std::endl;
    std::abort();
  }
  
  if(temp.empty()){
    std::cerr << "##### Error in SFPositionRes::GetSpectra()!" << std::endl;
    std::cerr << "Empty vector!" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFPositionRes::GetPositionResSeries(void){

  std::vector <double> temp;
  temp.push_back(fPosResSeries);
  temp.push_back(fPosResSeriesErr);
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFPositionRes::GetPositionResSeries()!";
    std::cerr << " Incorrect results!" << std::endl;
    std::cerr << "PosRes = " << temp[0] << " +/- " << temp[1] << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFPositionRes::GetPositionRes(void){

  if(fPosRes.empty()){
    std::cerr << "##### Error in SFPositionRes::GetPositionRes()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return fPosRes;
}
//------------------------------------------------------------------
std::vector <double> SFPositionRes::GetPositionResError(void){

  if(fPosResErr.empty()){
    std::cerr << "##### Error in SFPositionRes::GetPositionResError()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return fPosResErr;
}
//------------------------------------------------------------------
void SFPositionRes::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl; 
  std::cout << "This is print out of SFPositionRes class object" << std::endl;
  std::cout << "Experimental series number: " << fSeriesNo << std::endl;
  std::cout << "\n-------------------------------------------" << std::endl;   
}
//------------------------------------------------------------------
