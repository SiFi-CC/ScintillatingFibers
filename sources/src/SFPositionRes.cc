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
                                            fMLRPosResSeries(-1),
                                            fMLRPosResSeriesErr(-1),
                                            fYPosResSeries(-1),
                                            fYPosResSeriesErr(-1),
                                            fData(nullptr),
                                            fAtt(nullptr),
                                            fMLRMeanVsPos(nullptr),
                                            fMLRPosResVsPos(nullptr),
                                            fYMeanVsPos(nullptr),
                                            fYPosResVsPos(nullptr) {
                                                
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
    
  if(fAtt != nullptr) delete fAtt;
  if(fData != nullptr) delete fData;
}
//------------------------------------------------------------------
bool SFPositionRes::LoadRatios(void){
 
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  
  std::vector <double> customNum(2);
  std::vector <double> a0 = fAtt->GetA0();
  customNum[0] = fAtt->GetAttLength();
  customNum[1] = a0[0];
  
  double mean, sigma;
  double xminCh0, xmaxCh0;
  double xminCh1, xmaxCh1;
  TString cut = "";
  
  TF1 *fun_QRatio = new TF1("fun_QRatio", "gaus", -100, 100);
  TF1 *fun_QRatioCorr = new TF1("fun_QRatioCorr", "gaus", -100, 100);
  std::vector <SFPeakFinder*> peakFinCh0;
  std::vector <SFPeakFinder*> peakFinCh1;
  fSpecCh0 = fData->GetSpectra(0, SFSelectionType::PE, 
                               "ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  fSpecCh1 = fData->GetSpectra(1, SFSelectionType::PE, 
                               "ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0"); 
  for(int i=0; i<npoints; i++){
    peakFinCh0.push_back(new SFPeakFinder(fSpecCh0[i], false));
    peakFinCh1.push_back(new SFPeakFinder(fSpecCh1[i], false));
    peakFinCh0[i]->FindPeakRange(xminCh0, xmaxCh0);
    peakFinCh1[i]->FindPeakRange(xminCh1, xmaxCh1);
    cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f", xminCh0, xmaxCh0, xminCh1, xmaxCh1);
    
    fQRatios.push_back(fData->GetCustomHistogram(SFSelectionType::LogSqrtPERatio, 
                                                 cut, positions[i], customNum));
     
    mean = fQRatios[i]->GetMean();
    sigma = fQRatios[i]->GetRMS();
    fQRatios[i]->Fit(fun_QRatio, "Q", "", mean-5*sigma, mean+5*sigma);
    
    fQRatiosCorr.push_back(fData->GetCustomHistogram(SFSelectionType::MLRRatioCorrected, 
                                                     cut, positions[i], customNum));
     
    mean = fQRatiosCorr[i]->GetMean();
    sigma = fQRatiosCorr[i]->GetRMS();
    fQRatiosCorr[i]->Fit(fun_QRatioCorr, "Q", "", mean-5*sigma, mean+5*sigma);
  }
  
  return true;  
}
//------------------------------------------------------------------
bool SFPositionRes::AnalyzePositionResMLR(void){
 
  std::cout << "\n\n----- Position resolution analysis from M_LR distribution " << std::endl;
  std::cout << "----- Series: " << fSeriesNo << std::endl;  
  
  if(fQRatios.empty())
    LoadRatios();  
    
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  fMLRMeanVsPos = new TGraphErrors(npoints);
  fMLRMeanVsPos->SetMarkerStyle(4);
  fMLRMeanVsPos->GetXaxis()->SetTitle("source position [mm]");
  fMLRMeanVsPos->GetYaxis()->SetTitle("mean of log(#sqrt{Q1/Q0}) distribution");
  fMLRMeanVsPos->SetName(Form("MLRMeanVsPos_S%i", fSeriesNo)); 
  fMLRMeanVsPos->SetTitle(Form("Mean of log(#sqrt{Q1/Q0}) S%i", fSeriesNo));
  
  fMLRPosResVsPos = new TGraphErrors(npoints);
  fMLRPosResVsPos->SetMarkerStyle(4);
  fMLRPosResVsPos->GetXaxis()->SetTitle("source position [mm]");
  fMLRPosResVsPos->GetYaxis()->SetTitle("position resolution [mm]");
  fMLRPosResVsPos->SetName(Form("MLRPosResVsPos_S%i", fSeriesNo));
  fMLRPosResVsPos->SetTitle(Form("Position resolution from M_{LR} S%i", fSeriesNo));
  
  double mean, sigma;
  double meanErr, sigmaErr;
  double posResAv = 0;
  double posResAvErr = 0;
  std::vector <double> attlen = fAtt->GetAttenuation();
  
  for(int i=0; i<npoints; i++){ 
    mean = fQRatios[i]->GetFunction("fun_QRatio")->GetParameter(1);
    meanErr = fQRatios[i]->GetFunction("fun_QRatio")->GetParError(1);
    sigma = fQRatios[i]->GetFunction("fun_QRatio")->GetParameter(2);
    sigmaErr = fQRatios[i]->GetFunction("fun_QRatio")->GetParError(2);
    
    fMLRMeanVsPos->SetPoint(i, positions[i], mean);
    fMLRMeanVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), meanErr);
    
    fMLRPosRes.push_back(sigma*attlen[0]);
    fMLRPosResErr.push_back(attlen[0]*sigmaErr + sigma*attlen[1]);
    fMLRPosResVsPos->SetPoint(i, positions[i], fMLRPosRes[i]);
    fMLRPosResVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), fMLRPosResErr[i]);
    
    posResAv += fMLRPosRes[i] * (1./pow(fMLRPosResErr[i], 2));
    posResAvErr += (1./pow(fMLRPosResErr[i], 2));
  }
  
  fMLRPosResSeries = posResAv/posResAvErr;
  fMLRPosResSeriesErr = sqrt(1./posResAvErr);
  
  std::cout << "Average position resolution for this series is: ";
  std::cout << fMLRPosResSeries << " +/- " << fMLRPosResSeriesErr << " mm"  << std::endl;
    
  return true;  
}
//------------------------------------------------------------------
bool SFPositionRes::AnalyzePositionResY(void){
    
  std::cout << "\n\n ----- Position resolution analysis from Y distribution " << std::endl;
  std::cout << "----- Series: " << fSeriesNo << std::endl;
  
  if(fQRatiosCorr.empty())
    LoadRatios();
  
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  fYMeanVsPos = new TGraphErrors(npoints);
  fYMeanVsPos->SetMarkerStyle(4);
  fYMeanVsPos->GetXaxis()->SetTitle("true source position [mm]");
  fYMeanVsPos->GetYaxis()->SetTitle("reconstructed source position [mm]");
  fYMeanVsPos->SetName(Form("YMeanVsPos_S%i", fSeriesNo));
  fYMeanVsPos->SetTitle(Form("Reconstructed source position S%i", fSeriesNo));
  
  fYPosResVsPos = new TGraphErrors(npoints);
  fYPosResVsPos->SetMarkerStyle(4);
  fYPosResVsPos->GetXaxis()->SetTitle("source position [mm]");
  fYPosResVsPos->GetYaxis()->SetTitle("position resolution [mm]");
  fYPosResVsPos->SetName(Form("YPosResVsPos_S%i", fSeriesNo));
  fYPosResVsPos->SetTitle(Form("Position resolution form Y S%i", fSeriesNo));
  
  double mean, meanErr;
  double sigma, sigmaErr;
  double posResAv = 0;
  double posResAvErr = 0;
  
  std::vector <double> attlen = fAtt->GetAttenuation();
  std::vector <double> a0 = fAtt->GetA0();
  
  for(int i=0; i<npoints; i++){
    mean = fQRatiosCorr[i]->GetFunction("fun_QRatioCorr")->GetParameter(1);
    meanErr = fQRatiosCorr[i]->GetFunction("fun_QRatioCorr")->GetParError(1);
    sigma = fQRatiosCorr[i]->GetFunction("fun_QRatioCorr")->GetParameter(2);
    sigmaErr = fQRatiosCorr[i]->GetFunction("fun_QRatioCorr")->GetParError(2);
    
    fYMeanVsPos->SetPoint(i, positions[i], mean);
    fYMeanVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), meanErr);
    
    fYPosRes.push_back(sigma);
    fYPosResErr.push_back(sigmaErr);
    fYPosResVsPos->SetPoint(i, positions[i], sigma);
    fYPosResVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), sigmaErr);
    
    posResAv +=fYPosRes[i] * (1./pow(fYPosResErr[i], 2));
    posResAvErr += (1./pow(fYPosResErr[i], 2));
  }
  
  fYPosResSeries = posResAv/posResAvErr;
  fYPosResSeriesErr = sqrt(1./posResAvErr);
  
  std::cout << "Average position resolution for this series is: ";
  std::cout << fYPosResSeries << " +/- " << fYPosResSeriesErr << " mm"  << std::endl;
  
  return true;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetMeanGraph(TString opt){
  
  TGraphErrors *g;
  
  if(opt == "MLR") 
      g = fMLRMeanVsPos;
  else if(opt == "Y") 
      g = fYMeanVsPos;
  else{
      std::cerr << "##### Error in SFPositionRes::GetMeanGraph()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl; 
      std::abort();
  }
  
  if(g == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetMeanGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  
  return g;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetPositionResGraph(TString opt){
    
  TGraphErrors *g;
  
  if(opt == "MLR") 
      g = fMLRPosResVsPos;
  else if(opt == "Y") 
      g = fYPosResVsPos;
  else{
      std::cerr << "##### Error in SFPositionRes::GetPositionResGraph()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl;  
      std::abort();
  }
  
  if(g == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetSigmaGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  return g;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFPositionRes::GetRatios(TString opt){

  std::vector <TH1D*> temp;  

  if(opt == "MLR") 
      temp = fQRatios;
  else if(opt == "Y") 
      temp = fQRatiosCorr;
  else{
      std::cerr << "##### Error in SFPositionRes::GetRatios()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl;  
      std::abort();
  }
  
  if(temp.empty()){
    std::cerr << "##### Error in SFPositionRes::GetRatios()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return temp;
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
std::vector <double> SFPositionRes::GetPositionResSeries(TString opt){

  std::vector <double> temp;
  
  if(opt == "MLR"){ 
      temp.push_back(fMLRPosResSeries);
      temp.push_back(fMLRPosResSeriesErr);
  }
  else if(opt == "Y"){ 
      temp.push_back(fYPosResSeries);
      temp.push_back(fYPosResSeriesErr);
  }
  else{
      std::cerr << "##### Error in SFPositionRes::GetPositionResSeries()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl;  
      std::abort();
  }
  
  if(temp[0]==-1 || temp[1]==-1){
    std::cerr << "##### Error in SFPositionRes::GetPositionResSeries()!";
    std::cerr << " Incorrect results!" << std::endl;
    std::cerr << "PosRes = " << temp[0] << " +/- " << temp[1] << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFPositionRes::GetPositionRes(TString opt){

  std::vector <double> temp;
    
  if(opt == "MLR")
      temp = fMLRPosRes;
  else if(opt == "Y")
      temp = fYPosRes;
  else{
      std::cerr << "##### Error in SFPositionRes::GetPositionRes()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl;  
      std::abort();
  }
    
  if(temp.empty()){
    std::cerr << "##### Error in SFPositionRes::GetPositionRes()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
std::vector <double> SFPositionRes::GetPositionResError(TString opt){

  std::vector <double> temp;
    
  if(opt == "MLR")
      temp = fMLRPosResErr;
  else if(opt == "Y")
      temp = fYPosResErr;
  else{
      std::cerr << "##### Error in SFPositionRes::GetPositionResError()!" << std::endl;
      std::cerr << "Incorrect option. Possible options are: MLR or Y" << std::endl;  
      std::abort();
  }
  
  if(temp.empty()){
    std::cerr << "##### Error in SFPositionRes::GetPositionResError()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
void SFPositionRes::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl; 
  std::cout << "This is print out of SFPositionRes class object" << std::endl;
  std::cout << "Experimental series number: " << fSeriesNo << std::endl;
  std::cout << "\n-------------------------------------------" << std::endl;   
}
//------------------------------------------------------------------
