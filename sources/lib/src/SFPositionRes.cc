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
                                            fData(nullptr),
                                            fAtt(nullptr),
                                            fPosRecoVsPos(nullptr),
                                            fPosResVsPos(nullptr),
                                            fMLRvsPos(nullptr) {
                                                
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
  
  fSpecAv  = fData->GetCustomHistograms(SFSelectionType::PEAverage, 
                                        "ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fPE>0 && ch_1.fPE>0");
}
//------------------------------------------------------------------
SFPositionRes::~SFPositionRes(){
    
  if(fAtt != nullptr) delete fAtt;
  if(fData != nullptr) delete fData;
}
//------------------------------------------------------------------
bool SFPositionRes::AnalyzePositionRes(void){
 
  std::cout << "\n\n----- Position Resolution Analysis" << std::endl;
  std::cout << "----- Series: " << fSeriesNo << std::endl;  
  
  double xmin, xmax;
  TString cut = "";
  
  int npoints = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measurementsIDs = fData->GetMeasurementsIDs();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  
  std::vector <SFPeakFinder*> peakFinCh0;
  std::vector <SFPeakFinder*> peakFinCh1;
  std::vector <SFPeakFinder*> peakFinAv;
    
  fPosRecoVsPos = new TGraphErrors(npoints);
  fPosRecoVsPos->SetMarkerStyle(4);
  fPosRecoVsPos->GetXaxis()->SetTitle("source position [mm]");
  fPosRecoVsPos->GetYaxis()->SetTitle("reconstructed source position [mm]");
  fPosRecoVsPos->SetName(Form("PosRecoVsPos_S%i", fSeriesNo)); 
  fPosRecoVsPos->SetTitle(Form("Reconstructed source position S%i", fSeriesNo));
  
  fPosResVsPos = new TGraphErrors(npoints);
  fPosResVsPos->SetMarkerStyle(4);
  fPosResVsPos->GetXaxis()->SetTitle("source position [mm]");
  fPosResVsPos->GetYaxis()->SetTitle("position resolution [mm]");
  fPosResVsPos->SetName(Form("MLRPosResVsPos_S%i", fSeriesNo));
  fPosResVsPos->SetTitle(Form("Position resolution S%i", fSeriesNo));
  
  double mean, sigma;
  double meanErr, sigmaErr;
  double MLR, pos;
  double posResAv = 0;
  double posResAvErr = 0;
  
  fAtt->AttAveragedCh();
  fAtt->Fit3rdOrder();
  fMLRvsPos = fAtt->GetAttGraph();
  TF1 *fPol3 = (TF1*)fMLRvsPos->GetFunction("fpol3");
  
  if(fPol3 == nullptr){
    std::cerr << "##### Error in SFPositionRes::AnalyzePositionRes()" << std::endl;
    std::cerr << "Attenuation function was not found!" << std::endl;
    return false;
  }
  
  std::vector <TF1*> funGaus;
  std::vector <TTree*> trees;
  std::vector <double> FWHM;
  int nentries = 0;
  
  DDSignal *sig_ch0 = new DDSignal();
  DDSignal *sig_ch1 = new DDSignal();
  
  for(int i=0; i<npoints; i++){ 
    
    std::cout << "\t Analyzing position " << positions[i] << " mm..." << std::endl;  
      
    //----- geting tree
    trees.push_back(fData->GetTree(measurementsIDs[i]));
    nentries = trees[i]->GetEntries();
    trees[i]->SetBranchAddress("ch_0", &sig_ch0);
    trees[i]->SetBranchAddress("ch_1", &sig_ch1);

    //----- setting energy cut
    peakFinAv.push_back(new SFPeakFinder(fSpecAv[i], false));
    peakFinAv[i]->FindPeakRange(xmin, xmax);

    //----- setting histogram
    TString hname = Form("hPosReco_S%i_pos%.1f", fSeriesNo, positions[i]);
    fPosRecoDist.push_back(new TH1D(hname, hname, 500, -50, 150));
    
    //----- filling histogram
    for(int ii=0; ii<nentries; ii++){
      trees[i]->GetEntry(ii);
      if(sig_ch0->GetT0()>0 && sig_ch1->GetT0()>0 &&
         sqrt(sig_ch0->GetPE()*sig_ch1->GetPE())>xmin &&
         sqrt(sig_ch0->GetPE()*sig_ch1->GetPE())<xmax){  
        MLR = log(sqrt(sig_ch1->GetPE()/sig_ch0->GetPE()));
        pos = fPol3->GetX(MLR);
        fPosRecoDist[i]->Fill(pos);
      }
    }
    
    //----- fitting histogram and calculating position resolution
    mean  = fPosRecoDist[i]->GetMean();
    sigma = fPosRecoDist[i]->GetRMS();
    funGaus.push_back(new TF1("funGaus", "gaus", mean-5*sigma, mean+5*sigma));
    fPosRecoDist[i]->Fit(funGaus[i], "QR");
    mean     = funGaus[i]->GetParameter(1);
    meanErr  = funGaus[i]->GetParError(1);
    FWHM = SFTools::GetFWHM(fPosRecoDist[i]);
    
    fResults.fPosReco.push_back(mean);
    fResults.fPosRecoErr.push_back(meanErr);
    fPosRecoVsPos->SetPoint(i, positions[i], fResults.fPosReco[i]);
    fPosRecoVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), fResults.fPosRecoErr[i]);
    
    fResults.fPosResAll.push_back(FWHM[0]);
    fResults.fPosResAllErr.push_back(FWHM[1]);
    fPosResVsPos->SetPoint(i, positions[i], fResults.fPosResAll[i]);
    fPosResVsPos->SetPointError(i, SFTools::GetPosError(collimator, testBench), fResults.fPosResAllErr[i]);
    
    posResAv += fResults.fPosResAll[i] * (1./pow(fResults.fPosResAllErr[i], 2));
    posResAvErr += (1./pow(fResults.fPosResAllErr[i], 2));
  }
  
  fResults.fPosRes = posResAv/posResAvErr;
  fResults.fPosResErr = sqrt(1./posResAvErr);
  
  std::cout << "Average position resolution for this series is: ";
  std::cout << fResults.fPosRes << " +/- " << fResults.fPosResErr << " mm\n\n"  << std::endl;
    
  return true;  
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetPositionRecoGraph(void){
  
  if(fPosRecoVsPos == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetPosRecoGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  
  return fPosRecoVsPos;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetPositionResGraph(void){
  
  if(fPosResVsPos == nullptr){
    std::cerr << "##### Error in SFPositionRes::GetPositionResGraph()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  
  return fPosResVsPos;
}
//------------------------------------------------------------------
TGraphErrors* SFPositionRes::GetAttenuationCurve(void){
  
  if(fMLRvsPos==nullptr){
    std::cerr << "##### Error in SFPositionRes::GetAttenuationCurve()!" << std::endl; 
    std::cerr << "Requested graph doesn't exist!" << std::endl;  
    std::abort();
  }
  
  return fMLRvsPos;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFPositionRes::GetPositionRecoDist(void){
  
  if(fPosRecoDist.empty()){
    std::cerr << "##### Error in SFPositionRes::GetPositionsRecoDist()! Empty vector!" << std::endl;
    std::abort();
  }
  
  return fPosRecoDist;
}
//------------------------------------------------------------------
std::vector <TH1D*> SFPositionRes::GetSpectra(void){
  
  if(fSpecAv.empty()){
    std::cerr << "##### Error in SFPositionRes::GetSpectra()!" << std::endl;
    std::cerr << "Empty vector!" << std::endl;
    std::abort();
  }
  
  return fSpecAv;
}
//------------------------------------------------------------------
void SFPositionRes::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl; 
  std::cout << "This is print out of SFPositionRes class object" << std::endl;
  std::cout << "Experimental series number: " << fSeriesNo << std::endl;
  std::cout << "\n-------------------------------------------" << std::endl;   
}
//------------------------------------------------------------------
