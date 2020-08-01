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

const double ampMax = 660;

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
  
  double s = SFTools::GetSigmaBL(fData->GetSiPM());
  std::vector <double> sigmas = {s, s};
  TString cut = SFDrawCommands::GetCut(SFCutType::CombCh0Ch1, sigmas);
  fSpecAv  = fData->GetCustomHistograms(SFSelectionType::PEAverage, cut);
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
  
  int npointsMax = fData->GetNpoints();
  std::vector <double> positions = fData->GetPositions();
  std::vector <int> measurementsIDs = fData->GetMeasurementsIDs();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  TString sipm = fData->GetSiPM();
  
  std::vector <SFPeakFinder*> peakFinAv;
    
  fPosRecoVsPos = new TGraphErrors(npointsMax);
  fPosRecoVsPos->SetMarkerStyle(4);
  fPosRecoVsPos->GetXaxis()->SetTitle("source position [mm]");
  fPosRecoVsPos->GetYaxis()->SetTitle("reconstructed source position [mm]");
  fPosRecoVsPos->SetName(Form("PosRecoVsPos_S%i", fSeriesNo)); 
  fPosRecoVsPos->SetTitle(Form("Reconstructed source position S%i", fSeriesNo));
  
  fPosResVsPos = new TGraphErrors(npointsMax);
  fPosResVsPos->SetMarkerStyle(4);
  fPosResVsPos->GetXaxis()->SetTitle("source position [mm]");
  fPosResVsPos->GetYaxis()->SetTitle("position resolution [mm]");
  fPosResVsPos->SetName(Form("MLRPosResVsPos_S%i", fSeriesNo));
  fPosResVsPos->SetTitle(Form("Position resolution S%i", fSeriesNo));
  
  double mean, sigma;
  double meanErr;
  double MLR, pos;
  double posResAv = 0;
  double posResAvErr = 0;  
  double xmin, xmax;
  
  //-----
  fAtt->AttAveragedCh();
  fAtt->Fit3rdOrder();
  TGraphErrors *tmp = fAtt->GetAttGraph();
  
  double *x = tmp->GetX();
  double *ex = tmp->GetEX();
  double *y = tmp->GetY();
  double *ey = tmp->GetEY();
  
  fMLRvsPos = new TGraphErrors(npointsMax, y, x, ey, ex);
  fMLRvsPos->SetName("PosVsMLR");
  fMLRvsPos->SetTitle(Form("Source position vs. ln(M_{LR}) S%i", fSeriesNo));
  fMLRvsPos->GetXaxis()->SetTitle("ln(M_{LR})");
  fMLRvsPos->GetYaxis()->SetTitle("source position [mm]");
  fMLRvsPos->SetMarkerStyle(4);
  fMLRvsPos->GetXaxis()->SetRangeUser(-1,1);
  
  TF1 *funPol3 = new TF1("funpol3", "pol3", -1, 1);
  funPol3->SetParLimits(3,0,100000);
  fMLRvsPos->Fit(funPol3, "R");

  //-----
  
  if(funPol3 == nullptr){
    std::cerr << "##### Error in SFPositionRes::AnalyzePositionRes()" << std::endl;
    std::cerr << "Attenuation function was not found!" << std::endl;
    return false;
  }
  
  std::vector <TF1*> funGaus;
  std::vector <SLoop*> trees;
  std::vector <double> FWHM;
  
  double BL_sigma_cut = SFTools::GetSigmaBL(sipm);

  for(int npoint=0; npoint<npointsMax; npoint++){ 
    
    std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;  
      
    //----- geting tree
    trees.push_back(fData->GetTree(measurementsIDs[npoint]));
    int nloopMax = trees[npoint]->getEntries();
    SCategory *tSig = SCategoryManager::getCategory(SCategory::CatDDSamples);
    SCategory *tCal = SCategoryManager::getCategory(SCategory::CatFibersStackCal);
    
    //nentries = trees[npoint]->GetEntries();
    //trees[npoint]->SetBranchAddress("ch_0", &sig_ch0);
    //trees[npoint]->SetBranchAddress("ch_1", &sig_ch1);

    //----- setting energy cut
    peakFinAv.push_back(new SFPeakFinder(fSpecAv[npoint], false));
    peakFinAv[npoint]->FindPeakRange(xmin, xmax);

    //----- setting histogram
    TString hname = Form("hPosReco_S%i_pos%.1f", fSeriesNo, positions[npoint]);
    fPosRecoDist.push_back(new TH1D(hname, hname, 500, -50, 150));
    
    //----- filling histogram
    for(int nloop=0; nloop<nloopMax; ++nloop){
        
      trees[npoint]->nextEvent();
      size_t tentriesMax = tSig->getEntries();
      
      assert(tSig->getEntries() == tCal->getEntries());
        
        for (int tentries=0; tentries<tentriesMax; ++tentries){

            int m, l, f;
            SDDSamples *samples    = (SDDSamples*)tSig->getObject(tentries);
            SFibersStackCal *calib = (SFibersStackCal*)tCal->getObject(tentries);
            samples->getAddress(m, l, f);
            
            if(m == 0){
                double t0Ch0 = calib->getTimeL();
                double t0Ch1 = calib->getTimeR(); 
                double peCh0 = calib->getQDCL(); 
                double peCh1 = calib->getQDCR();
                double blCh0 = samples->getSignalL()->GetBLSigma();
                double blCh1 = samples->getSignalR()->GetBLSigma();
                double totCh0 = samples->getSignalL()->GetTOT();
                double totCh1 = samples->getSignalR()->GetTOT();
                double ampCh0 = samples->getSignalL()->GetAmplitude();
                double ampCh1 = samples->getSignalR()->GetAmplitude();
                if(t0Ch0>0 && t0Ch1>0 && totCh0>0 && totCh1>0 &&
                    blCh0<BL_sigma_cut && blCh1<BL_sigma_cut && 
                    ampCh0<ampMax && ampCh1<ampMax && 
                    sqrt(peCh0*peCh1)>xmin && sqrt(peCh0*peCh1)<xmax)
                {
                  MLR = log(sqrt(peCh1/peCh0));
                  pos = funPol3->Eval(MLR);
                  fPosRecoDist[npoint]->Fill(pos);
                }
            }
        }
    }
    
    //----- fitting histogram and calculating position resolution
    mean  = fPosRecoDist[npoint]->GetMean();
    sigma = fPosRecoDist[npoint]->GetRMS();
    double xmin = fPosRecoDist[npoint]->GetBinCenter(2);
    funGaus.push_back(new TF1("funGaus", "gaus", xmin, 200));
    
    //if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
        fPosRecoDist[npoint]->Fit(funGaus[npoint], "QR");
        mean = funGaus[npoint]->GetParameter(1);
        meanErr  = funGaus[npoint]->GetParError(1);
        if(npoint==0) FWHM.resize(2);
        FWHM[0] = 2.35*funGaus[npoint]->GetParameter(2);
        FWHM[1] = 2.35*funGaus[npoint]->GetParError(2);
    //}
    //else{
    //  fPosRecoDist[npoint]->Fit(funGaus[npoint], "Q", "", mean-5*sigma, mean+5*sigma);
    //  mean     = funGaus[npoint]->GetParameter(1);
    //  meanErr  = funGaus[npoint]->GetParError(1);
    //  FWHM = SFTools::GetFWHM(fPosRecoDist[npoint]);
    //}
    
    fResults.fPosReco.push_back(mean);
    fResults.fPosRecoErr.push_back(meanErr);
    fPosRecoVsPos->SetPoint(npoint, positions[npoint], fResults.fPosReco[npoint]);
    fPosRecoVsPos->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), fResults.fPosRecoErr[npoint]);
    
    fResults.fPosResAll.push_back(FWHM[0]);
    fResults.fPosResAllErr.push_back(FWHM[1]);
    fPosResVsPos->SetPoint(npoint, positions[npoint], fResults.fPosResAll[npoint]);
    fPosResVsPos->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), fResults.fPosResAllErr[npoint]);
    
    posResAv += fResults.fPosResAll[npoint] * (1./pow(fResults.fPosResAllErr[npoint], 2));
    posResAvErr += (1./pow(fResults.fPosResAllErr[npoint], 2));
  }

  fResults.fPosRes = posResAv/posResAvErr;
  fResults.fPosResErr = sqrt(1./posResAvErr);
  
  TF1 *funpol1 = new TF1("funpol1", "pol1", 0, 100);
  fPosRecoVsPos->Fit(funpol1, "Q");
  
  fResiduals = new TGraphErrors(npointsMax);
  fResiduals->SetName(Form("PosRecoResiduals_S%i", fSeriesNo));
  fResiduals->SetTitle(Form("Reconstructed position residuals S%i", fSeriesNo));
  fResiduals->GetXaxis()->SetTitle("source position [mm]");
  fResiduals->GetYaxis()->SetTitle("residual [mm]");
  fResiduals->SetMarkerStyle(4);
  
  double res;
  double point_x, point_y;
  
  for(int npoint=0; npoint<npointsMax; npoint++){
    fPosRecoVsPos->GetPoint(npoint, point_x, point_y);
    res = point_y - funpol1->Eval(point_x);
    fResiduals->SetPoint(npoint, point_x, res);
  }
  
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
TGraphErrors* SFPositionRes::GetResiduals(void){
    
  if(fResiduals == nullptr){
      std::cerr << "#### Error in SFPositionRes::GetResiduals()!" << std::endl;
      std::cerr << "Requested graph doesn't exist!" << std::endl;
      std::abort();
  }
   
   return fResiduals;
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
