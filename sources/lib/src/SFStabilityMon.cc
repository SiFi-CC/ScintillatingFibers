// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFStabilityMon.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFStabilityMon.hh"

ClassImp(SFStabilityMon);

//------------------------------------------------------------------
SFStabilityMon::SFStabilityMon(int seriesNo): fSeriesNo(seriesNo),
                                              fData(nullptr),
                                              fCh0Graph(nullptr),
                                              fCh1Graph(nullptr),
                                              fCh0ResGraph(nullptr),
                                              fCh1ResGraph(nullptr) {
  
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;  
    throw "##### Exception in SFStabilityMon constructor!";
  }
  
  TString desc = fData->GetDescription();
  
  if(! desc.Contains("Stability monitoring")){
    std::cerr << "##### Error in SFStabilityMon constructor!" << std::endl;
    std::cerr << "Series " << fSeriesNo << " is NOT stability monitoring series!" << std::endl;
    throw "##### Exception in SFStabilityMon constructor!";
  }
  
  double s = SFTools::GetSigmaBL(fData->GetSiPM());
  std::vector <double> sigma = {s};
  TString cutCh0 = SFDrawCommands::GetCut(SFCutType::SpecCh0, sigma);
  TString cutCh1 = SFDrawCommands::GetCut(SFCutType::SpecCh1, sigma);
  
  fSpecCh0 = fData->GetSpectra(0, SFSelectionType::PE, cutCh0);
  fSpecCh1 = fData->GetSpectra(1, SFSelectionType::PE, cutCh1);
                                                  
}
//------------------------------------------------------------------
SFStabilityMon::~SFStabilityMon(){
  delete fData;
}
//------------------------------------------------------------------
bool SFStabilityMon::AnalyzeStability(int ch){
   
  std::cout << "\n\n----- Stability analysis " << std::endl; 
  std::cout << "----- Series: " << fSeriesNo << std::endl; 
  std::cout << "----- Channel: " << ch << std::endl; 
    
  int npoints = fData->GetNpoints();
    
  std::vector <SFPeakFinder*> peakFin;
  std::vector <TH1D*> spec;
  SFPeakParams  peakParams;
  std::vector <double> peakPositions;
  
  if(ch==0)
    spec = fSpecCh0;
  else if(ch==1)
    spec = fSpecCh1;
  else{
    std::cerr << "##### Error in SFStabilityMon::AnalyzeStability()!" << std::endl;  
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    return false;
  }
  
  TGraphErrors *gPeakPos = new TGraphErrors(npoints);
  gPeakPos->SetName(Form("511PeakPosition_S%i_ch%i", fSeriesNo, ch));
  gPeakPos->SetTitle(Form("511PeakPosition_S%i_ch%i", fSeriesNo, ch));
  gPeakPos->GetXaxis()->SetTitle("source position [mm]");
  gPeakPos->GetYaxis()->SetTitle("511 keV peak position [PE]");
  gPeakPos->SetMarkerStyle(4);
  
  TGraphErrors *gResiduals = new TGraphErrors(npoints);  
  gResiduals->SetName(Form("Residuals_S%i_ch%i", fSeriesNo, ch));
  gResiduals->SetTitle(Form("Residuals_S%i_ch%i", fSeriesNo, ch));
  gResiduals->GetXaxis()->SetTitle("source position [mm]");
  gResiduals->GetYaxis()->SetTitle("residual [PE]");
  gResiduals->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    peakFin.push_back(new SFPeakFinder(spec[i], false));
    peakFin[i]->FindPeakFit();
    peakParams = peakFin[i]->GetParameters();
    gPeakPos->SetPoint(i, i, peakParams.fPosition);
    gPeakPos->SetPointError(i, 0, peakParams.fPositionErr);
    peakPositions.push_back(peakParams.fPosition);
  }
    
  double mean = SFTools::GetMean(peakPositions);
  double stdDev = SFTools::GetStandardDev(peakPositions);
  
  TF1 *funPol0 = new TF1("funPol0", "pol0", 0, npoints);
  funPol0->FixParameter(0, mean);
  gPeakPos->Fit(funPol0, "Q");
  
  std::cout << "\tCalculation results: " << mean << "+/-" << stdDev << std::endl;
  
  double x, y;
  double res;
  
  for(int i=0; i<npoints; i++){
    gPeakPos->GetPoint(i, x, y);
    res = y - mean;
    gResiduals->SetPoint(i, i, res);
  }
  
  if(ch==0){
    fCh0Graph           = gPeakPos;
    fCh0ResGraph        = gResiduals;
    fResults.fCh0Mean   = mean;
    fResults.fCh0StdDev = stdDev;
  }
  else if(ch==1){
    fCh1Graph           = gPeakPos;
    fCh1ResGraph        = gResiduals;
    fResults.fCh1Mean   = mean;
    fResults.fCh1StdDev = stdDev;
  }
    
  return true;  
}
//------------------------------------------------------------------
TGraphErrors *SFStabilityMon::GetPeakPosGraph(int ch){
    
  TGraphErrors *g;

  if(ch==0)
    g = fCh0Graph;
  else if(ch==1)
    g = fCh1Graph;
  else{
    std::cerr << "##### Error in SFStabilityMon::GetPeakPosGraph()" << std::endl;  
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    std::abort();
  }
  
  if(g==nullptr){
    std::cerr << "##### Error in SFStabilityMon::GetPeakPosGraph()" << std::endl;  
    std::cerr << "Requested graph doesn't exist! Please check!" << std::endl;  
    std::abort();
  }
  
  return g;
}
//------------------------------------------------------------------
TGraphErrors *SFStabilityMon::GetResidualsGraph(int ch){
    
  TGraphErrors *g;

  if(ch==0)
    g = fCh0ResGraph;
  else if(ch==1)
    g = fCh1ResGraph;
  else{
    std::cerr << "##### Error in SFStabilityMon::GetResidualsGraph()" << std::endl;  
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    std::abort();
  }
  
  if(g==nullptr){
    std::cerr << "##### Error in SFStabilityMon::GetResidualsGraph()" << std::endl;  
    std::cerr << "Requested graph doesn't exist! Please check!" << std::endl;  
    std::abort();
  }
  
  return g;  
}
//------------------------------------------------------------------
std::vector <TH1D*> SFStabilityMon::GetSpectra(int ch){
    
  std::vector <TH1D*> tmp;
    
  if(ch==0)
    tmp = fSpecCh0;
  else if(ch==1)
    tmp = fSpecCh1;
  else{
    std::cerr << "##### Error in SFStabilityMon::GetSpectra()! " << std::endl;
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    std::abort();
  }
  
  if(tmp.empty()){
    std::cerr << "##### Error in SFStabilityMon::GetSpectra()! " << std::endl;
    std::cerr << "Requested spectra don't exist!" << std::endl; 
    std::abort(); 
  }
  
  return tmp;
}
//------------------------------------------------------------------
void SFStabilityMon::Print(void){
 std::cout << "\n-------------------------------------------" << std::endl;
 std::cout << "This is print out of SFStabilityMon class object" << std::endl;
 std::cout << "Experimental series number " << fSeriesNo << std::endl;
 std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
