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
                                              fCh1ResGraph(nullptr),
                                              fCh0StdDev(-1),
                                              fCh1StdDev(-1),
                                              fCh0Mean(-1),
                                              fCh1Mean(-1){
  
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
  
  fSpecCh0 = fData->GetSpectra(0, SFSelectionType::PE, "ch_0.fPE>0 && ch_0.fT0>0");
  fSpecCh1 = fData->GetSpectra(1, SFSelectionType::PE, "ch_1.fPE>0 && ch_1.fT0>0");
                                                  
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
    
  int npoints          = fData->GetNpoints();
    
  std::vector <SFPeakFinder*> peakFin;
  std::vector <TH1D*> spec;
  std::vector <double> peakParams;
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
  gPeakPos->SetTitle(Form("511PeakPosition_S%i_ch%i", fSeriesNo, ch));
  gPeakPos->GetXaxis()->SetTitle("source position [mm]");
  gPeakPos->GetYaxis()->SetTitle("511 keV peak position [PE]");
  gPeakPos->SetMarkerStyle(4);
  
  TGraphErrors *gResiduals = new TGraphErrors(npoints);  
  gResiduals->SetTitle(Form("Residuals_S%i_ch%i", fSeriesNo, ch));
  gResiduals->GetXaxis()->SetTitle("source position [mm]");
  gResiduals->GetYaxis()->SetTitle("residual [PE]");
  gResiduals->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    peakFin.push_back(new SFPeakFinder(spec[i], false));
    peakFin[i]->FindPeakFit();
    peakParams = peakFin[i]->GetParameters();
    gPeakPos->SetPoint(i, i, peakParams[0]);
    gPeakPos->SetPointError(i, 0, peakParams[2]);
    peakPositions.push_back(peakParams[0]);
  }
  
  TF1 *funPol0 = new TF1("funPol0", "pol0", 0, npoints);
  gPeakPos->Fit(funPol0, "Q");
    
  double mean = SFTools::GetMean(peakPositions);
  double stdDev = SFTools::GetStandardDev(peakPositions);
  
  std::cout << "\n\tFit parameters: " << funPol0->GetParameter(0) <<  "+/-"  
            << funPol0->GetParError(0) << std::endl;
  std::cout << "\tCalculation results: " << mean << "+/-" << stdDev << std::endl;
  
  double x, y;
  double res, resErr;
  
  for(int i=0; i<npoints; i++){
    gPeakPos->GetPoint(i, x, y);
    res = y - funPol0->GetParameter(0);
    resErr = gPeakPos->GetErrorY(i) + funPol0->GetParError(0);
    gResiduals->SetPoint(i, i, res);
    gResiduals->SetPointError(i, 0, resErr);
  }
  
  if(ch==0){
    fCh0Graph    = gPeakPos;
    fCh0ResGraph = gResiduals;
    fCh0Mean     = mean;
    fCh0StdDev   = stdDev;
  }
  else if(ch==1){
    fCh1Graph    = gPeakPos;
    fCh1ResGraph = gResiduals;
    fCh1Mean     = mean;
    fCh1StdDev   = stdDev;
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
double SFStabilityMon::GetStdDev(int ch){
    
  double stdDev = -1;
  
  if(ch==0)
    stdDev = fCh0StdDev;
  else if(ch==1)
    stdDev = fCh1StdDev;
  else{
    std::cerr << "##### Error in SFStabilityMon::GetStdDev()" << std::endl;  
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    std::abort();
  }
  
  if(fabs(stdDev+1)<1E-10){
    std::cerr << "##### Error in SFStabilityMon::GetStdDev()" << std::endl;  
    std::cerr << "Incorrect standard deviation value! Please check!" << std::endl;  
    std::abort();
  }
    
  return stdDev;
}
//------------------------------------------------------------------
double SFStabilityMon::GetMean(int ch){
    
  double mean = -1;
  
  if(ch==0)
    mean = fCh0Mean;
  else if(ch==1)
    mean = fCh1Mean;
  else{
    std::cerr << "##### Error in SFStabilityMon::GetMean()" << std::endl;  
    std::cerr << "Incorrect channel number! Please check!" << std::endl;  
    std::abort();
  }
  
  if(fabs(mean+1)<1E-10){
    std::cerr << "##### Error in SFStabilityMon::GetMean()" << std::endl;  
    std::cerr << "Incorrect mean value! Please check!" << std::endl;  
    std::abort();
  }
    
  return mean;
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
