// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPeakFinder.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFPeakFinder.hh"

ClassImp(SFPeakFinder);

//------------------------------------------------------------------
/// Default constructor.
SFPeakFinder::SFPeakFinder(): fSpectrum(nullptr),
                              fPeak(nullptr),
                              fPosition(-1),
                              fPosErr(-1),
                              fSigma(-1),
                              fSigErr(-1),
                              fChi2NDF(-1),
                              fVerbose(false),
                              fTests(false) {
                                  
  std::cout << "#### Warning in SFPeakFinder constructor!" << std::endl;
  std::cout << "You are using default constructor!" << std::endl;
  fVerbose = false;
}
//------------------------------------------------------------------
/// Standard constructor.
/// \param spectrum - analyzed spectrum
/// \param verbose - print-outs level
SFPeakFinder::SFPeakFinder(TH1D *spectrum, bool verbose): fSpectrum(spectrum),
                                                          fPeak(nullptr),
                                                          fPosition(-1),
                                                          fPosErr(-1),
                                                          fSigma(-1),
                                                          fSigErr(-1),
                                                          fChi2NDF(-1),
                                                          fVerbose(verbose),
                                                          fTests(false) { 
                                                              
}
//------------------------------------------------------------------
/// Standard constructor. Vebose level set by default to quiet. In order to change verbose
/// level use SetVerbLevel().
/// \param spectrum - analyzed spectrum
SFPeakFinder::SFPeakFinder(TH1D *spectrum): fSpectrum(spectrum),
                                            fPeak(nullptr),
                                            fPosition(-1),
                                            fPosErr(-1),
                                            fSigma(-1),
                                            fSigErr(-1),
                                            fChi2NDF(-1),
                                            fVerbose(false),
                                            fTests(false) {

  std::cout << "##### Warning in SFPeakFinder constructor. Quiet mode on, no print outs." << std::endl;
  std::cout << "##### To set verbose level use SetVerbLevel(bool)" << std::endl;
}
//------------------------------------------------------------------
/// Standard constructor.
/// \param spectrum - analyzed spectrum
/// \param verbose - print-outs level
/// \param tests - flag for testing mode
SFPeakFinder::SFPeakFinder(TH1D *spectrum, bool verbose, bool tests): fSpectrum(spectrum),
                                                                      fPeak(nullptr),
                                                                      fPosition(-1),
                                                                      fPosErr(-1),
                                                                      fSigma(-1),
                                                                      fSigErr(-1),
                                                                      fChi2NDF(-1),
                                                                      fVerbose(verbose),
                                                                      fTests(tests) { 
                                                              
}
//------------------------------------------------------------------
/// Default destructor.
SFPeakFinder::~SFPeakFinder(){
}
//------------------------------------------------------------------
/// If default constructor was used, sets analyzed spectrum histogram.
/// \param spectrum - histogram containing analyzed spectrum.
void SFPeakFinder::SetSpectrum(TH1D *spectrum){
  
  if(spectrum==nullptr){
    std::cerr << "##### Error in SFPeakFinder::SetSpectrum()! " << std::endl;
    std::cerr << "Spectrum cannot be null pointer!" << std::endl;
    std::abort();
  }
  fSpectrum = spectrum;
  return;
}
//------------------------------------------------------------------
/// Finds 511 keV peak in the analyzed spectrum using ROOT's TSpectrum
/// object. Simple fit of Gaussian function is then performed in narrow
/// range in order to better describe peak parameters. Results can be 
/// obtained via SFPeakFinder::GetParameters() function and other defined
/// getter functions. fPeak histogram is not filled.
bool SFPeakFinder::FindPeakSpectrum(void){

  // Searching peaks with TSpectrum
  TSpectrum *spec= new TSpectrum(10);
  int npeaks = spec->Search(fSpectrum, 20, "goff", 0.1);
  double *peaksX = spec->GetPositionX();
  double peak = TMath::MaxElement(npeaks, peaksX);

  // Setting fitting option based on verbose level
  TString opt;
  if(fVerbose) opt = "0RS";
  else opt = "Q0RS";  

  // Fitting Gaussian function to get peak position and sigma
  TF1 *fun = new TF1("fun", "gaus", peak-50, peak+50);
  TFitResultPtr fitRes = fSpectrum->Fit("fun", opt);
  int counter = 0;
  
  while(fitRes!=0 && counter<20){
   fitRes = fSpectrum->Fit("fun", opt);
   counter++;
  }
  
  if(fVerbose){
    std::cout << "Fit counter: " << counter << std::endl;
  }
  
  fun = fSpectrum->GetFunction("fun");
  
  fPosition = fun->GetParameter(1);
  fPosErr   = fun->GetParError(1);
  fSigma    = fun->GetParameter(2);
  fSigErr   = fun->GetParError(2);
  fChi2NDF  = fitRes->Chi2()/fitRes->Ndf();

  if(fPosition<0 || fSigma<0){
    std::cerr << "##### Error in SFPeakFinder::FindPeakSpectrum(). Negative peak position or sigma!" << std::endl;
    std::cerr << "position = " << fPosition << "\t sigma = " << fSigma << std::endl;
    return false;
  }
  
  return true;
}
//------------------------------------------------------------------
/// Finds range of the 511 keV peak. Range is returned as references.
/// For measurements with lead collimator range is defined as position
/// +/- sigma and for measurements with electronic collimator range is
/// defined as position +/- 2 sigma.
bool SFPeakFinder::FindPeakRange(double &min, double &max){
  
  // Getting series attributes
  int seriesNo = SFTools::GetSeriesNo(fSpectrum->GetName());
  SFData *data;
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in SFPeakFinder::FindPeakRange()!" << std::endl;
    std::abort();
  }
  
  TString type = data->GetCollimator();
  delete data;
  
  // Calculating peak range
  double pos, sigma;
  bool stat = FindPeakSpectrum();
  
  if(stat){
    pos = fPosition;
    sigma = fSigma;
  }
  else{
    return false;
  }
  
  if(type=="Lead"){
    min = pos - sigma;
    max = pos + 1.5*sigma;
  }
  else if(type=="Electronic"){
    min = pos - 2*sigma;
    max = pos + 2*sigma;
  }

  if(max<min || fabs(min+1)<1E-8 || fabs(max+1)<1E-8){
  std::cerr << "##### Error in SFPeakFinder::FindPeakRange(). Incorrect range." << std::endl;
  std::cerr << "min = " << min << "\t max = " << max << std::endl;
   return false;
  }
  
  return true;
}
//------------------------------------------------------------------
/// Fits function describing spectrum in the area where 511 keV peak
/// is visible. Fitted function: expo + pol0 + gaus. Results can be 
/// accessed via SFPeakFinder::GetParameters() function and other 
/// defined getters. fPeak histogram is not filled.
bool SFPeakFinder::FindPeakFit(void){

  // setting fitting option based on verbose level
  TString opt;
  if(fVerbose) opt = "S+";
  else opt = "SQ+";
  
  // setting fit parameters and fitting 
  double mean, sigma;
  bool stat = FindPeakSpectrum();
  
  if(stat){
    mean = fPosition;
    sigma = fSigma;
  }
  else{
    return false;   
  }
  
  double fit_min = 0;
  double fit_max = 0;
  TString name = fSpectrum->GetName();
  
  int seriesNo = SFTools::GetSeriesNo(name);
  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << "##### Error in SFPeakFinder::FindPeakFit()!" << std::endl;
    std::cerr << message << std::endl;
    return false;
  }
  TString collimator = data->GetCollimator();
  
  if(name.Contains("Sum")){
    fit_min = mean - 5*sigma;
    fit_max = mean - 4*sigma;
  }
  else{
    fit_min = mean - 3*sigma;
    fit_max = mean - 2*sigma;
  }

  TF1 *fun_expo = new TF1("fun_expo", "[0]*TMath::Exp((x-[1])*[2])", 0, 100);
  fSpectrum->Fit("fun_expo", opt, "", fit_min, fit_max);
  
  TF1 *fun_pol0 = new TF1("fun_pol0", "pol1", 0, 100);
  fSpectrum->Fit("fun_pol0", opt, "", mean+3*sigma, mean+4*sigma);
  
  if(name.Contains("Sum")){
    fit_min = mean - 3.5*sigma;
    fit_max = mean + 5*sigma;
  }
  else{
    fit_min = mean - 2.5*sigma;
    fit_max = mean + 4*sigma;
  }
  
  TF1 *fun_bg = new TF1("fun_bg", "[0]*TMath::Exp((x-[1])*[2])+pol1(3)+gaus(5)", 0, 100);
  int counter = 20;

  TFitResultPtr fitRes;
  float parlim = 100;

  printf("%f  %f  %f\n", fSpectrum->GetBinContent(fSpectrum->FindBin(mean)), mean, sigma);
  while (true) {
    fit_min = mean - 5*sigma;
    fun_bg->SetParameter(0, fun_expo->GetParameter(0));
    fun_bg->SetParameter(1, fun_expo->GetParameter(1));
    fun_bg->SetParameter(2, fun_expo->GetParameter(2));
    fun_bg->SetParameter(3, fun_pol0->GetParameter(0));
    fun_bg->SetParameter(4, fun_pol0->GetParameter(1));
    fun_bg->SetParameter(5, fSpectrum->GetBinContent(fSpectrum->FindBin(mean)));
    fun_bg->SetParameter(6, mean);
    fun_bg->SetParLimits(6, mean-2*sigma, mean+2*sigma);
    fun_bg->SetParameter(7, sigma);
    fun_bg->SetParLimits(7, 0, parlim);
    fitRes = fSpectrum->Fit("fun_bg", opt, "", fit_min, fit_max );

    if (fitRes == 0 or counter == 0) break;

    if(fVerbose) 
      std::cout << "FitStatus: " << fitRes <<  std::endl;
    --counter;
    parlim = 200;
  };
  std::cout << "Fit converged in " << (20 - counter + 1) << " passes." << std::endl;
  
  if(fVerbose){
      std::cout << "Fit counter: " << counter << std::endl;
  }
  
  fun_bg = fSpectrum->GetFunction("fun_bg");
  
  // getting peak parameters
  fPosition = fun_bg->GetParameter(6);
  fPosErr   = fun_bg->GetParError(6);
  fSigma    = fun_bg->GetParameter(7);
  fSigErr   = fun_bg->GetParError(7);
  fChi2NDF  = fitRes->Chi2()/fitRes->Ndf(); 
  
  // for tests
  if(fTests){
    opt = opt+"B";
  
    TF1 *fun_expo_clone = new TF1("fun_expo_clone","expo",fit_min,fit_max); // FIXME [0]*TMath::Exp((x-[1])*[2])
    fun_expo_clone->SetLineColor(kGreen+3);
    std::cout << fun_expo_clone->GetParameter(0) << "\t" << fun_expo_clone->GetParameter(1) << std::endl;
    fun_expo_clone->FixParameter(0,fun_bg->GetParameter(0));
    fun_expo_clone->FixParameter(1,fun_bg->GetParameter(1));
    fSpectrum->Fit("fun_expo_clone",opt);
    std::cout << fun_bg->GetParameter(0) << "\t" << fun_bg->GetParameter(1) << std::endl;
    std::cout << fun_expo_clone->GetParameter(0) << "\t" << fun_expo_clone->GetParameter(1) << std::endl;
  
    TF1 *fun_pol0_clone = new TF1("fun_pol0_clone","pol0",fit_min,fit_max);
    fun_pol0_clone->SetLineColor(kBlue);
    std::cout << fun_pol0_clone->GetParameter(0) << std::endl;
    fun_pol0_clone->FixParameter(0,fun_bg->GetParameter(2));
    fSpectrum->Fit("fun_pol0_clone",opt);
    std::cout << fun_bg->GetParameter(2) << std::endl;
    std::cout << fun_pol0_clone->GetParameter(0) << std::endl;
  
    TF1 *fun_gaus_clone = new TF1("fun_gaus_clone","gaus",fit_min,fit_max);
    fun_gaus_clone->SetLineColor(kMagenta);
    std::cout << fun_gaus_clone->GetParameter(0) << "\t" << fun_gaus_clone->GetParameter(1) <<  "\t" << fun_gaus_clone->GetParameter(2) <<std::endl;
    fun_gaus_clone->FixParameter(0,fun_bg->GetParameter(3));
    fun_gaus_clone->FixParameter(1,fun_bg->GetParameter(4));
    fun_gaus_clone->FixParameter(2,fun_bg->GetParameter(5));
    fSpectrum->Fit("fun_gaus_clone",opt);
    std::cout << fun_bg->GetParameter(3) << "\t" << fun_bg->GetParameter(4) <<  "\t" << fun_bg->GetParameter(5) <<std::endl;
    std::cout << fun_gaus_clone->GetParameter(0) << "\t" << fun_gaus_clone->GetParameter(1) <<  "\t" << fun_gaus_clone->GetParameter(2) <<std::endl;
   
    TFile *file = new TFile("PeakFinderTests.root","UPDATE");
    fSpectrum->Write();
    file->Close();  
  }
  
  if(fPosition<0 || fSigma<0){
    std::cerr << "##### Error in SFPeakFinder::FitPeak(). Position and Sigma cannot be negative!" << std::endl;
    std::cerr << "position = " << fPosition << "\t fSigma = " << fSigma << std::endl;
    return false;
  }

  return true;
}
//------------------------------------------------------------------
/// Finds 511 keV peak via SFPeakFinder::FindPeakSoectrum() method and
/// performs background subtraction. Exponential function is fitted on 
/// the left sige of the peak and pol0 function - on the right side. For 
/// the background subtraction functions are sewed together in the peak 
/// region. fPeak histogram is filled with the background-subtracted peak.
/// Subsequently Gaussian function is fitted in order to describe peak shape.
/// Results are accessible via SFPeakFinder::GetParameters() and other 
/// definded getter functions.
/// IMPORTANT - if you want to use this function, all the parameters of 
/// the fits need to be adjusted!
bool SFPeakFinder::FindPeakNoBackground(void){
 
  // setting background-subtracted histogram
  TString tmp = fSpectrum->GetName();
  TString pname = tmp.Append("_peak");
  fPeak = (TH1D*) fSpectrum->Clone(pname);
  fPeak->Reset();
  
  // setting fitting option based on verbose level
  TString opt;
  if(fVerbose) opt = "SR+";
  else opt = "SQR+";
  
  // setting fit parameters and fitting
  double mean, sigma;
  bool stat = FindPeakSpectrum();
  
  if(stat){
    mean = fPosition;
    sigma = fSigma;
  }
  else{
    return false;
  }
  
  double peak_min = mean-4*sigma;
  double peak_max = mean+4*sigma;
  double fit_min = 0;
  double fit_max = 0;
  TString name = fSpectrum->GetName();
  
  if(name.Contains("Sum")){
    fit_min = mean - 5*sigma;
    fit_max = mean - 4*sigma;
  }
  else{
    fit_min = mean - 4*sigma;
    fit_max = mean - 3*sigma;
  }
  
  TF1 *fun_pol0 = new TF1("fun_pol0", "pol0", mean+3*sigma, mean+5*sigma);
  fSpectrum->Fit("fun_pol0", opt);
  
  TF1 *fun_expo = new TF1("fun_expo", "pol0(0)+expo(1)", fit_min, fit_max);
  fun_expo->FixParameter(0, fun_pol0->GetParameter(0));
  fSpectrum->Fit("fun_expo", opt);
  
  // background subtraction
  double x, y;
  
  for(int i=0; i<fSpectrum->GetNbinsX(); i++){
    x = fSpectrum->GetBinCenter(i);
    if(x>peak_min && x<peak_max){
      if(fun_expo->Eval(x) > fun_pol0->Eval(x))
        y = fSpectrum->GetBinContent(i) - fun_expo->Eval(x);
      else
        y = fSpectrum->GetBinContent(i) - fun_pol0->Eval(x);
    }
    else 
        y = 0;
    fPeak->SetBinContent(i,y);
  }
  
  TF1 *fun_gaus = new TF1("fun_gaus", "gaus", peak_min, peak_max);
  fun_gaus->SetParameters(fPeak->GetBinContent(fPeak->GetMaximumBin()),
                          fPeak->GetBinCenter(fPeak->GetMaximumBin()),
                          30);
  TFitResultPtr fitRes = fPeak->Fit("fun_gaus", opt);
  
  fun_gaus = fPeak->GetFunction("fun_gaus");
  
  fPosition = fun_gaus->GetParameter(1);
  fPosErr   = fun_gaus->GetParError(1);
  fSigma    = fun_gaus->GetParameter(2);
  fSigErr   = fun_gaus->GetParError(2);
  fChi2NDF  = fitRes->Chi2()/fitRes->Ndf();
  
  // for tests
  if(fTests){
    TFile *file = new TFile("PeakFinderTests.root", "UPDATE");
    fSpectrum->Write();
    fPeak->Write();
    file->Close();
  }
  
  return true;
}
//------------------------------------------------------------------
/// Returns vector containing peak position and sigma along with their errors. 
/// Order in the vector: position, sigma, position error, sigma error.
std::vector <double> SFPeakFinder::GetParameters(){
  std::vector <double> temp;
  temp.push_back(fPosition);
  temp.push_back(fSigma);
  temp.push_back(fPosErr);
  temp.push_back(fSigErr);
  return temp;
}
//------------------------------------------------------------------
/// Prints details of SFPeakFinder class object.
void SFPeakFinder::Print(void){
 std::cout << "\n------------------------------------------------" << std::endl;
 std::cout << "This is print out of SFPeakFinder class object" << std::endl;
 std::cout << "Analyzed spectrum: " << std::endl;
 if(fSpectrum!=nullptr) fSpectrum->Print();
 else std::cout << "NULL" << std::endl;
}
//------------------------------------------------------------------
