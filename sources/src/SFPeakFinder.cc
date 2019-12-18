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
/// \param tests - flag for testing mode
SFPeakFinder::SFPeakFinder(TH1D *spectrum, bool verbose, bool tests): fSpectrum(spectrum),
                                                                      fPeak(nullptr),
                                                                      fPosition(-1),
                                                                      fPosErr(-1),
                                                                      fSigma(-1),
                                                                      fSigErr(-1),
                                                                      fVerbose(verbose),
                                                                      fTests(tests) { 
                                                              
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
                                            fVerbose(false),
                                            fTests(false) {

  std::cout << "##### Warning in SFPeakFinder constructor. Quiet mode on, no print outs." << std::endl;
  std::cout << "##### To set verbose level use SetVerbLevel(bool)" << std::endl;
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
    std::cerr << "Spectrum cannot be a null pointer!" << std::endl;
    std::abort();
  }
  fSpectrum = spectrum;
  return;
}
//------------------------------------------------------------------
TString SFPeakFinder::Init(void){
    
  TString hname = fSpectrum->GetName();
  int ID = SFTools::GetMeasurementID(hname);  
  int seriesNo = SFTools::GetSeriesNo(hname);
  
  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Error in SFPeakFinder::Init()!" << std::endl;
  }
  
  std::vector <TString> names = data->GetNames();
  std::vector <int> measureID = data->GetMeasurementsIDs();
  std::vector <double> positions = data->GetPositions();
  int index = SFTools::GetIndex(measureID, ID);
  TString dir_name = names[index];
  TString full_path = SFTools::FindData(dir_name);
  
  TString conf_name = "/fitconfig.txt";
  
  TString functions = "gaus(0) pol1(3)+[5]*TMath::Exp((x-[6])*[7])";
  TString hnames[4] = {Form("S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[index], ID),
                       Form("S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[index], ID),
                       Form("S%i_pos%.1f_ID%i_PEAverage", seriesNo, positions[index], ID),
                       Form("S%i_pos%.1f_ID%i_PEAttCorrectedSum", seriesNo, positions[index], ID)};
  
  std::fstream test(full_path+conf_name, std::ios::in);
  
  if(test.fail()){
    std::cout << "Fitting config for " << full_path << " doesn't exist..." << std::endl;
    std::cout << "Creating new config file..." << std::endl;
    
    TSpectrum *spec = new TSpectrum(10);
    int npeaks = spec->Search(fSpectrum, 20, "goff", 0.1);
    double *peaksX = spec->GetPositionX();
    double peak = TMath::MaxElement(npeaks, peaksX);
    
    TString opt;
    if(fVerbose) opt = "0R";
    else opt = "Q0R";  
    
    TF1 *fun_gaus = new TF1("fun_gaus", "gaus", peak-100, peak+100);
    fSpectrum->Fit("fun_gaus", opt);
    
    double par0 = fun_gaus->GetParameter(0);
    double par1 = fun_gaus->GetParameter(1);
    double par2 = fun_gaus->GetParameter(2);
    double par2_min = 0; double par2_max = 300;
    
    TF1 *fun_pol1 = new TF1("fun_pol1", "pol1", par1+3*par2, par1+4*par2);
    fSpectrum->Fit("fun_pol1", opt);
    
    double par3 = fun_pol1->GetParameter(0);
    double par4 = fun_pol1->GetParameter(1);
    
    TF1 *fun_expo = new TF1("fun_expo", "[0]*TMath::Exp((x-[1])*[2])", par1-5*par2, par1-4*par2);
    fSpectrum->Fit("fun_expo", opt);
    
    double par5 = fun_expo->GetParameter(0);
    double par6 = fun_expo->GetParameter(1);
    double par7 = fun_expo->GetParameter(2);
    
    double xmin = par1 - par2*4;
    double xmax = par1 + par2*5;
    
    std::fstream config(full_path+conf_name, std::ios::out);
      for(int i=0; i<4; i++){
        config << " " << hnames[i] << " " << functions << " " 
               << 0 << " " << xmin << " " << xmax << " " 
               << par0 << " " << par1 << " " << par2 << " : " 
               << par2_min << " " << par2_max << " " << par3 
               << " " << par4 << " " << par5 << " " << par6 
               << " " << par7 << "\n";
      }
    config.close();
  }
  else{
    std::cout << "Fitting config for " << full_path << " exists!" << std::endl;
    test.close();
  }
  
  return full_path;  
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
  const double delta = 1E-8;
  
  //if(fabs(fSigma+1)<delta || fabs(fPosition+1)<delta){
    FindPeakFit();
  //}

  if(type=="Lead"){
    min = fPosition - fSigma;
    max = fPosition + 1.5*fSigma;
  }
  else if(type=="Electronic"){
    min = fPosition - 2*fSigma;
    max = fPosition + 2*fSigma;
  }

  if(max<min || fabs(min+1)<delta || fabs(max+1)<delta){
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
  
  TString data_path = Init();

  FitterFactory fitter;
  fitter.initFactoryFromFile((data_path+"/fitconfig.txt").Data(), 
                             (data_path+"/fitparams.out").Data());

  HistFitParams histFP;
  fitter.findParams(fSpectrum->GetName(), histFP);
  fitter.fit(histFP, fSpectrum);
  fitter.updateParams(fSpectrum, histFP);
  fitter.exportFactoryToFile();

  TString fun_name = std::string("f_")+fSpectrum->GetName();
  TF1 *fun_fitted = histFP.funSum;
  
  if(fun_fitted == nullptr){
    std::cerr << "##### Error in SFPeakFinder::FindPeak()! Function is null pointer" << std::endl;
    std::abort();
  }

  fPosition = fun_fitted->GetParameter(1);
  fPosErr = fun_fitted->GetParError(1);
  fSigma = fun_fitted->GetParameter(2);
  fSigErr = fun_fitted->GetParError(2);
    
  // for tests
  if(fTests){
  
    double fit_min = 0;
    double fit_max = 0;
    fun_fitted->GetRange(fit_min, fit_max);
      
    TString opt;
    if(fVerbose) opt = "BS+";
    else opt = "BSQ+";
  
    TF1 *fun_expo_clone = new TF1("fun_expo_clone", "[0]*TMath::Exp((x-[1])*[2])", fit_min, fit_max); 
    fun_expo_clone->SetLineColor(kGreen+3);
    fun_expo_clone->FixParameter(0, fun_fitted->GetParameter(5));
    fun_expo_clone->FixParameter(1, fun_fitted->GetParameter(6));
    fun_expo_clone->FixParameter(2, fun_fitted->GetParameter(7));
    fSpectrum->Fit("fun_expo_clone", opt);
    
    if(fVerbose){
      std::cout << "Fitted functions parameters: \n" 
                << fun_fitted->GetParameter(5) << "\t" 
                << fun_fitted->GetParameter(6) << "\t" 
                << fun_fitted->GetParameter(7) << std::endl;
                
      std::cout << "Expo clone parameters: \n" 
                << fun_expo_clone->GetParameter(0) << "\t" 
                << fun_expo_clone->GetParameter(1) << "\t" 
                << fun_expo_clone->GetParameter(2) << std::endl;
    }
  
    TF1 *fun_pol1_clone = new TF1("fun_pol1_clone", "pol1", fit_min, fit_max);
    fun_pol1_clone->SetLineColor(kBlue);
    fun_pol1_clone->FixParameter(0, fun_fitted->GetParameter(3));
    fun_pol1_clone->FixParameter(1, fun_fitted->GetParameter(4));
    fSpectrum->Fit("fun_pol1_clone", opt);
    
    if(fVerbose){
      std::cout << "Fitted function parameters: \n"
                << fun_fitted->GetParameter(3) << "\t" 
                << fun_fitted->GetParameter(4) << std::endl;
                
      std::cout << "Pol1 clone parameters: \n" 
                << fun_pol1_clone->GetParameter(0) << "\t" 
                << fun_pol1_clone->GetParameter(1) << std::endl;
    }
  
    TF1 *fun_gaus_clone = new TF1("fun_gaus_clone", "gaus", fit_min, fit_max);
    fun_gaus_clone->SetLineColor(kMagenta);
    fun_gaus_clone->FixParameter(0, fun_fitted->GetParameter(0));
    fun_gaus_clone->FixParameter(1, fun_fitted->GetParameter(1));
    fun_gaus_clone->FixParameter(2, fun_fitted->GetParameter(2));
    fSpectrum->Fit("fun_gaus_clone", opt);
    
    if(fVerbose){
      std::cout << "Fitted function parameters: \n" 
                << fun_fitted->GetParameter(0) << "\t" 
                << fun_fitted->GetParameter(1) <<  "\t" 
                << fun_fitted->GetParameter(2) << std::endl;
      
      std::cout << "Gaus clone parameters: \n" 
                << fun_gaus_clone->GetParameter(0) << "\t" 
                << fun_gaus_clone->GetParameter(1) <<  "\t" 
                << fun_gaus_clone->GetParameter(2) << std::endl;
    }
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
/*
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
  
  // for tests
  if(fTests){
    TFile *file = new TFile("PeakFinderTests.root", "UPDATE");
    fSpectrum->Write();
    fPeak->Write();
    file->Close();
  }
  
  return true;
}*/
//------------------------------------------------------------------
/// Returns vector containing peak position and sigma along with their errors. 
/// Order in the vector: position, sigma, position error, sigma error.
std::vector <double> SFPeakFinder::GetParameters(void){
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
