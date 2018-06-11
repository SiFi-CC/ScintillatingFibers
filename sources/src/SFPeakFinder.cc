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
///Default constructor.
SFPeakFinder::SFPeakFinder(){
  cout << "#### Warning in SFPeakFinder constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  Clear();
}
//------------------------------------------------------------------
///Standard constructor.
///\param spectrum - analyzed spectrum
///\param peakID - flag to identify the peak for analysis. Possible options are: 511 or 1270.
///\param verbose - print-outs level
SFPeakFinder::SFPeakFinder(TH1D *spectrum, TString peakID, bool verbose){
  SetSpectrum(spectrum,peakID);
  fVerbose = verbose;
}
//------------------------------------------------------------------
///Standard constructor. Vebose level set by default to quiet. In order to change verbose
///level use SetVerbLevel().
///\param spectrum - analyzed spectrum
///\param peakID - flag to identify the peak for analysis. Possible options are: 511 or 1270.
SFPeakFinder::SFPeakFinder(TH1D *spectrum, TString peakID){
  SetSpectrum(spectrum,peakID);
  cout << "##### Warning in SFPeakFinder constructor. Quiet mode on, no print outs." << endl;
  cout << "##### To set verbose level use SetVerbLevel(bool)" << endl;
}
//------------------------------------------------------------------
///Default destructor.
SFPeakFinder::~SFPeakFinder(){
}
//------------------------------------------------------------------
///Sets analyzed histogram and peak ID. Calls private function Fit().
///\param spectrum - analyzed spectrum
///\param peakID - flag to identify the peak for analysis
bool SFPeakFinder::SetSpectrum(TH1D *spectrum, TString peakID){
  Clear();
  fSpectrum = spectrum;
  fPeakID = peakID;
  
  if(!(peakID=="511" || peakID=="1270")){
    cout << "##### Error in SFPeakFinder::SetSpectrum(). Incorrect peak ID." << endl;
    cout << "Possible options are: 511 or 1270" << endl;
    return false;
  }
  
  Fit();
  fSpectrum->GetXaxis()->UnZoom();
  
  return true;
}
//------------------------------------------------------------------
///Finds ranges of the analyzed peak. Ranges are returned by reference.
bool SFPeakFinder::FindPeakRange(double &min, double &max){
 
  min = -1;
  max = -1;
  double search_min,search_max;
  
  if(fPeakID=="511"){
    search_min = 120;
    search_max = 450;//320
  }
  else if(fPeakID=="1270"){
    search_min = 400;
    search_max = 1000;
  }
  
  int bin_min = fSpectrum->FindBin(search_min);
  int bin_max = fSpectrum->FindBin(search_max);
  int bin_delta = 0;
  int bin_peak = 0;
  int step = 10;
  
  while(bin_delta<step){
   fSpectrum->GetXaxis()->SetRange(bin_min,bin_max);
   bin_peak = fSpectrum->GetMaximumBin();
   bin_delta = fabs(bin_peak-bin_min);
   //cout << "bin_min = " << bin_min << "\t bin_delta = " << bin_delta <<endl;
   bin_min+=step;
  }
  
  double peak = fSpectrum->GetBinCenter(bin_peak);
  
  //setting fitting option based on verbose level
  TString opt;
  if(fVerbose) opt = "0R";
  else opt = "Q0R";
  
  TF1 *fun = new TF1("fun","gaus",peak-30,peak+30);
  fSpectrum->Fit(fun,opt);
  
  min = peak-fun->GetParameter(2);
  max = peak+fun->GetParameter(2);
  //cout << "min peak: " << min << " max peak: " << peak << endl;
  
  if(max<min || fabs(min+1)<1E-8 || fabs(max+1)<1E-8){
   cout << "##### Error in SFPeakFinder::FindFitRange(). Incorrect range." << endl;
   cout << "min = " << min << "\t max = " << max << endl;
   return false;
  }
  
  fSpectrum->GetXaxis()->UnZoom();
  
  return true;
}
//------------------------------------------------------------------
///Fits background to the analyzed spectrum and performs background subtraction. 
///Inside this function peak position and sigma along with their errors are found. 
///Additionally histogram fPeak is created and filled.
bool SFPeakFinder::Fit(void){
  
  //setting background-subtracted histogram
  TString tmp = fSpectrum->GetName();
  TString pname = tmp.Append("_peak");
  fPeak = (TH1D*) fSpectrum->Clone(pname);
  fPeak->Reset();
  
  //setting fitting option based on verbose level
  TString opt;
  if(fVerbose) opt = "R+";
  else opt = "QR+";
  
  //fitting background function
  double peak_min, peak_max;
  FindPeakRange(peak_min,peak_max);
  //cout << "peak min = " << peak_min << "\t peak max = " << peak_max << endl;
  double fit_min = peak_min-10;
  double fit_max = peak_max+70;
  BGFit *bg = new BGFit(peak_min,peak_max);
  TF1 *bg_fun = new TF1("bg_fun",bg,&BGFit::EvaluateExpo,fit_min,fit_max,2,"BGFit","Evaluate");
  fSpectrum->Fit("bg_fun",opt);
  
  //background subtraction
  int nbins = fSpectrum->GetXaxis()->GetNbins();
  double x, y;
  
  for(int i=1; i<nbins+1; i++){
   x = fSpectrum->GetBinCenter(i);
   if(x>fit_min && x<fit_max){
    y = fSpectrum->GetBinContent(i) - bg_fun->Eval(x);
    if(y<0) y=0;
    fPeak->SetBinContent(i,y);
   }
   else
     fPeak->SetBinContent(i,0);
  }
  
  //fitting Gauss to the peak
  TF1 *gaus_fun = new TF1("fun_gaus","gaus",peak_min,peak_max);
  gaus_fun->SetParameters(100,(peak_max-peak_min)/2.,10);
  fPeak->Fit(gaus_fun,opt);
  
  fPosition = gaus_fun->GetParameter(1);
  fPosErr = gaus_fun->GetParError(1);
  fSigma = gaus_fun->GetParameter(2);
  fSigErr = gaus_fun->GetParError(2);
  
  if(fPosition<0 || fSigma<0){
   cout << "##### Error in SFPeakFinder(). Position and Sigma cannot be negative." << endl;
   cout << "position = " << fPosition << "\t fSigma = " << fSigma << endl;
   return false;
  }
  
  return true;
}
//------------------------------------------------------------------
///Returns vector containing peak position and sigma along with their errors. 
///Order in the vector: position, sigma, position error, sigma error.
vector<double> SFPeakFinder::GetParameter(){
	vector<double> temp;
	temp.push_back(fPosition);
	temp.push_back(fSigma);
	temp.push_back(fPosErr);
	temp.push_back(fSigErr);
	
	return temp;
}
//------------------------------------------------------------------
///Sets members of the class to their default values.
void SFPeakFinder::Clear(void){
 fSpectrum = NULL;
 fPeak     = NULL;
 fVerbose  = false;
 fPeakID   = "dummy";
 fPosition = -100;
 fPosErr   = -100;
 fSigma    = -100;
 fSigErr   = -100; 
}
//------------------------------------------------------------------
///Prints details of SFPeakFinder class object.
void SFPeakFinder::Print(void){
 cout << "\n------------------------------------------------" << endl;
 cout << "This is print out of SFPeakFinder class object" << endl;
 cout << "Analyzed spectrum: " << endl;
 if(fSpectrum!=NULL) fSpectrum->Print();
 else cout << "NULL" << endl;
 cout << "Analyzed peak: " << fPeakID << endl << endl;
}
//------------------------------------------------------------------
