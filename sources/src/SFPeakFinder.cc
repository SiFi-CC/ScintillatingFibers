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
SFPeakFinder::SFPeakFinder(TH1D *spectrum, TString peakID){
  SetSpectrum(spectrum,peakID);
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
  
  return true;
}
//------------------------------------------------------------------
bool SFPeakFinder::FindPeakRange(double &min, double &max){
 
  min = -1;
  max = -1;
  double search_min,search_max;
  
  if(fPeakID=="511"){
    search_min = 220;
    search_max = 400;
  }
  else if(fPeakID=="1270"){
    min = -1;
    max = -1;
  }
  
  int bin_min = fSpectrum->FindBin(search_min);
  int bin_max = fSpectrum->FindBin(search_max);
  fSpectrum->GetXaxis()->SetRange(bin_min,bin_max);
  double peak = fSpectrum->GetBinCenter(fSpectrum->GetMaximumBin());
  
  min = peak-50;
  max = peak+50;
  cout << "max peak: " << peak << endl;
  
  if(max<min || fabs(min+1)<1E-8 || fabs(max+1)<1E-8){
   cout << "##### Error in SFPeakFinder::FindFitRange(). Incorrect range." << endl;
   cout << "min = " << min << "\t max = " << max << endl;
   return false;
  }
  
  
  return true;
}
//------------------------------------------------------------------
bool SFPeakFinder::Fit(void){
  
  //setting background-subtracted histogram
  TString tmp = fSpectrum->GetName();
  TString pname = tmp.Append("_peak");
  fPeak = (TH1D*) fSpectrum->Clone(pname);
  fPeak->Reset();
  
  //fitting background function
  double peak_min, peak_max;
  FindPeakRange(peak_min,peak_max);
  double fit_min = peak_min-20;
  double fit_max = peak_max+50;
  BGFit *bg = new BGFit(peak_min,peak_max);
  TF1 *bg_fun = new TF1("bg_fun",bg,&BGFit::Evaluate,fit_min,fit_max,4,"BGFit","Evaluate");
  fSpectrum->Fit("bg_fun","R");
  
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
  fPeak->Fit(gaus_fun,"R");
  
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
///Sets members of the class to their default values.
void SFPeakFinder::Clear(void){
 fSpectrum = NULL;
 fPeak     = NULL;
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