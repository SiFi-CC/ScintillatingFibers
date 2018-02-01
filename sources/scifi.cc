// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               scifi.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
using namespace std;

int main(){
 
  SFData *data = new SFData(1);
  data->Print();
  
  TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fPE>5",30.0);
  TH1D *h2 = data->GetSpectrum(0,"fCharge","",30.0);
  
  int n = data->GetNpoints();
  TH1D **spectra = data->GetSpectra(0,"fAmp","");
  
  TProfile *prof = data->GetSignalAverage(0,50.0,20,true);
  TProfile *prof2 = data->GetSignalAverage(0,50.0,20,false);
  
  TH1D *sig = data->GetSignal(0,50.0,10,true);
  TH1D *sig2 = data->GetSignal(0,50.0,14,false);
  
  TFile *f = new TFile("test.root","RECREATE");
  
  for(int i=0; i<n; i++){
    spectra[i]->Write();
  }
  
  prof->Write();
  prof2->Write();
  sig->Write();
  sig2->Write();
  h1->Write();
  h2->Write();
  f->Close();
  
  delete data;
  
  return 1;
}