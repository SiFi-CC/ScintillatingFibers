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
  
  TFile *f = new TFile("test.root","RECREATE");
  
  for(int i=0; i<n; i++){
    spectra[i]->Write();
  }
  
  h1->Write();
  h2->Write();
  f->Close();
  
  delete data;
  
  return 1;
}