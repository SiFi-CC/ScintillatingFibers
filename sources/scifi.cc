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
#include "TTimeStamp.h"
#include "TSystem.h"
using namespace std;

int main(){
 
  TTimeStamp starttime;
  gSystem->Exec("date");
  
  SFData *data = new SFData(1);
  data->Print();
  
  TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fPE>5",30.0);
  TH1D *h2 = data->GetSpectrum(0,"fCharge","",30.0);
  
  int n = data->GetNpoints();
  TH1D **spectra = data->GetSpectra(0,"fAmp","");
  
  TProfile *prof = data->GetSignalAverage(0,50.0,"",20,true);
  TProfile *prof2 = data->GetSignalAverage(0,50.0,"",20,false);
  TProfile *prof3 = data->GetSignalAverage(0,50.0,"ch_0.fPE>99.99 && ch_0.fPE100.01",50,true);
  
  TH1D *sig = data->GetSignal(0,50.0,"",10,true);
  TH1D *sig2 = data->GetSignal(0,50.0,"",14,false);
  TH1D *sig3 = data->GetSignal(0,50,"fAmp>100 && fAmp<200",1,true);
  
  
  
  TFile *f = new TFile("test.root","RECREATE");
  
  for(int i=0; i<n; i++){
    spectra[i]->Write();
  }
  
  prof->Write();
  prof2->Write();
  prof3->Write();
  sig->Write();
  sig2->Write();
  sig3->Write();
  h1->Write();
  h2->Write();
 
  
  TTimeStamp stoptime;
  
  cout << "this took " << stoptime.GetSec()-starttime.GetSec() << " s" << endl;
  
  SFData *dd = new SFData(9);
  
  TH1D **darkc = dd->GetSpectra(0,"fCharge","");
  
  for(int i=0; i<n; i++){
    darkc[i]->Write(); 
  }
  
  f->Close();
 
  delete data;
  delete dd;
  
  return 1;
}