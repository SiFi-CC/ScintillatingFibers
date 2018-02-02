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
 
  int n =0;
  
  SFData *data = new SFData(1);
  data->Print();
  
  n = data->GetNpoints();
  
  TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fT0>0",10);
  TH1D *h2 = data->GetSpectrum(0,"fCharge","",20);
  
  TFile *f = new TFile("test.root","RECREATE");
  
  TH1D** hh1 = data->GetSpectra(0,"fAmp","");
  TH1D** hh2 = data->GetSpectra(1,"fT0","ch_1.fT0!=-100");
  
  f->cd();
  for(int i=0; i<n; i++){
    hh1[i]->Write();
  }
  
  for(int i=0; i<n; i++){
    hh2[i]->Write();
  }
  
  TProfile *p1 = data->GetSignalAverage(1,30,"",20,true);
  TProfile *p2 = data->GetSignalAverage(1,40,"",20,false);
  TProfile *p3 = data->GetSignalAverage(1,50,"ch_0.fPE>59.99 && ch_0.fPE60.01",50,true);
  
  TH1D *s1 = data->GetSignal(0,60,"",10,true);
  TH1D *s2 = data->GetSignal(0,70,"",14,false);
  TH1D *s3 = data->GetSignal(0,80,"fAmp>100 && fAmp<200",1,true);
  
  f->cd();
  h1->Write();
  h2->Write();
  p1->Write();
  p2->Write();
  p3->Write();
  s1->Write();
  s2->Write();
  s3->Write();
  
  f->Close();

  delete data;

  
  return 1;
}