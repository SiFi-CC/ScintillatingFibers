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
#include "TCanvas.h"

int main(){
  
  TFile *f = new TFile("test.root","RECREATE");
  
  SFData *data = new SFData(1);
  data->Print();
  
  int n = data->GetNpoints();
  
  TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fT0>0",10);
  TH1D *h2 = data->GetSpectrum(0,"fCharge","",20);
  
  vector <TH1D*> hh1 = data->GetSpectra(0,"fAmp","");
  vector <TH1D*> hh2 = data->GetSpectra(1,"fT0","ch_1.fT0!=-100");
  
  f->cd();
  for(int i=0; i<n; i++){
    hh1[i]->Write();
    hh2[i]->Write();
  }
  
  TProfile *p1 = data->GetSignalAverage(1,30,"",20,true);
  TProfile *p2 = data->GetSignalAverage(1,40,"",20,false);
  TProfile *p3 = data->GetSignalAverage(1,50,"ch_0.fPE>59.99 && ch_0.fPE<60.01",50,true);
  
  TH1D *s1 = data->GetSignal(0,60,"",10,true);
  TH1D *s2 = data->GetSignal(0,70,"",14,false);
  TH1D *s3 = data->GetSignal(0,80,"fAmp<100",1,true);
  
  //----- this part is to see if signals with low 
  //----- cut on PE are negative after bl subtraction
  TCanvas *can = new TCanvas("can","can",1000,1000);
  can->Divide(3,3);
  TH1D *sig[9];
  
  for(int i=0; i<9; i++){
    sig[i] = data->GetSignal(0,50,"ch_0.fPE>19.99 && ch_0.fPE<20.01",i+1,true);
    can->cd(i+1);
    sig[i]->Draw();
  }
  
  f->cd();
  can->Write();
  //-----
  
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