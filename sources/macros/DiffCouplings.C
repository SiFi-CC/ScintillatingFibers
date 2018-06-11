R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFData.hh" 

bool DiffCouplings(void){
 
  SFData *data = new SFData(6);
  data->Print();
  int npoints = data->GetNpoints();
  
  vector <TH1D*> hCh0 = data->GetSpectra(0,"fPE","ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0");
  vector <TH1D*> hCh1 = data->GetSpectra(1,"fPE","ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0");
  
  double maxCh0, maxCh1;
  
  for(int i=0; i<npoints; i++){
   hCh0[i]->Scale(1./hCh0[i]->GetEntries());
   hCh1[i]->Scale(1./hCh1[i]->GetEntries());
   hCh0[i]->SetLineColor(i+1);
   hCh1[i]->SetLineColor(i+1);
   hCh0[i]->SetStats(0);
   hCh1[i]->SetStats(0);
   if(i==0){
     hCh0[i]->GetXaxis()->SetRangeUser(0,800);
     hCh1[i]->GetXaxis()->SetRangeUser(0,800);
     hCh0[i]->GetYaxis()->SetMaxDigits(2);
     hCh1[i]->GetYaxis()->SetMaxDigits(2);
     hCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
     hCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
     hCh0[i]->GetYaxis()->SetTitle("counts normalized");
     hCh1[i]->GetYaxis()->SetTitle("counts normalized");
     hCh0[i]->SetTitle("Influence of coupling - Ch0");
     hCh1[i]->SetTitle("Influence of coupling - Ch1");
     maxCh0 = hCh0[i]->GetBinContent(hCh0[i]->GetMaximumBin());
     maxCh1 = hCh1[i]->GetBinContent(hCh1[i]->GetMaximumBin());
   }
   if(hCh0[i]->GetBinContent(hCh0[i]->GetMaximumBin())>maxCh0)
     maxCh0 = hCh0[i]->GetBinContent(hCh0[i]->GetMaximumBin());
   if(hCh1[i]->GetBinContent(hCh1[i]->GetMaximumBin())>maxCh1)
     maxCh1 = hCh1[i]->GetBinContent(hCh1[i]->GetMaximumBin());
  }
  
  hCh0[0]->GetYaxis()->SetRangeUser(0,maxCh0+0.1*maxCh0);
  hCh1[0]->GetYaxis()->SetRangeUser(0,maxCh1+0.1*maxCh1);
  
  TCanvas *can = new TCanvas("diff_couplings","diff_couplings",1500,600);
  can->Divide(2,1);
  
  for(int i=0; i<npoints; i++){
    if(i==0){
     can->cd(1);
     gPad->SetGrid(1,1);
     hCh0[i]->DrawClone();
     can->cd(2);
     gPad->SetGrid(1,1);
     hCh1[i]->DrawClone();
   }
   else{
     can->cd(1);
     hCh0[i]->DrawClone("same");
     can->cd(2);
     hCh1[i]->DrawClone("same");  
    }
  }
  
  can->cd(1);
  TPad *padCh0 = new TPad("padCh0","padCh0",0.38,0.38,0.88,0.88);
  padCh0->Draw();
  gPad->SetGrid(1,1);
  
  can->cd(2);
  TPad *padCh1 = new TPad("padCh1","padCh1",0.38,0.38,0.88,0.88);
  padCh1->Draw();
  
  for(int i=0; i<npoints; i++){
    if(i==0){
     hCh0[i]->GetXaxis()->SetRangeUser(120,350);
     hCh1[i]->GetXaxis()->SetRangeUser(120,350);
     hCh0[i]->GetYaxis()->SetRangeUser(0,6E-3);
     hCh1[i]->GetYaxis()->SetRangeUser(0,6E-3);
     hCh0[i]->GetXaxis()->SetTitle("");
     hCh1[i]->GetXaxis()->SetTitle("");
     hCh0[i]->GetYaxis()->SetTitle("");
     hCh1[i]->GetYaxis()->SetTitle("");
     hCh0[i]->SetTitle("");
     hCh1[i]->SetTitle("");
     padCh0->cd();
     gPad->SetGrid(1,1);
     hCh0[i]->Draw();
     padCh1->cd();
     gPad->SetGrid(1,1);
     hCh1[i]->Draw();
    }
    padCh0->cd();
    hCh0[i]->Draw("same");
    padCh1->cd();
    hCh1[i]->Draw("same");
  }
  
  
  return true;
}
