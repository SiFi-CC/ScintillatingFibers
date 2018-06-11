R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFData.hh"

bool InternalAct(void){
 
  //----- Internal activity of fiber #1
  
  SFData *data_S7 = new SFData(7);
  data_S7->Print();
  
  TH1D *hFib1_ch0 = data_S7->GetSpectrum(0,"fPE","ch_0.fPE>0",1);
  TH1D *hFib1_ch1 = data_S7->GetSpectrum(1,"fPE","ch_1.fPE>0",1);
  TH1D *hFib1_av  = data_S7->GetCustomHistogram("sqrt(ch_0.fPE*ch_1.fPE)","ch_0.fPE>0 && ch_1.fPE>0",1);
  
  //----- Internal activity of fiber #2
  
  SFData *data_S8 = new SFData(8);
  data_S8->Print();
  
  TH1D *hFib2_ch0 = data_S8->GetSpectrum(0,"fPE","ch_0.fPE>0",1);
  TH1D *hFib2_ch1 = data_S8->GetSpectrum(1,"fPE","ch_1.fPE>0",1);
  TH1D *hFib2_av  = data_S8->GetCustomHistogram("sqrt(ch_0.fPE*ch_1.fPE)","ch_0.fPE>0 && ch_1.fPE>0",1);
  
  //----- Fiber #1 - drawing
  
  TCanvas *can_S7 = new TCanvas("int_act_fib1","int_act_fib1",1200,600);
  can_S7->Divide(2,1);
  
  can_S7->cd(1);
  gPad->SetGrid(1,1);
  hFib1_ch0->SetTitle("Internal activity spectrum - Fiber 1, Ch0");
  hFib1_ch0->GetXaxis()->SetTitle("charge [P.E.]");
  hFib1_ch0->GetYaxis()->SetTitle("counts");
  hFib1_ch0->GetXaxis()->SetRangeUser(0,1000);
  hFib1_ch0->SetStats(0);
  hFib1_ch0->Draw();
  
  can_S7->cd(2);
  gPad->SetGrid(1,1);
  hFib1_ch1->SetTitle("Internal activity spectrum - Fiber 1, Ch1");
  hFib1_ch1->GetXaxis()->SetTitle("charge [P.E.]");
  hFib1_ch1->GetYaxis()->SetTitle("counts");
  hFib1_ch1->GetXaxis()->SetRangeUser(0,1000);
  hFib1_ch1->SetStats(0);
  hFib1_ch1->Draw();
  
  //----- Fiber #2 - drawing
  
  TCanvas *can_S8 = new TCanvas("int_act_fib2","int_act_fib2",1200,600);
  can_S8->Divide(2,1);
  
  can_S8->cd(1);
  gPad->SetGrid(1,1);
  hFib2_ch0->SetTitle("Internal activity spectrum - Fiber 2, Ch0");
  hFib2_ch0->GetXaxis()->SetTitle("charge [P.E.]");
  hFib2_ch0->GetYaxis()->SetTitle("counts");
  hFib2_ch0->GetXaxis()->SetRangeUser(0,1000);
  hFib2_ch0->SetStats(0);
  hFib2_ch0->Draw();
  
  can_S8->cd(2);
  gPad->SetGrid(1,1);
  hFib2_ch1->SetTitle("Internal activity spectrum - Fiber 2, Ch1");
  hFib2_ch1->GetXaxis()->SetTitle("charge [P.E.]");
  hFib2_ch1->GetYaxis()->SetTitle("counts");
  hFib2_ch1->GetXaxis()->SetRangeUser(0,1000);
  hFib2_ch1->SetStats(0);
  hFib2_ch1->Draw();
  
  //----- Fiber #1 and #2 - comparison
  
  TCanvas *can = new TCanvas("can","can",800,600);
  can->cd();
  gPad->SetGrid(1,1);
  
  hFib1_av->SetLineColor(kRed);
  hFib1_av->GetXaxis()->SetRangeUser(0,1000);
  hFib1_av->GetXaxis()->SetTitle("sqrt(ch0*ch1) [P.E.]");
  hFib1_av->GetYaxis()->SetTitle("counts");
  hFib1_av->SetTitle("Fiber 1 and 2 - internal actvity comparison");
  hFib1_av->SetStats(0);
  hFib1_av->Draw();
  hFib2_av->SetStats(0);
  hFib2_av->Draw("same");
  
  TLegend *leg = new TLegend(0.6,0.75,0.78,0.87);
  leg->AddEntry(hFib1_av,"Fiber 1","PEL");
  leg->AddEntry(hFib2_av,"Fiber 2","PEL");
  leg->Draw();
  
  return true;
}