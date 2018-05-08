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
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include "SFTimingRes.hh"
#include "SFTimeConst.hh"
#include "TCanvas.h"
#include "TF1.h"

int main(){
  
  int seriesNo = 3;
  SFData *data = new SFData(seriesNo);
  data->Print();
  int n = data->GetNpoints();
  double *positions = data->GetPositions();
  
  //----- tests of SFData class
  
  //TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fT0>0",10);
  //TH1D *s1 = data->GetSignal(0,60,"",10,true);
  //TH1D *s2 = data->GetSignal(0,70,"",14,false);
  //TH1D *s3 = data->GetSignal(0,80,"fAmp<100",1,true);
  
  //vector <TH1D*> hh1 = data->GetSpectra(0,"fAmp","");
  //vector <TH1D*> hh2 = data->GetSpectra(1,"fT0","ch_1.fT0!=-100");
  //vector <TH1D*> rr1 = data->GetRatios("log(ch_0.fPE/ch_1.fPE)","ch_0.fT0<590 && ch_0.fPE>0 && ch_1.fT0<590 && ch_1.fPE>0");
  
  TProfile *p1 = data->GetSignalAverage(0,50,"ch_0.fPE>99.5 && ch_0.fPE<100.5",20,true);
  TProfile *p2 = data->GetSignalAverage(1,50,"ch_0.fPE>99.5 && ch_0.fPE<100.5",20,true);
  TH1D *sig1   = data->GetSignal(0,50,"fAmp>200",1,true);
  TH1D *sig2   = data->GetSignal(1,50,"fAmp>50 && fAmp<70",1,true);
  TH1D *hist1  = data->GetSpectrum(0,"fPE","ch_0.fT0>0 && ch_0.fT0<590 && ch_0.fPE>0",10);
  
  TCanvas *can_sig = new TCanvas("can_sig","can_sig",1200,1200);
  can_sig->Divide(2,2);
  
  can_sig->cd(1);
  gPad->SetGrid(1,1);
  p1->GetXaxis()->SetTitle("time [ns]");
  p1->GetYaxis()->SetTitle("amplitude [mV]");
  p1->Draw();
  can_sig->cd(2);
  gPad->SetGrid(1,1);
  p2->GetXaxis()->SetTitle("time [ns]");
  p2->GetYaxis()->SetTitle("amplitude [mV]");
  p2->Draw();
  can_sig->cd(3);
  gPad->SetGrid(1,1);
  sig1->GetXaxis()->SetTitle("time [ns]");
  sig1->GetYaxis()->SetTitle("amplitude [mV]");
  sig1->Draw();
  can_sig->cd(4);
  gPad->SetGrid(1,1);
  sig2->GetXaxis()->SetTitle("time [ns]");
  sig2->GetYaxis()->SetTitle("amplitude [mV]");
  sig2->Draw();
  
  TCanvas *can_spectrum = new TCanvas("can_spectrum","can_spectrum",1500,500);
  can_spectrum->Divide(2,1);
  can_spectrum->cd(1);
  gPad->SetGrid(1,1);
  hist1->GetXaxis()->SetRangeUser(0,1000);
  hist1->GetXaxis()->SetTitle("PE");
  hist1->GetYaxis()->SetTitle("counts");
  hist1->Draw();
  
  can_spectrum->cd(2);
  gPad->SetLogy();
  gPad->SetGrid();
  hist1->Draw();
 
  TFile *fileSFData = new TFile("../results/SFData_test.root","RECREATE");
  fileSFData->cd();
  p1->Write();
  p2->Write();
  sig1->Write();
  sig2->Write();
  hist1->Write();
  can_spectrum->Write();
  can_sig->Write();
  /*
  h1->Write();
  s1->Write();
  s2->Write();
  s3->Write();
  p1->Write();
  p2->Write();
  for(int i=0; i<n; i++){
    hh1[i]->Write();
    hh2[i]->Write();
    rr1[i]->Write();
  }
  */
  fileSFData->Close();
 
  //-----
  
  
  //----- this part is to see if signals with low 
  //----- cut on PE are negative after BL subtraction
  /*
  TCanvas *can = new TCanvas("can","can",1000,1000);
  can->Divide(3,3);
  TH1D *sig[10];
  
  for(int i=0; i<n; i++){
    sig[i] = data->GetSignal(0,50,"ch_0.fPE>19.99 && ch_0.fPE<20.01",i+1,true);
    can->cd(i+1);
    sig[i]->Draw();
  }
  
  TFile *fileSignals = new TFile("../results/Signals_test.root","RECREATE");
  fileSignals->cd();
  
  for(int i=0; i<n; i++){
   sig[i]->Write();
  }
  fileSignals->Close();
  */
  //-----
  
  
  //----- this part is for tests of SFPeakFinder class
  /*
  std::vector <TH1D*> h1;
  std::vector <TH1D*> h2;
  std::vector <SFPeakFinder*> peakfin;
  
  
  for(int i=0; i<n; i++){
	h1.push_back(data->GetSpectrum(1,"fPE","ch_1.fT0>0 && ch_1.fT0<590 && ch_1.fPE>0",positions[i]));
	peakfin.push_back(new SFPeakFinder(h1[i],"511",false));
	peakfin[i]->Print();
	h2.push_back(peakfin[i]->GetPeak());
  }
  
  TCanvas *can = new TCanvas("can","can",1200,1000);
  can->Divide(3,3);
  
  for(int i=0; i<n; i++){
        can->cd(i+1);
        h1[i]->Draw();
	h2[i]->Draw("same");
  }
  
  TFile *fileSFPeakFinder = new TFile("../results/SFPeakFinder_tests.root","RECREATE");
  fileSFPeakFinder->cd();
  can->Write();
  fileSFPeakFinder->Close();
  */
  //-----
  
  
  //----- this part is for tests of SFAttenuation class
  /*
  SFAttenuation *att;
  try{
    att = new SFAttenuation(5);
  } 
  catch(const char* message){
    cout << message << endl;
    return 0;
  }
  
  att->Print();
  
  //att->AttAveragedCh();
  //vector <TH1D*> attRatios = att->GetRatios();
  //TGraphErrors *attGraph = att->GetAttGraph();
  //vector <double> attValue = att->GetAttenuation();
  //cout << "Attenuation is: " << attValue[0] << "\t" << attValue[1] << endl;
 
  int ch = 1;
  att->AttSeparateCh(ch);
  vector <TH1D*> spectra = att->GetSpectra(ch);
  vector <TH1D*> peaks = att->GetPeaks(ch);
  TGraphErrors *attGraphCh = att->GetAttGraph(ch);
  vector <double> attValueCh = att->GetAttenuation(ch);
  cout << "Attenuation for ch " << ch << " is " << attValueCh[0] << "\t" << attValueCh[1] << endl;
  
  TFile *fileSFAttenuation = new TFile("../results/SFAttenuation_tests.root","RECREATE");
  fileSFAttenuation->cd();
  //attGraph->Write();
  attGraphCh->Write();
  for(int i=0; i<9; i++){
    //attRatios[i]->Write();
    spectra[i]->Write();
    peaks[i]->Write();
  }
  fileSFAttenuation->Close();
  */
  //-----

  
  //----- this part is for tests of SFTimingRes class
  /*
  SFTimingRes *tim;
  try{
    tim = new SFTimingRes(1,"ft","no cut");
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  tim->Print();
  vector <TH1D*> T0Diff = tim->GetT0Diff();
  vector <double> tres = tim->GetTimingResolutions();
  vector <double> treserr = tim->GetTimingResErrors();
  TGraphErrors *graph = tim->GetT0Graph();
  
  SFTimingRes *tim_cut;
  try{
    tim_cut = new SFTimingRes(1,"ft","with cut");
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  tim_cut->Print();
  vector <TH1D*> T0Diff_cut = tim_cut->GetT0Diff();
  vector <double> tres_cut = tim_cut->GetTimingResolutions();
  vector <double> treserr_cut = tim_cut->GetTimingResErrors();
  TGraphErrors *graph_cut = tim_cut->GetT0Graph();
  
  TFile *fileSFTiming = new TFile("../results/SFTimingRes_tests.root","RECREATE");
  fileSFTiming->cd();
  for(int i=0; i<n; i++){
    cout << "Energy resolution for position: " << positions[i] << "\t" << tres[i] << " +/- " << treserr[i] << endl;
    cout << "With energy cut: " << tres_cut[i] << " +/- " << treserr_cut[i] << endl;
    T0Diff[i]->Write();
    T0Diff_cut[i]->Write();
  }
  graph->Write();
  graph_cut->Write();
  fileSFTiming->Close();
  */
  //-----
  
  
  //----- this part is to test SFTimeConst class
  /*
  SFTimeConst *tconst1;
  try{
    tconst1 = new SFTimeConst(seriesNo,60,false);
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  tconst1->Print();
  tconst1->FitAllSignals("double");
  vector <TProfile*> signals1 = tconst1->GetAllSignals();
  vector <SFFitResults*> results1 = tconst1->GetAllResults();
 
  
  TCanvas *can1 = new TCanvas("can","can",1000,400);
  can1->Divide(5,2);
  
  for(int i=0; i<9; i++){
   can1->cd(i+1);
   signals1[i]->Draw();
   results1[i]->GetDecayFunction()->Draw("same");
  }

  TFile *fileSFTimeConst = new TFile("../results/SFTimeConst_tests.root","RECREATE");
  fileSFTimeConst->cd();

  for(int i=0; i<n; i++){
    signals1[i]->Write(); 

  }
  can1->Write();
  
  fileSFTimeConst->Close();
  */
  //----- 
  
  //delete data;
  //delete peakfin;
  //delete att;
  //delete tim;
  //delete tim_cut;
  //delete tconst1;
 
  return 1;
}