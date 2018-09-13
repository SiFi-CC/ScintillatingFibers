// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            attenuation.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFAttenuation.hh"
#include "TCanvas.h"
#include "TLatex.h"

int main(int argc, char **argv){
  
  if(argc!=2){
    cout << "to run type: ./attenuation seriesNo" << endl;
    return 0;
  }
 
  int seriesNo = atoi(argv[1]);

  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char* message){
   cout << message << endl;
   cout << "##### Exception in attenuation.cc!" << endl;
   return 0;
  }
 
  data->Print();
 
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    cout << "##### Error in attenuation.cc! This is not regular series!" << endl;
    cout << "Series number: " << seriesNo << endl;
    cout << "Description: " << desc << endl;
    return 0;
  }
 
  TString mtype= data->GetMeasureType();
 
  int npoints = data->GetNpoints();
  vector <double> positions = data->GetPositions();
  data->Print();
 
  SFAttenuation *att;
  try{
    att = new SFAttenuation(seriesNo);
  }
  catch(const char* message){
    cout << message << endl;
    cout << "##### Exception in attenuation.cc!" << endl;
    return 0;
  }
  
  //----- averaged channels method
  att->AttAveragedCh();
  vector <TH1D*> attRatios = att->GetRatios();
  TGraphErrors *attGraph   = att->GetAttGraph();
  vector <double> attlen   = att->GetAttenuation();
  
  //----- separate channels method
  att->AttSeparateCh(0);
  TGraphErrors *attGraphCh0  = att->GetAttGraph(0);
  vector <TH1D*> spectraCh0  = att->GetSpectra(0);
  vector <TH1D*> peaksCh0    = att->GetPeaks(0);
  vector <double> attlenCh0  = att->GetAttenuation(0);
  
  att->AttSeparateCh(1);
  TGraphErrors *attGraphCh1 = att->GetAttGraph(1);
  vector <TH1D*> spectraCh1 = att->GetSpectra(1);
  vector <TH1D*> peaksCh1   = att->GetPeaks(1);
  vector <double> attlenCh1 = att->GetAttenuation(1);    
  
  //-----drawing averaged channels
  TLatex text;
  text.SetNDC(true);
  
  TCanvas *can_averaged_ch = new TCanvas("can_averaged_ch","can_averaged_ch",700,500);
  can_averaged_ch->cd();
  gPad->SetGrid(1,1);
  attGraph->SetTitle(Form("Series %i, attenuation curve",seriesNo));
  attGraph->Draw();
  text.SetTextSize(0.04);
  text.DrawLatex(0.2,0.8,Form("L_{att} = (%.2f +/- %.2f) mm",attlen[0],attlen[1]));
  
  TCanvas *can_ratios = new TCanvas("can_ratios","can_ratios",1200,1200);
  can_ratios->Divide(3,3);
  TString htitle;
  text.SetTextSize(0.02);
  
  for(int i=0; i<npoints; i++){
   can_ratios->cd(i+1);
   gPad->SetGrid(1,1);
   attRatios[i]->SetTitle(Form("ln(#sqrt{ch1/ch0}), source position %.2f mm",positions[i]));
   attRatios[i]->GetXaxis()->SetRangeUser(-2,2);
   attRatios[i]->GetXaxis()->SetTitle("ln(#sqrt{ch1/ch0})");
   attRatios[i]->Draw();
  }
  
  //----- drawing separate channels 
  text.SetTextSize(0.04);
  TCanvas *can_separate_ch = new TCanvas("can_separate_ch","can_separate_ch",1000,500);
  can_separate_ch->Divide(2,1);
  can_separate_ch->cd(1);
  gPad->SetGrid(1,1);
  attGraphCh0->SetTitle(Form("Series %i channel 0, attenuation curve",seriesNo));
  attGraphCh0->GetYaxis()->SetTitleSize(0.03);
  attGraphCh0->Draw();
  text.DrawLatex(0.3,0.8,Form("L_{att} = (%.2f +/- %.2f) mm",attlenCh0[0],attlenCh0[1]));
  can_separate_ch->cd(2);
  gPad->SetGrid(1,1);
  attGraphCh1->SetTitle(Form("Series %i channel 1, attenuation curve",seriesNo));
  attGraphCh1->GetYaxis()->SetTitleSize(0.03);
  attGraphCh1->Draw();
  text.DrawLatex(0.3,0.8,Form("L_{att} = (%.2f +/- %.2f) mm",attlenCh1[0],attlenCh1[1]));
  
  TCanvas *can_spectra_ch0 = new TCanvas("can_spectra_ch0","can_spectra_ch0",1200,1200);
  if(mtype.Contains("Lead")) can_spectra_ch0->Divide(3,3);
  else if(mtype.Contains("Electric")) can_spectra_ch0->Divide(2,3);
  
  TCanvas *can_spectra_ch1 = new TCanvas("can_spectra_ch1","can_spectra_ch1",1200,1200);
  if(mtype.Contains("Lead")) can_spectra_ch1->Divide(3,3);
  else if(mtype.Contains("Electric")) can_spectra_ch1->Divide(2,3);
  
  for(int i=0; i<npoints; i++){
   can_spectra_ch0->cd(i+1);
   gPad->SetGrid(1,1);
   if(mtype.Contains("Lead")){
	spectraCh0[i]->SetStats(false);
	spectraCh0[i]->GetXaxis()->SetRangeUser(0,800);
	spectraCh0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm",seriesNo,positions[i]));
	spectraCh0[i]->GetXaxis()->SetTitle("P.E.");
	spectraCh0[i]->GetYaxis()->SetTitle("counts");
	spectraCh0[i]->GetYaxis()->SetMaxDigits(2);
	spectraCh0[i]->Draw();
	peaksCh0[i]->Draw("same");
   }
   else if(mtype.Contains("Electric")){
	peaksCh0[i]->SetStats(false);
	peaksCh0[i]->GetXaxis()->SetRangeUser(0,800);
	peaksCh0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm",seriesNo,positions[i]));
	peaksCh0[i]->GetXaxis()->SetTitle("P.E.");
	peaksCh0[i]->GetYaxis()->SetTitle("counts");
	peaksCh0[i]->GetYaxis()->SetMaxDigits(2);
   	peaksCh0[i]->Draw();
   }
   can_spectra_ch1->cd(i+1);
   gPad->SetGrid(1,1);
   if(mtype.Contains("Lead")){
	spectraCh1[i]->SetStats(false);
	spectraCh1[i]->GetXaxis()->SetRangeUser(0,800);
	spectraCh1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm",seriesNo,positions[i]));
	spectraCh1[i]->GetXaxis()->SetTitle("P.E.");
	spectraCh1[i]->GetYaxis()->SetTitle("counts");
	spectraCh1[i]->GetYaxis()->SetMaxDigits(2);
	spectraCh1[i]->Draw();
	peaksCh1[i]->Draw("same");
   }
   else if(mtype.Contains("Electric")){
	peaksCh1[i]->SetStats(false);
	peaksCh1[i]->GetXaxis()->SetRangeUser(0,800);
	peaksCh1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm",seriesNo,positions[i]));
	peaksCh1[i]->GetXaxis()->SetTitle("P.E.");
	peaksCh1[i]->GetYaxis()->SetTitle("counts");
	peaksCh1[i]->GetYaxis()->SetMaxDigits(2);
   	peaksCh1[i]->Draw();
   }
  }
  
  //----- saving
  TString fname = Form("../results/attenuation_series%i.root",seriesNo);
  TFile *file = new TFile(fname,"RECREATE");
  can_averaged_ch->Write();
  can_separate_ch->Write();
  can_ratios->Write();
  can_spectra_ch0->Write();
  can_spectra_ch1->Write();
  file->Close();
  
  delete data;
  delete att;
  
  return 1;
}
