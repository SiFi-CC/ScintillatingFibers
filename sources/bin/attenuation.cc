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
#include <sys/stat.h> 
#include <sys/types.h> 
#include "common_options.h"

int main(int argc, char **argv){
  
  if(argc<2 || argc>6){
    std::cout << "to run type: ./attenuation seriesNo";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }
 
  int seriesNo = atoi(argv[1]);

  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char* message){
   std::cerr << message << std::endl;
   std::cerr << "##### Exception in attenuation.cc!" << std::endl;
   return 1;
  }
 
  data->Print();
 
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in attenuation.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cout << "Description: " << desc << std::endl;
    return 1;
  }
 
  int npoints = data->GetNpoints();
  TString collimator = data->GetCollimator();
  std::vector <double> positions = data->GetPositions();
  data->Print();
  
  SFAttenuation *att;
  try{
    att = new SFAttenuation(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in attenuation.cc!" << std::endl;
    return 1;
  }
  
  //----- averaged channels method
  att->AttAveragedCh();
  att->Fit3rdOrder();
  std::vector <TH1D*> attRatios = att->GetRatios();
  TGraphErrors *attGraph        = att->GetAttGraph();
  TGraphErrors *attGraphPol3    = (TGraphErrors*)attGraph->Clone("attGraphPol3");
  TGraphErrors *sigGraph        = att->GetSigmaGraph();
  
  //----- separate channels method
  att->AttSeparateCh(0);
  TGraphErrors *attGraphCh0       = att->GetAttGraph(0);
  std::vector <TH1D*> spectraCh0  = att->GetSpectra(0);
  
  att->AttSeparateCh(1);
  TGraphErrors *attGraphCh1      = att->GetAttGraph(1);
  std::vector <TH1D*> spectraCh1 = att->GetSpectra(1);

  //----- numeric results
  AttenuationResults results = att->GetResults();
  
  //-----drawing averaged channels
  TLatex text;
  text.SetNDC(true);
  
  TCanvas *can_averaged_ch = new TCanvas("can_averaged_ch", "can_averaged_ch", 1000, 700);
  can_averaged_ch->Divide(2,1);
  
  can_averaged_ch->cd(1);
  gPad->SetGrid(1,1);
  attGraph->SetTitle(Form("Series %i, attenuation curve", seriesNo));
  attGraph->GetFunction("fpol3")->Delete();
  attGraph->Draw("AP");
  text.SetTextSize(0.04);
  text.DrawLatex(0.2, 0.8, Form("L_{att} = (%.2f +/- %.2f) mm", 
                 results.fAttCombPol1, results.fAttCombPol1Err));
  
  can_averaged_ch->cd(2);
  gPad->SetGrid(1,1);
  attGraphPol3->SetTitle(Form("Series %i, attenuation curve", seriesNo));
  attGraphPol3->GetFunction("fpol1")->Delete();
  attGraphPol3->Draw("AP");
  TF1 *funPol3 = (TF1*)attGraphPol3->GetFunction("fpol3");
  text.SetTextSize(0.03);
  text.DrawLatex(0.2, 0.8, Form("A_{0} = %.4e +/- %.4e", 
                funPol3->GetParameter(0), funPol3->GetParError(0)));
  text.DrawLatex(0.2, 0.75, Form("A_{1} = %.4e +/- %.4e",  
                funPol3->GetParameter(1), funPol3->GetParError(1)));
  text.DrawLatex(0.2, 0.70, Form("A_{2} = %.4e +/- %.4e",  
                funPol3->GetParameter(2), funPol3->GetParError(2)));
  text.DrawLatex(0.2, 0.65, Form("A_{3} = %.4e +/- %.4e",  
                funPol3->GetParameter(3), funPol3->GetParError(3)));
  text.DrawLatex(0.2, 0.50, Form("L_{att} = (%.2f +/- %.2f) mm", 
                 results.fAttCombPol3, results.fAttCombPol3Err));

  TCanvas *can_sig = new TCanvas("can_sig", "can_sig", 600, 600);
  gPad->SetGrid(1,1);
  sigGraph->SetTitle(Form("Sigma of M_{LR} distribution S%i", seriesNo));
  sigGraph->Draw("AP");
  
  TCanvas *can_ratios = new TCanvas("can_ratios", "can_ratios", 1200, 1200);
  can_ratios->Divide(3,3);
  TString htitle;
  text.SetTextSize(0.02);
  
  TF1 *fun;
  TF1 *fthin = new TF1("fthin", "gaus", -1, 1);
  TF1 *fthick = new TF1("fthick", "gaus", -1, 1);
  
  for(int i=0; i<npoints; i++){
   can_ratios->cd(i+1);
   gPad->SetGrid(1,1);
   attRatios[i]->SetTitle(Form("ln(#sqrt{ch1/ch0}), source position %.2f mm", positions[i]));
   attRatios[i]->GetXaxis()->SetRangeUser(-2,2);
   attRatios[i]->GetXaxis()->SetTitle("ln(#sqrt{ch1/ch0})");
   attRatios[i]->Draw();
   if(collimator=="Lead"){
     fun = attRatios[i]->GetFunction("fun");
     fthin->SetParameters(fun->GetParameter(0),
                          fun->GetParameter(1),
                          fun->GetParameter(2));
     fthick->SetParameters(fun->GetParameter(3),
                           fun->GetParameter(4),
                           fun->GetParameter(5));
     fthin->SetLineColor(kMagenta);
     fthick->SetLineColor(kMagenta-10);
     fthin->DrawClone("same");
     fthick->DrawClone("same");
   }
  }
  
  //----- drawing separate channels 
  text.SetTextSize(0.025);
  TCanvas *can_separate_ch = new TCanvas("can_separate_ch", "can_separate_ch", 1000, 1000);
  gPad->SetGrid(1,1);
  
  attGraphCh0->SetTitle(Form("Series %i channel 0, attenuation curve", seriesNo));
  attGraphCh0->GetYaxis()->SetTitleSize(0.03);
  attGraphCh0->SetMarkerColor(kRed);
  attGraphCh0->SetLineColor(kRed);
  attGraphCh0->Draw("AP");
  text.SetTextColor(kRed);
  text.DrawLatex(0.3, 0.8, Form("L_{att Ch0} = (%.2f +/- %.2f) mm", 
                 results.fAttCh0, results.fAttCh0Err));
  
  attGraphCh1->SetTitle(Form("Series %i channel 1, attenuation curve", seriesNo));
  attGraphCh1->GetYaxis()->SetTitleSize(0.03);
  attGraphCh1->SetMarkerColor(kGreen+3);
  attGraphCh1->SetLineColor(kGreen+3);
  attGraphCh1->Draw("P");
  attGraphCh1->GetFunction("fexp")->SetLineColor(kGreen+3);
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.3, 0.7, Form("L_{att Ch1} = (%.2f +/- %.2f) mm", 
                 results.fAttCh1, results.fAttCh1Err));
  
  double *yCh0 = attGraphCh0->GetY();
  double *yCh1 = attGraphCh1->GetY();
  double yminCh0 = TMath::MinElement(npoints, yCh0);
  double yminCh1 = TMath::MinElement(npoints, yCh1);
  double ymin = TMath::Min(yminCh0, yminCh1);
  double ymaxCh0 = TMath::MaxElement(npoints, yCh0);
  double ymaxCh1 = TMath::MaxElement(npoints, yCh1);
  double ymax = TMath::Max(ymaxCh0, ymaxCh1);
  
  attGraphCh0->GetYaxis()->SetRangeUser(ymin-0.2*ymin, ymax+0.1*ymax);
  
  TCanvas *can_spectra_ch0 = new TCanvas("can_spectra_ch0", "can_spectra_ch0", 1200, 1200);
  can_spectra_ch0->DivideSquare(npoints);
  
  TCanvas *can_spectra_ch1 = new TCanvas("can_spectra_ch1", "can_spectra_ch1", 1200, 1200);
  can_spectra_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
   can_spectra_ch0->cd(i+1);
   gPad->SetGrid(1,1);
   spectraCh0[i]->SetStats(false);
   spectraCh0[i]->GetXaxis()->SetRangeUser(0,1200); 
   spectraCh0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm", seriesNo, positions[i]));spectraCh0[i]->GetXaxis()->SetTitle("P.E.");
   spectraCh0[i]->GetYaxis()->SetTitle("counts");
   spectraCh0[i]->GetYaxis()->SetMaxDigits(2);
   spectraCh0[i]->Draw();
   
   can_spectra_ch1->cd(i+1);
   gPad->SetGrid(1,1);
   spectraCh1[i]->SetStats(false);
   spectraCh1[i]->GetXaxis()->SetRangeUser(0,1200);
   spectraCh1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm", seriesNo, positions[i]));
   spectraCh1[i]->GetXaxis()->SetTitle("P.E.");
   spectraCh1[i]->GetYaxis()->SetTitle("counts");
   spectraCh1[i]->GetYaxis()->SetMaxDigits(2);
   spectraCh1[i]->Draw();
  }
  
  //----- saving
  
  TString fname = Form("attenuation_series%i.root", seriesNo);
  TString outdir;
  TString dbase;

  int ret = parse_common_options(argc, argv, outdir, dbase);
  if(ret != 0) 
    exit(ret);
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in attenuation.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  can_averaged_ch->Write();
  can_sig->Write();
  can_separate_ch->Write();
  can_ratios->Write();
  can_spectra_ch0->Write();
  can_spectra_ch1->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "ATTENUATION_LENGTH";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ATT_CH0, ATT_CH0_ERR, ATT_CH1, ATT_CH1_ERR, ATT_COMB, ATT_COMB_ERR, ATT_COMB_POL3, ATT_COMB_POL3_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), results.fAttCh0, results.fAttCh0Err, results.fAttCh1, results.fAttCh1Err, results.fAttCombPol1, results.fAttCombPol1Err, results.fAttCombPol3, results.fAttCombPol3Err);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete att;
  
  return 0;
}
