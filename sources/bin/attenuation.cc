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
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv){
  
  if(argc<2 || argc>6){
    std::cout << "to run type: ./attenuation seriesNo";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 0;
  }
 
  int seriesNo = atoi(argv[1]);

  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char* message){
   std::cerr << message << std::endl;
   std::cerr << "##### Exception in attenuation.cc!" << std::endl;
   return 0;
  }
 
  data->Print();
 
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in attenuation.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cout << "Description: " << desc << std::endl;
    return 0;
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
    return 0;
  }
  
  //----- averaged channels method
  att->AttAveragedCh();
  std::vector <TH1D*> attRatios = att->GetRatios();
  TGraphErrors *attGraph        = att->GetAttGraph();
  std::vector <double> attlen   = att->GetAttenuation();
  
  //----- separate channels method
  att->AttSeparateCh(0);
  TGraphErrors *attGraphCh0       = att->GetAttGraph(0);
  std::vector <TH1D*> spectraCh0  = att->GetSpectra(0);
  std::vector <double> attlenCh0  = att->GetAttenuation(0);
  std::vector <TH1D*> peaksCh0;   
  
  att->AttSeparateCh(1);
  TGraphErrors *attGraphCh1      = att->GetAttGraph(1);
  std::vector <TH1D*> spectraCh1 = att->GetSpectra(1);
  std::vector <double> attlenCh1 = att->GetAttenuation(1);  
  std::vector <TH1D*> peaksCh1;  
  
  if(collimator=="Lead"){
      peaksCh0 = att->GetPeaks(0);
      peaksCh1 = att->GetPeaks(1);
  }
  
  //-----drawing averaged channels
  TLatex text;
  text.SetNDC(true);
  
  TCanvas *can_averaged_ch = new TCanvas("can_averaged_ch", "can_averaged_ch", 700, 500);
  can_averaged_ch->cd();
  gPad->SetGrid(1,1);
  attGraph->SetTitle(Form("Series %i, attenuation curve", seriesNo));
  attGraph->Draw("AP");
  text.SetTextSize(0.04);
  text.DrawLatex(0.2, 0.8, Form("L_{att} = (%.2f +/- %.2f) mm", attlen[0], attlen[1]));
  
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
  text.SetTextSize(0.04);
  TCanvas *can_separate_ch = new TCanvas("can_separate_ch", "can_separate_ch", 1000, 1000);
  gPad->SetGrid(1,1);
  attGraphCh0->SetTitle(Form("Series %i channel 0, attenuation curve", seriesNo));
  attGraphCh0->GetYaxis()->SetTitleSize(0.03);
  attGraphCh0->SetMarkerColor(kRed);
  attGraphCh0->SetLineColor(kRed);
  attGraphCh0->Draw("AP");
  text.SetTextColor(kRed);
  text.DrawLatex(0.3, 0.8, Form("L_{att Ch0} = (%.2f +/- %.2f) mm", attlenCh0[0], attlenCh0[1]));
  attGraphCh1->SetTitle(Form("Series %i channel 1, attenuation curve", seriesNo));
  attGraphCh1->GetYaxis()->SetTitleSize(0.03);
  attGraphCh1->SetMarkerColor(kGreen+3);
  attGraphCh1->SetLineColor(kGreen+3);
  attGraphCh1->Draw("P");
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.3, 0.7, Form("L_{att Ch1} = (%.2f +/- %.2f) mm", attlenCh1[0], attlenCh1[1]));
  
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
   if(collimator=="Lead"){
       peaksCh0[i]->Draw("same");
   }
   can_spectra_ch1->cd(i+1);
   gPad->SetGrid(1,1);
   spectraCh1[i]->SetStats(false);
   spectraCh1[i]->GetXaxis()->SetRangeUser(0,1200);
   spectraCh1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm", seriesNo, positions[i]));
   spectraCh1[i]->GetXaxis()->SetTitle("P.E.");
   spectraCh1[i]->GetYaxis()->SetTitle("counts");
   spectraCh1[i]->GetYaxis()->SetMaxDigits(2);
   spectraCh1[i]->Draw();
   if(collimator=="Lead"){
     peaksCh1[i]->Draw("same");
   }
  }
  
  //----- saving
  TString path = std::string(getenv("SFPATH"));
  TString fname = Form("attenuation_series%i.root", seriesNo);

  CmdLineOption cmd_outdir("Output directory", "-out", "Output directory (string), default: $SFPATH/results", path+"results");
  
  CmdLineOption cmd_dbase("Data base", "-db", "Data base name (string), default: ScintFibRes.db", "ScintFibRes.db");
  
  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
  
  TString outdir = CmdLineOption::GetStringValue("Output directory");
  TString dbase = CmdLineOption::GetStringValue("Data base");
  
  if(!gSystem->ChangeDirectory(outdir)){
    std::cout << "Creating new directory... " << std::endl;
    std::cout << outdir << std::endl;
    int stat = mkdir(outdir, 0777);
    if(stat==-1){
      std::cerr << "##### Error in attenuation.cc! Unable to create new direcotry!" << std::endl;
      return 0;
    }
  }
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in attenuation.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 0;
  }
  
  can_averaged_ch->Write();
  can_separate_ch->Write();
  can_ratios->Write();
  can_spectra_ch0->Write();
  can_spectra_ch1->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "ATTENUATION_LENGTH";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ATT_CH0, ATT_CH0_ERR, ATT_CH1, ATT_CH1_ERR, ATT_COMB, ATT_COMB_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), attlenCh0[0], attlenCh0[1], attlenCh1[0], attlenCh1[1], attlen[0], attlen[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete att;
  
  return 1;
}
