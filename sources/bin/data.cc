// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               data.cc                 *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************  

#include "SFData.hh"
#include "SFDrawCommands.hh"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"

int main(int argc, char **argv){
 
  if(argc!=2){
    std::cout << "to run type: ./data seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 0;
  }
  
  int seriesNo = atoi(argv[1]);
  
  SFData *data;
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in data.cc!" << std::endl;
    return 0;
  }
  
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in data.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cerr << "Description: " << desc << std::endl;
    return 0;
  }
  
  int npoints  = data->GetNpoints();
  std::vector <double> positions = data->GetPositions();
  data->Print(); 
  
  //----- accessing spectra
  std::vector <TH1D*> hAmpCh0 = data->GetSpectra(0, SFSelectionType::Amplitude, "");
  std::vector <TH1D*> hAmpCh1 = data->GetSpectra(1, SFSelectionType::Amplitude, "");
  
  std::vector <TH1D*> hChargeCh0 = data->GetSpectra(0, SFSelectionType::PE, "ch_0.fPE>0");
  std::vector <TH1D*> hChargeCh1 = data->GetSpectra(1, SFSelectionType::PE, "ch_1.fPE>0");
   
  std::vector <TH1D*> hT0Ch0 = data->GetSpectra(0, SFSelectionType::T0, "ch_0.fT0>0");
  std::vector <TH1D*> hT0Ch1 = data->GetSpectra(1, SFSelectionType::T0, "ch_1.fT0>0");
  
  std::vector <TH1D*> hTOTCh0 = data->GetSpectra(0, SFSelectionType::TOT, "ch_0.fTOT>0");
  std::vector <TH1D*> hTOTCh1 = data->GetSpectra(1, SFSelectionType::TOT, "ch_1.fTOT>0");
  
  std::vector <TH2D*> hCorrAmp = data->GetCorrHistograms(SFSelectionType::AmplitudeCorrelation, "");
  std::vector <TH2D*> hCorrPE  = data->GetCorrHistograms(SFSelectionType::PECorrelation, "ch_0.fPE>0 && ch_1.fPE>0");
  std::vector <TH2D*> hCorrT0  = data->GetCorrHistograms(SFSelectionType::T0Correlation, "ch_0.fT0>0 && ch_1.fT0>0");
  
  //----- accessing signals
  const int nsig = 6;
  int number = 0;
  std::vector <TH1D*> hSigCh0(nsig);
  std::vector <TH1D*> hSigCh1(nsig);
  
  for(int i=0; i<nsig/2; i++){
    number = 100*(i+1);
    hSigCh0[i]          = data->GetSignal(0, 30, "", number, true);
    hSigCh0[i+(nsig/2)] = data->GetSignal(0, 90, "", number, true);
    hSigCh1[i]          = data->GetSignal(1, 30, "", number, true);
    hSigCh1[i+(nsig/2)] = data->GetSignal(1, 90, "", number, true);
  }
  
  const int nsigav = 3;
  std::vector <TProfile*> hSigAvCh0(nsigav);
  std::vector <TProfile*> hSigAvCh1(nsigav);
  
  hSigAvCh0[0] = data->GetSignalAverage(0, 50, "ch_0.fPE>59.5 && ch_0.fPE<60.5", 20, true);
  hSigAvCh0[1] = data->GetSignalAverage(0, 50, "ch_0.fPE>199.5 && ch_0.fPE<200.5", 20, true);
  hSigAvCh0[2] = data->GetSignalAverage(0, 50, "ch_0.fPE>399.5 && ch_0.fPE<400.5", 20, true); 
  
  hSigAvCh1[0] = data->GetSignalAverage(1, 50, "ch_1.fPE>59.5 && ch_1.fPE<60.5", 20, true);
  hSigAvCh1[1] = data->GetSignalAverage(1, 50, "ch_1.fPE>199.5 && ch_1.fPE<200.5", 20, true);
  hSigAvCh1[2] = data->GetSignalAverage(1, 50, "ch_1.fPE>399.5 && ch_1.fPE<400.5", 20, true);
  
  //----- drawing spectra
  TCanvas *can_ampl = new TCanvas("can_ampl", "can_ampl", 1500, 1200);
  can_ampl->DivideSquare(npoints);
  
  TCanvas *can_charge = new TCanvas("can_charge", "can_charge", 1500, 1200);
  can_charge->DivideSquare(npoints);
  
  TCanvas *can_t0 = new TCanvas("can_t0", "can_t0", 1500, 1200);
  can_t0->DivideSquare(npoints);
  
  TCanvas *can_tot = new TCanvas("can_tot", "can_tot", 1500, 1200);
  can_tot->DivideSquare(npoints);
  
  TCanvas *can_ampl_corr = new TCanvas("can_ampl_corr", "can_ampl_corr", 1500, 1200);
  can_ampl_corr->DivideSquare(npoints);
  
  TCanvas *can_charge_corr = new TCanvas("can_charge_corr", "can_charge_corr", 1500, 1200);
  can_charge_corr->DivideSquare(npoints);
  
  TCanvas *can_t0_corr = new TCanvas("can_t0_corr", "can_t0_corr", 1500, 1200);
  can_t0_corr->DivideSquare(npoints);
  
  TString stringCh0;
  TString stringCh1;
  TString string;
  TLatex textCh0;
  TLatex textCh1;
  TLatex text;
  textCh0.SetTextSize(0.033);
  textCh0.SetTextColor(kBlue);
  textCh0.SetTextFont(42);
  textCh0.SetNDC(true);
  textCh1.SetTextSize(0.033);
  textCh1.SetTextColor(kRed);
  textCh1.SetTextFont(42);
  textCh1.SetNDC(true);
  text.SetTextSize(0.031);
  text.SetTextColor(kGray+2);
  text.SetTextFont(42);
  text.SetNDC(true);
  
  std::vector <TPaveStats*> paves;
  
  for(int i=0; i<npoints; i++){
    can_ampl->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hAmpCh0[i]->GetTitle();
    stringCh1 = hAmpCh1[i]->GetTitle();
    double maxCh0 = hAmpCh0[i]->GetBinContent(hAmpCh0[i]->GetMaximumBin());
    double maxCh1 = hAmpCh1[i]->GetBinContent(hAmpCh1[i]->GetMaximumBin());
    double maxYaxis = std::max(maxCh0 ,maxCh1);
    maxYaxis += maxYaxis*0.1;
    hAmpCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    hAmpCh0[i]->SetTitle(Form("Amplitude spectrum, source position %.2f mm", positions[i]));
    hAmpCh0[i]->GetXaxis()->SetTitle("signal amplitude [mV]");
    hAmpCh0[i]->GetYaxis()->SetTitle("counts");
    hAmpCh0[i]->GetYaxis()->SetMaxDigits(2);
    hAmpCh0[i]->SetStats(false);
    hAmpCh0[i]->SetLineColor(kBlue);
    hAmpCh1[i]->SetLineColor(kRed);
    hAmpCh1[i]->SetStats(false);
    hAmpCh0[i]->Draw();
    hAmpCh1[i]->Draw("same");
    textCh0.DrawLatex(0.3, 0.8, stringCh0);
    textCh1.DrawLatex(0.3, 0.75, stringCh1);
    
    can_charge->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hChargeCh0[i]->GetTitle();
    stringCh1 = hChargeCh1[i]->GetTitle();
    maxCh0 = hChargeCh0[i]->GetBinContent(hChargeCh0[i]->GetMaximumBin());
    maxCh1 = hChargeCh1[i]->GetBinContent(hChargeCh1[i]->GetMaximumBin());
    maxYaxis = std::max(maxCh0, maxCh1);
    maxYaxis += maxYaxis*0.1;
    hChargeCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    hChargeCh0[i]->SetTitle(Form("Charge spectrum, source position %.2f mm", positions[i]));
    hChargeCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hChargeCh0[i]->GetYaxis()->SetTitle("counts");
    hChargeCh0[i]->GetYaxis()->SetMaxDigits(2);
    hChargeCh0[i]->GetXaxis()->SetRangeUser(0, 1300);
    hChargeCh0[i]->SetStats(false);
    hChargeCh0[i]->SetLineColor(kBlue);
    hChargeCh1[i]->SetStats(false);
    hChargeCh1[i]->SetLineColor(kRed);
    hChargeCh0[i]->Draw();
    hChargeCh1[i]->Draw("same");
    textCh0.DrawLatex(0.3, 0.8, stringCh0);
    textCh1.DrawLatex(0.3, 0.75, stringCh1);
    
    can_t0->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hT0Ch0[i]->GetTitle();
    stringCh1 = hT0Ch0[i]->GetTitle();
    maxCh0 = hT0Ch0[i]->GetBinContent(hT0Ch0[i]->GetMaximumBin());
    maxCh1 = hT0Ch1[i]->GetBinContent(hT0Ch1[i]->GetMaximumBin());
    maxYaxis = std::max(maxCh0, maxCh1);
    maxYaxis += maxYaxis*0.1;
    hT0Ch0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    hT0Ch0[i]->SetTitle(Form("T0 spectrum, source position %.2f mm", positions[i]));
    hT0Ch0[i]->GetXaxis()->SetTitle("time [ns]");
    hT0Ch0[i]->GetXaxis()->SetRangeUser(-50, 400);
    hT0Ch0[i]->SetLineColor(kBlue);
    hT0Ch1[i]->SetLineColor(kRed);
    hT0Ch0[i]->Draw();
    gPad->Update();
    paves.push_back((TPaveStats*)hT0Ch0[i]->FindObject("stats"));
    if(paves[i]==NULL) std::cout << "Warning " << i << std::endl;
    paves[i]->SetY1NDC(0.6);
    paves[i]->SetY2NDC(0.75);
    hT0Ch1[i]->Draw("sames");
    gPad->Update();
    textCh0.DrawLatex(0.4, 0.3, stringCh0);
    textCh1.DrawLatex(0.4, 0.25, stringCh1);
    
    can_tot->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hTOTCh0[i]->GetTitle();
    stringCh1 = hTOTCh1[i]->GetTitle();
    maxCh0 = hTOTCh0[i]->GetBinContent(hTOTCh0[i]->GetMaximumBin());
    maxCh1 = hTOTCh1[i]->GetBinContent(hTOTCh1[i]->GetMaximumBin());
    maxYaxis = std::max(maxCh0, maxCh1);
    maxYaxis += maxYaxis*0.1; 
    hTOTCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    hTOTCh0[i]->SetTitle(Form("TOT spectrum, source position %.2f mm", positions[i]));
    hTOTCh0[i]->GetXaxis()->SetTitle("time [ns]");
    hTOTCh0[i]->SetStats(false);
    hTOTCh0[i]->SetLineColor(kBlue);
    hTOTCh1[i]->SetStats(false);
    hTOTCh1[i]->SetLineColor(kRed);
    hTOTCh0[i]->Draw();
    hTOTCh1[i]->Draw("same");
    textCh0.DrawLatex(0.3, 0.8, stringCh0);
    textCh1.DrawLatex(0.3, 0.75, stringCh1);
    
    can_ampl_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrAmp[i]->GetTitle();
    hCorrAmp[i]->SetTitle(Form("Amplitude correlation spectrum, source position %.2f mm", positions[i]));
    hCorrAmp[i]->GetXaxis()->SetTitle("Ch0 amplitude [mV]");
    hCorrAmp[i]->GetYaxis()->SetTitle("Ch1 amplitude [mV]");
    hCorrAmp[i]->SetStats(false);
    hCorrAmp[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
    
    can_charge_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrPE[i]->GetTitle();
    hCorrPE[i]->SetTitle(Form("Charge correlation spectrum, source position %.2f mm", positions[i]));
    hCorrPE[i]->GetXaxis()->SetTitle("Ch0 charge [P.E.]");
    hCorrPE[i]->GetYaxis()->SetTitle("Ch1 charge [P.E.]");
    hCorrPE[i]->GetXaxis()->SetRangeUser(0, 1300);
    hCorrPE[i]->GetYaxis()->SetRangeUser(0, 1300);
    hCorrPE[i]->SetStats(false);
    hCorrPE[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
    
    can_t0_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrT0[i]->GetTitle();
    hCorrT0[i]->SetTitle(Form("T0 correlation spectrum, source position %.2f mm", positions[i]));
    hCorrT0[i]->GetXaxis()->SetTitle("Ch0 T0 [ns]");
    hCorrT0[i]->GetYaxis()->SetTitle("Ch1 T0 [ns]");
    hCorrT0[i]->SetStats(false);
    hCorrT0[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
  }
  
  //----- drawing signals
  TCanvas *can_sig = new TCanvas("can_sig","can_sig",1200,800);
  can_sig->Divide(3,2);
  
  for(int i=0; i<nsig; i++){
    can_sig->cd(i+1);
    gPad->SetGrid(1,1);
    double maxCh0 = hSigCh0[i]->GetBinContent(hSigCh0[i]->GetMaximumBin());
    double maxCh1 = hSigCh1[i]->GetBinContent(hSigCh1[i]->GetMaximumBin());
    double maxYaxis = std::max(maxCh0,maxCh1) + 10.;
    stringCh0 = hSigCh0[i]->GetTitle();
    stringCh1 = hSigCh1[i]->GetTitle();
    hSigCh0[i]->SetLineColor(kBlue);
    hSigCh0[i]->SetTitle(" ");
    hSigCh0[i]->GetXaxis()->SetTitle("time [ns]");
    hSigCh0[i]->GetYaxis()->SetTitle("amplitude [mV]");
    hSigCh0[i]->GetYaxis()->SetRangeUser(-2,maxYaxis);
    hSigCh0[i]->SetStats(false);
    hSigCh1[i]->SetLineColor(kRed);
    hSigCh1[i]->SetStats(false);
    hSigCh0[i]->Draw();
    hSigCh1[i]->Draw("same");
    textCh0.DrawLatex(0.4,0.8,stringCh0);
    textCh1.DrawLatex(0.4,0.75,stringCh1);
  }
  
  TCanvas *can_sigav = new TCanvas("can_sigav","can_sigav",1200,800);
  can_sigav->Divide(3,2);
  
  textCh0.SetTextColor(kGray+2);
  textCh0.SetTextSize(0.025);
  textCh1.SetTextColor(kGray+2);
  textCh1.SetTextSize(0.025);
  
  for(int i=0; i<nsigav; i++){
    can_sigav->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hSigAvCh0[i]->GetTitle();
    hSigAvCh0[i]->SetTitle(" ");
    hSigAvCh0[i]->GetXaxis()->SetTitle("time [ns]");
    hSigAvCh0[i]->GetYaxis()->SetTitle("amplitude [mV]");
    double maxYaxis = hSigAvCh0[i]->GetBinContent(hSigAvCh0[i]->GetMaximumBin());
    maxYaxis = maxYaxis+0.2*maxYaxis;
    hSigAvCh0[i]->GetYaxis()->SetRangeUser(-10,maxYaxis);
    hSigAvCh0[i]->SetStats(false);
    hSigAvCh0[i]->Draw();
    textCh0.DrawLatex(0.15,0.85,stringCh0);
    
    can_sigav->cd(i+1+nsigav);
    gPad->SetGrid(1,1);
    stringCh1 = hSigAvCh1[i]->GetTitle();
    hSigAvCh1[i]->SetTitle(" ");
    hSigAvCh1[i]->GetXaxis()->SetTitle("time [ns]");
    hSigAvCh1[i]->GetYaxis()->SetTitle("amplitude [mV]");
    maxYaxis = hSigAvCh1[i]->GetBinContent(hSigAvCh1[i]->GetMaximumBin());
    maxYaxis = maxYaxis+0.2*maxYaxis;
    hSigAvCh1[i]->GetYaxis()->SetRangeUser(-10,maxYaxis);
    hSigAvCh1[i]->SetStats(false);
    hSigAvCh1[i]->Draw();
    textCh1.DrawLatex(0.15,0.85,stringCh1);
  }
  
  //----- saving
  TString path = std::string(getenv("SFPATH"));
  TString fname = Form("data_series%i.root", seriesNo);
  
  CmdLineOption cmd_outdir("Output directory", "-out", "Output directory (string), default: $SFPATH/results" , path+"results");
  
  CmdLineOption cmd_dbase("Data base", "-db", "Data base name (string), default: ScintFibRes.db", "ScintFibRes.db");
  
  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
  
  TString outdir = CmdLineOption::GetStringValue("Output directory");
  TString dbase = CmdLineOption::GetStringValue("Data base");
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");

  if(!file->IsOpen()){
    std::cerr << "##### Error in data.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 0;
  }
  
  can_ampl->Write();
  can_charge->Write();
  can_t0->Write();
  can_tot->Write();
  can_ampl_corr->Write();
  can_charge_corr->Write();
  can_t0_corr->Write();
  can_sig->Write();
  can_sigav->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "DATA";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES (%i, '%s')", table.Data(), seriesNo, fname_full.Data());
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  
  return 1;
}
