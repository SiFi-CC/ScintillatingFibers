// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               data.cc                 *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************  

#include "common_options.h"
#include "SFData.hh"
#include "SFDrawCommands.hh"

#include <DistributionContext.h>

#include <TCanvas.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TF1.h>

#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv){
    
  TString outdir;
  TString dbase;
  int seriesNo = -1;

  int ret = parse_common_options(argc, argv, outdir, dbase, seriesNo);
  if(ret != 0) 
    exit(ret);
 
  if(argc<2){
    std::cout << "to run type: ./data seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }
  
  SFData *data;
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in data.cc!" << std::endl;
    return 1;
  }
  
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Warning in data.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cerr << "Description: " << desc << std::endl;
  }
  
  int npoints  = data->GetNpoints();
  int anaGroup = data->GetAnalysisGroup();
  std::vector <double> positions = data->GetPositions();
  std::vector <int> measurementsIDs = data->GetMeasurementsIDs();
  TString collimator = data->GetCollimator();
  TString fiber = data->GetFiber();
  data->Print(); 
  
  DistributionContext ctx;
  ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

  ctx.dim = DIM1;
  ctx.x.min = 0;
  ctx.x.max = 100;
  ctx.y.min = 0;
  ctx.y.max = 1000;
  
  //gStyle->SetPalette(53);
  
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
  std::vector <TH2D*> hAmpPECh0 = data->GetCorrHistograms(SFSelectionType::AmpPECorrelation, "ch_0.fPE>0", 0);
  std::vector <TH2D*> hAmpPECh1 = data->GetCorrHistograms(SFSelectionType::AmpPECorrelation, "ch_1.fPE>0", 1);
  
  std::vector <TH1D*> hChargeCh2;
  std::vector <TH2D*> hChargeCh0Ch2;
  std::vector <TH2D*> hChargeCh1Ch2;
  
  if(collimator.Contains("Electronic")){
    hChargeCh2 = data->GetSpectra(2, SFSelectionType::Charge, "ch_2.fCharge>0"); 
    hChargeCh0Ch2 = data->GetCorrHistograms(SFSelectionType::PEvsPEch2Correlation, "", 0);
    hChargeCh1Ch2 = data->GetCorrHistograms(SFSelectionType::PEvsPEch2Correlation, "", 1);
  }
  
  //----- accessing signals
  const int nsig = 6;
  int number = 0;
  std::vector <TH1D*> hSigCh0(nsig);
  std::vector <TH1D*> hSigCh1(nsig);
  
  for(int i=0; i<nsig/2; i++){
    number = 100*(i+1);
    hSigCh0[i]          = data->GetSignal(0, measurementsIDs[2], "", number, true);
    hSigCh0[i+(nsig/2)] = data->GetSignal(0, measurementsIDs[npoints-1], "", number, true);
    hSigCh1[i]          = data->GetSignal(1, measurementsIDs[2], "", number, true);
    hSigCh1[i+(nsig/2)] = data->GetSignal(1, measurementsIDs[npoints-1], "", number, true);
  }
  
  const int nsigav = 3;
  std::vector <TProfile*> hSigAvCh0(nsigav);
  std::vector <TProfile*> hSigAvCh1(nsigav);
  
  double PE[3];
  if(collimator.Contains("Electronic")){
    if(fiber.Contains("LYSO") || fiber.Contains("GAGG")){
      PE[0] = 150.; 
      PE[1] = 200.;
      PE[2] = 300.;
    }
    else if(fiber.Contains("LuAG")){
      PE[0] = 50.;
      PE[1] = 100.;
      PE[2] = 150.;
    }
    else{
        std::cerr << "##### Error in data.cc! Unknown fiber material!" << std::endl;
        std::abort();
    }
  }
  else if(collimator.Contains("Lead")){
    PE[0] = 60.;
    PE[1] = 20.;
    PE[2] = 400.;
  }
  
  int ID = SFTools::GetMeasurementID(seriesNo, 50.0);
  //int ID = SFTools::GetMeasurementID(seriesNo, 12.0);
  
  hSigAvCh0[0] = data->GetSignalAverage(0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", 
                                        PE[0]-0.5, PE[0]+0.5), 20, true);
  hSigAvCh0[1] = data->GetSignalAverage(0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", 
                                        PE[1]-0.5, PE[1]+0.5), 20, true);
  hSigAvCh0[2] = data->GetSignalAverage(0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", 
                                        PE[2]-0.5, PE[2]+0.5), 20, true); 
  
  hSigAvCh1[0] = data->GetSignalAverage(1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", 
                                        PE[0]-0.5, PE[0]+0.5), 20, true);
  hSigAvCh1[1] = data->GetSignalAverage(1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", 
                                        PE[1]-0.5, PE[1]+0.5), 20, true);
  hSigAvCh1[2] = data->GetSignalAverage(1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", 
                                        PE[2]-0.5, PE[2]+0.5), 20, true);
  
  //----- drawing spectra
  TCanvas *can_ampl = new TCanvas("data_ampl", "data_ampl", 2000, 1200);
  can_ampl->DivideSquare(npoints);
  
  TCanvas *can_charge = new TCanvas("data_charge", "data_charge", 2000, 1200);
  can_charge->DivideSquare(npoints);
  
  TCanvas *can_t0 = new TCanvas("data_t0", "data_t0", 2000, 1200);
  can_t0->DivideSquare(npoints);
  
  TCanvas *can_tot = new TCanvas("data_tot", "data_tot", 2000, 1200);
  can_tot->DivideSquare(npoints);
  
  TCanvas *can_ampl_corr = new TCanvas("data_ampl_corr", "data_ampl_corr", 2000, 1200);
  can_ampl_corr->DivideSquare(npoints);
  
  TCanvas *can_charge_corr = new TCanvas("data_charge_corr", "data_charge_corr", 2000, 1200);
  can_charge_corr->DivideSquare(npoints);
  
  TCanvas *can_amp_pe_ch0 = new TCanvas("data_amp_pe_ch0", "data_amp_pe_ch0", 2000, 1200);
  can_amp_pe_ch0->DivideSquare(npoints);
  
  TCanvas *can_amp_pe_ch1 = new TCanvas("data_amp_pe_ch1", "data_amp_pe_ch1", 2000, 1200);
  can_amp_pe_ch1->DivideSquare(npoints);
  
  TCanvas *can_t0_corr = new TCanvas("data_t0_corr", "data_t0_corr", 2000, 1200);
  can_t0_corr->DivideSquare(npoints);
  
  TCanvas *can_ref;
  TCanvas *can_ref_ch0;
  TCanvas *can_ref_ch1;
  
  if(collimator.Contains("Electronic")){ 
      
    can_ref = new TCanvas("data_ref", "data_ref", 2000, 1200);
    can_ref->DivideSquare(npoints);
    
    can_ref_ch0 = new TCanvas("data_ref_ch0", "data_ref_ch0", 2000, 1200);
    can_ref_ch0->DivideSquare(npoints);
    
    can_ref_ch1 = new TCanvas("data_ref_ch1", "data_ref_ch1", 2000, 1200);
    can_ref_ch1->DivideSquare(npoints);
  }
  
  int colCh0 = kPink-8;
  int colCh1 = kAzure-6;
  
  TString stringCh0;
  TString stringCh1;
  TString string;
  TLatex textCh0;
  TLatex textCh1;
  TLatex text;
  textCh0.SetTextSize(0.033);
  textCh0.SetTextColor(colCh0);
  textCh0.SetTextFont(42);
  textCh0.SetNDC(true);
  textCh1.SetTextSize(0.033);
  textCh1.SetTextColor(colCh1);
  textCh1.SetTextFont(42);
  textCh1.SetNDC(true);
  text.SetTextSize(0.031);
  text.SetTextColor(kGray+2);
  text.SetTextFont(42);
  text.SetNDC(true);
  
  std::vector <TPaveStats*> paves;
  
  TF1 *fdiag; 
  
  for(int i=0; i<npoints; i++){
    can_ampl->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hAmpCh0[i]->GetTitle();
    stringCh1 = hAmpCh1[i]->GetTitle();
    //double maxCh0 = hAmpCh0[i]->GetBinContent(hAmpCh0[i]->GetMaximumBin());
    //double maxCh1 = hAmpCh1[i]->GetBinContent(hAmpCh1[i]->GetMaximumBin());
    //double maxYaxis = std::max(maxCh0 ,maxCh1);
    //maxYaxis += maxYaxis*0.1;
    //hAmpCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    ctx.configureFromJson("hAmpCh0");
    ctx.print();
    hAmpCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hAmpCh0[i]->SetTitle(Form("Amplitude spectrum, source position %.2f mm", positions[i]));
    hAmpCh0[i]->GetXaxis()->SetTitle("signal amplitude [mV]");
    hAmpCh0[i]->GetYaxis()->SetTitle("counts");
    hAmpCh0[i]->GetYaxis()->SetMaxDigits(2);
    hAmpCh0[i]->SetStats(false);
    hAmpCh0[i]->SetLineColor(colCh0);
    hAmpCh1[i]->SetLineColor(colCh1);
    hAmpCh1[i]->SetStats(false);
    hAmpCh0[i]->Draw();
    hAmpCh1[i]->Draw("same");
    textCh0.DrawLatex(0.3, 0.8, stringCh0);
    textCh1.DrawLatex(0.3, 0.75, stringCh1);
    
    can_charge->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hChargeCh0[i]->GetTitle();
    stringCh1 = hChargeCh1[i]->GetTitle();
    //maxCh0 = hChargeCh0[i]->GetBinContent(hChargeCh0[i]->GetMaximumBin());
    //maxCh1 = hChargeCh1[i]->GetBinContent(hChargeCh1[i]->GetMaximumBin());
    //maxYaxis = std::max(maxCh0, maxCh1);
    //maxYaxis += maxYaxis*0.1;
    //hChargeCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    ctx.configureFromJson("hChargeCh0");
    ctx.print();
    hChargeCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hChargeCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hChargeCh0[i]->SetTitle(Form("Charge spectrum, source position %.2f mm", positions[i]));
    hChargeCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hChargeCh0[i]->GetYaxis()->SetTitle("counts");
    hChargeCh0[i]->GetYaxis()->SetMaxDigits(2);
    hChargeCh0[i]->SetStats(false);
    hChargeCh0[i]->SetLineColor(colCh0);
    hChargeCh1[i]->SetStats(false);
    hChargeCh1[i]->SetLineColor(colCh1);
    hChargeCh0[i]->Draw();
    hChargeCh1[i]->Draw("same");
    textCh0.DrawLatex(0.3, 0.8, stringCh0);
    textCh1.DrawLatex(0.3, 0.75, stringCh1);
  
    if(collimator.Contains("Electronic")){
      can_ref->cd(i+1);
      gPad->SetGrid(1,1);
      hChargeCh2[i]->SetStats(false);
      string = hChargeCh2[i]->GetTitle();
      hChargeCh2[i]->GetXaxis()->SetTitle("charge [a.u.]");
      hChargeCh2[i]->GetYaxis()->SetTitle("counts");
      hChargeCh2[i]->GetXaxis()->SetRangeUser(0, 120E3);
      hChargeCh2[i]->SetTitle(Form("Charge spectrum, reference detector, position %.2f mm", positions[i]));
      hChargeCh2[i]->Draw();
      text.DrawLatex(0.3, 0.8, string);
      
      can_ref_ch0->cd(i+1);
      gPad->SetGrid(1,1);
      string = hChargeCh0Ch2[i]->GetTitle();
      hChargeCh0Ch2[i]->SetTitle(Form("Charge correlation spectrum Ch2 vs. Ch0, source position %.2f mm", positions[i]));
      hChargeCh0Ch2[i]->GetXaxis()->SetTitle("Ch2 charge [PE]");
      hChargeCh0Ch2[i]->GetYaxis()->SetTitle("Ch0 charge [a.u.]");
      ctx.configureFromJson("hChargeChXCh2");
      ctx.print();
      hChargeCh0Ch2[i]->GetXaxis()->SetRangeUser(0, 120E3);
      hChargeCh0Ch2[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
      hChargeCh0Ch2[i]->SetStats(false);
      hChargeCh0Ch2[i]->Draw("colz");
      text.DrawLatex(0.15, 0.85, string);
      
      can_ref_ch1->cd(i+1);
      gPad->SetGrid(1,1);
      string = hChargeCh1Ch2[i]->GetTitle();
      hChargeCh1Ch2[i]->SetTitle(Form("Charge correlation spectrum Ch2 vs. Ch1, source position %.2f mm", positions[i]));
      hChargeCh1Ch2[i]->GetXaxis()->SetTitle("Ch2 charge [PE]");
      hChargeCh1Ch2[i]->GetYaxis()->SetTitle("Ch1 charge [a.u.]");
      hChargeCh1Ch2[i]->GetXaxis()->SetRangeUser(0, 120E3);
      hChargeCh1Ch2[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
      hChargeCh1Ch2[i]->SetStats(false);
      hChargeCh1Ch2[i]->Draw("colz");
      text.DrawLatex(0.15, 0.85, string);
    }
    
    can_t0->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hT0Ch0[i]->GetTitle();
    stringCh1 = hT0Ch0[i]->GetTitle();
    double maxCh0 = hT0Ch0[i]->GetBinContent(hT0Ch0[i]->GetMaximumBin());
    double maxCh1 = hT0Ch1[i]->GetBinContent(hT0Ch1[i]->GetMaximumBin());
    double maxYaxis = std::max(maxCh0, maxCh1);
    maxYaxis += maxYaxis*0.1;
    hT0Ch0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    hT0Ch0[i]->SetTitle(Form("T_{0} spectrum, source position %.2f mm", positions[i]));
    hT0Ch0[i]->GetXaxis()->SetTitle("time [ns]");
    hT0Ch0[i]->GetYaxis()->SetTitle("counts");
    hT0Ch0[i]->GetXaxis()->SetRangeUser(150, 350);
    hT0Ch0[i]->SetLineColor(colCh0);
    hT0Ch1[i]->SetLineColor(colCh1);
    hT0Ch0[i]->Draw();
    gPad->Update();
    paves.push_back((TPaveStats*)hT0Ch0[i]->FindObject("stats"));
    if(paves[i]==nullptr) std::cout << "Warning " << i << std::endl;
    paves[i]->SetY1NDC(0.55);
    paves[i]->SetY2NDC(0.71);
    hT0Ch1[i]->Draw("sames");
    gPad->Update();
    textCh0.DrawLatex(0.2, 0.3, stringCh0);
    textCh1.DrawLatex(0.2, 0.25, stringCh1);
    
    can_tot->cd(i+1);
    gPad->SetGrid(1,1);
    stringCh0 = hTOTCh0[i]->GetTitle();
    stringCh1 = hTOTCh1[i]->GetTitle();
    //maxCh0 = hTOTCh0[i]->GetBinContent(hTOTCh0[i]->GetMaximumBin());
    //maxCh1 = hTOTCh1[i]->GetBinContent(hTOTCh1[i]->GetMaximumBin());
    //maxYaxis = std::max(maxCh0, maxCh1);
    //maxYaxis += maxYaxis*0.1; 
    //hTOTCh0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
    ctx.configureFromJson("hTOTCh0");
    ctx.print();
    hTOTCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hTOTCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hTOTCh0[i]->SetTitle(Form("TOT spectrum, source position %.2f mm", positions[i]));
    hTOTCh0[i]->GetXaxis()->SetTitle("time [ns]");
    hTOTCh0[i]->GetYaxis()->SetTitle("counts");
    hTOTCh0[i]->SetStats(false);
    hTOTCh0[i]->SetLineColor(colCh0);
    hTOTCh1[i]->SetStats(false);
    hTOTCh1[i]->SetLineColor(colCh1);
    hTOTCh0[i]->Draw();
    hTOTCh1[i]->Draw("same");
    textCh0.DrawLatex(0.4, 0.8, stringCh0);
    textCh1.DrawLatex(0.4, 0.75, stringCh1);
    
    can_ampl_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrAmp[i]->GetTitle();
    hCorrAmp[i]->SetTitle(Form("Amplitude correlation spectrum, source position %.2f mm", positions[i]));
    hCorrAmp[i]->GetXaxis()->SetTitle("Ch1 amplitude [mV]");
    hCorrAmp[i]->GetYaxis()->SetTitle("Ch0 amplitude [mV]");
    hCorrAmp[i]->SetStats(false);
    hCorrAmp[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
    
    can_charge_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrPE[i]->GetTitle();
    hCorrPE[i]->SetTitle(Form("Charge correlation spectrum, source position %.2f mm", positions[i]));
    hCorrPE[i]->GetXaxis()->SetTitle("Ch1 charge [P.E.]");
    hCorrPE[i]->GetYaxis()->SetTitle("Ch0 charge [P.E.]");
    ctx.configureFromJson("hCorrPE");
    ctx.print();
    hCorrPE[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hCorrPE[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hCorrPE[i]->SetStats(false);
    hCorrPE[i]->Draw("colz");
    fdiag = new TF1("fdiag","x[0]",0,ctx.x.max);
    fdiag->Draw("same");
    text.DrawLatex(0.15, 0.85, string);
    
    can_t0_corr->cd(i+1);
    gPad->SetGrid(1,1);
    string = hCorrT0[i]->GetTitle();
    hCorrT0[i]->SetTitle(Form("T0 correlation spectrum, source position %.2f mm", positions[i]));
    hCorrT0[i]->GetXaxis()->SetTitle("Ch1 T0 [ns]");
    hCorrT0[i]->GetYaxis()->SetTitle("Ch0 T0 [ns]");
    hCorrT0[i]->GetXaxis()->SetRangeUser(0,400);
    hCorrT0[i]->GetYaxis()->SetRangeUser(0,400);
    hCorrT0[i]->SetStats(false);
    hCorrT0[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
    
    can_amp_pe_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    string = hAmpPECh0[i]->GetTitle();
    hAmpPECh0[i]->SetTitle(Form("Amplitude vs. Charge correlation spectrum Ch0, source position %.2f mm", positions[i]));
    hAmpPECh0[i]->GetXaxis()->SetTitle("Charge [PE]");
    hAmpPECh0[i]->GetYaxis()->SetTitle("Amplitude [mV]");
    hAmpPECh0[i]->GetXaxis()->SetRangeUser(-10, ctx.x.max);
    hAmpPECh0[i]->GetYaxis()->SetRangeUser(-10, 800);
    hAmpPECh0[i]->SetStats(false);
    hAmpPECh0[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
    
    can_amp_pe_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    string = hAmpPECh1[i]->GetTitle();
    hAmpPECh1[i]->SetTitle(Form("Amplitude vs. Charge correlation spectrum Ch1, source position %.2f mm", positions[i]));
    hAmpPECh1[i]->GetXaxis()->SetTitle("Charge [PE]");
    hAmpPECh1[i]->GetYaxis()->SetTitle("Amplitude [mV]");
    hAmpPECh1[i]->GetXaxis()->SetRangeUser(-10, ctx.x.max);
    hAmpPECh1[i]->GetYaxis()->SetRangeUser(-10, 800);
    hAmpPECh1[i]->SetStats(false);
    hAmpPECh1[i]->Draw("colz");
    text.DrawLatex(0.15, 0.85, string);
  }
  
  //----- drawing signals
  TCanvas *can_sig = new TCanvas("data_sig","data_sig",1800,800);
  can_sig->Divide(3,2);
  
  int polarity = 0;
  double min = hSigAvCh0[0]->GetBinContent(hSigAvCh0[0]->GetMinimumBin());
  double max = hSigAvCh0[0]->GetBinContent(hSigAvCh0[0]->GetMaximumBin());
  
  if(fabs(min)>max){
      polarity = -1;
      std::cout << "Signals are negative, am I right?" << std::endl;
  }
  else{
      polarity = 1;
      std::cout << "Signals are positive, am I right?" << std::endl;
  }
  
  for(int i=0; i<nsig; i++){
    can_sig->cd(i+1);
    gPad->SetGrid(1,1);
    
    stringCh0 = hSigCh0[i]->GetTitle();
    stringCh1 = hSigCh1[i]->GetTitle();
    hSigCh0[i]->SetLineColor(colCh0);
    hSigCh0[i]->SetTitle(" ");
    hSigCh0[i]->GetXaxis()->SetTitle("time [ns]");
    hSigCh0[i]->GetYaxis()->SetTitle("amplitude [mV]");
    hSigCh0[i]->SetStats(false);
    hSigCh1[i]->SetLineColor(colCh1);
    hSigCh1[i]->SetStats(false);
    hSigCh0[i]->Draw();
    hSigCh1[i]->Draw("same");
    if(polarity==1){
      double maxCh0 = hSigCh0[i]->GetBinContent(hSigCh0[i]->GetMaximumBin());
      double maxCh1 = hSigCh1[i]->GetBinContent(hSigCh1[i]->GetMaximumBin());
      double maxYaxis = std::max(maxCh0, maxCh1) + 10.;
      hSigCh0[i]->GetYaxis()->SetRangeUser(-2, maxYaxis);
    }
    else if(polarity==-1){
      double minCh0 = hSigCh0[i]->GetBinContent(hSigCh0[i]->GetMinimumBin());
      double minCh1 = hSigCh1[i]->GetBinContent(hSigCh1[i]->GetMinimumBin());
      double minYaxis = std::min(minCh0, minCh1) - 10.;
      hSigCh0[i]->GetYaxis()->SetRangeUser(minYaxis, 20);
    }
    textCh0.DrawLatex(0.5, 0.6, stringCh0);
    textCh1.DrawLatex(0.5, 0.55, stringCh1);
  }
  
  TCanvas *can_sigav = new TCanvas("data_sigav", "data_sigav", 1800, 800);
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
    hSigAvCh0[i]->SetStats(false);
    hSigAvCh0[i]->Draw();
    if(polarity==1){
      double maxYaxis = hSigAvCh0[i]->GetBinContent(hSigAvCh0[i]->GetMaximumBin());
      maxYaxis = maxYaxis+0.2*maxYaxis;
      hSigAvCh0[i]->GetYaxis()->SetRangeUser(-10, maxYaxis);
    }
    else if(polarity==-1){
      double minYaxis = hSigAvCh0[i]->GetBinContent(hSigAvCh0[i]->GetMinimumBin());
      minYaxis = minYaxis+0.2*minYaxis;
      hSigAvCh0[i]->GetYaxis()->SetRangeUser(minYaxis, 100);
    }
    textCh0.DrawLatex(0.15, 0.20, stringCh0);
    
    can_sigav->cd(i+1+nsigav);
    gPad->SetGrid(1,1);
    stringCh1 = hSigAvCh1[i]->GetTitle();
    hSigAvCh1[i]->SetTitle(" ");
    hSigAvCh1[i]->GetXaxis()->SetTitle("time [ns]");
    hSigAvCh1[i]->GetYaxis()->SetTitle("amplitude [mV]");
    if(polarity==1){
      double maxYaxis = hSigAvCh1[i]->GetBinContent(hSigAvCh1[i]->GetMaximumBin());
      maxYaxis = maxYaxis+0.2*maxYaxis;
      hSigAvCh1[i]->GetYaxis()->SetRangeUser(-10, maxYaxis);
    }
    else if(polarity==-1){
      double minYaxis = hSigAvCh1[i]->GetBinContent(hSigAvCh1[i]->GetMinimumBin());
      minYaxis = minYaxis+0.2*minYaxis;
      hSigAvCh1[i]->GetYaxis()->SetRangeUser(minYaxis, 100);
    }
    hSigAvCh1[i]->SetStats(false);
    hSigAvCh1[i]->Draw();
    textCh1.DrawLatex(0.15, 0.20, stringCh1);
  }
  
  //----- saving
  TString fname = Form("data_series%i.root", seriesNo);
  TString fname_full = outdir + fname;
  TString dbname_full = outdir + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");

  if(!file->IsOpen()){
    std::cerr << "##### Error in data.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  std::cout << "----- Saving results in ROOT file: " << fname_full << std::endl;
  
  can_ampl->Write();
  can_charge->Write();
  can_t0->Write();
  can_tot->Write();
  can_ampl_corr->Write();
  can_charge_corr->Write();
  can_t0_corr->Write();
  can_amp_pe_ch0->Write();
  can_amp_pe_ch1->Write();
  can_sig->Write();
  can_sigav->Write();
  if(collimator.Contains("Electronic")){
    can_ref->Write();
    can_ref_ch0->Write();
    can_ref_ch1->Write();
  }
  file->Close();
  
  //----- writing results to the data base
  TString table = "DATA";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES (%i, '%s')", table.Data(), seriesNo, fname_full.Data());
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  
  return 0;
}
