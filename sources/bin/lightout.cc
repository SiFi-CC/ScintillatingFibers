// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             lightout.cc               *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "common_options.h"
#include "SFData.hh"
#include "SFLightOutput.hh"
#include "SFTools.hh"

#include <TSystem.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <DistributionContext.h>

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
    std::cout << "to run type: ./lightout seriesNo";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }
  
  SFData *data; 
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in lightout.cc!" << std::endl;
    return 1;
  }

  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in lightout.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cerr << "Description: " << desc << std::endl;
    return 1;    
  }
  
  int npoints = data->GetNpoints();
  int anaGroup = data->GetAnalysisGroup();
  TString collimator = data->GetCollimator();
  std::vector <double> positions = data->GetPositions();
  data->Print();

  SFLightOutput* lout;
  try{
    lout = new SFLightOutput(seriesNo);
  }
  catch(const char *message){
    std:: cerr << message << std::endl;
    std:: cerr << "##### Exception in lightout.cc!" << std::endl;
    return 1;
  }
  
  DistributionContext ctx;
  ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

  ctx.dim = DIM1;
  ctx.x.min = 0;
  ctx.x.max = 100;
  ctx.y.min = 0;
  ctx.y.max = 1000;

  TCanvas *can = lout->GetInputData();
  
  //----- 
  lout->CalculateLightOut(0);
  lout->CalculateLightOut(1);
  lout->CalculateLightOut();
  
  TGraphErrors *gLightOutCh0 = lout->GetLightOutputGraph(0);
  TGraphErrors *gLightOutCh1 = lout->GetLightOutputGraph(1);
  TGraphErrors *gLightOut    = lout->GetLightOutputGraph();
  LightResults LOresults     = lout->GetLOResults();
  std::vector <TH1D*>  specCh0 = lout->GetSpectra(0);
  std::vector <TH1D*>  specCh1 = lout->GetSpectra(1);
  
  //-----
  lout->CalculateLightCol(0);
  lout->CalculateLightCol(1);
  lout->CalculateLightCol();
  
  TGraphErrors *gLightColCh0 = lout->GetLightColGraph(0);
  TGraphErrors *gLightColCh1 = lout->GetLightColGraph(1);
  TGraphErrors *gLightCol    = lout->GetLightColGraph();
  LightResults LCresults     = lout->GetLCResults();
  
  //----- drawing
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
  
  //----- light output channels 0 and 1
  TCanvas *can_lout_ch = new TCanvas("lo_lout_ch", "lo_lout_ch", 1200, 600);
  can_lout_ch->Divide(2,1);
  
  can_lout_ch->cd(1);
  gPad->SetGrid(1,1);
  gLightOutCh0->GetYaxis()->SetMaxDigits(2);
  gLightOutCh0->GetXaxis()->SetLabelSize(0.035);
  gLightOutCh0->GetYaxis()->SetLabelSize(0.035);
  gLightOutCh0->GetYaxis()->SetTitleOffset(1.4);
  gLightOutCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LOresults.fResCh0, LOresults.fResCh0Err));
  
  can_lout_ch->cd(2);
  gPad->SetGrid(1,1);
  gLightOutCh1->GetYaxis()->SetMaxDigits(2);
  gLightOutCh1->GetXaxis()->SetLabelSize(0.035);
  gLightOutCh1->GetYaxis()->SetLabelSize(0.035);
  gLightOutCh1->GetYaxis()->SetTitleOffset(1.4);
  gLightOutCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LOresults.fResCh1, LOresults.fResCh1Err));
  
  //----- light output summed
  TCanvas *can_lout = new TCanvas("lo_lout","lo_lout", 700, 500);
  can_lout->cd();
  gPad->SetGrid(1,1);
  gLightOut->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LOresults.fRes, LOresults.fResErr));

   //----- light collection channels 0 and 1
  TCanvas *can_lcol_ch = new TCanvas("lo_lcol_ch", "lo_lcol_ch", 1200, 600);
  can_lcol_ch->Divide(2,1);
  
  can_lcol_ch->cd(1);
  gPad->SetGrid(1,1);
  gLightColCh0->GetYaxis()->SetMaxDigits(2);
  gLightColCh0->GetXaxis()->SetLabelSize(0.035);
  gLightColCh0->GetYaxis()->SetLabelSize(0.035);
  gLightColCh0->GetYaxis()->SetTitleOffset(1.4);
  gLightColCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LCresults.fResCh0, LCresults.fResCh0Err));
  
  can_lcol_ch->cd(2);
  gPad->SetGrid(1,1);
  gLightColCh1->GetYaxis()->SetMaxDigits(2);
  gLightColCh1->GetXaxis()->SetLabelSize(0.035);
  gLightColCh1->GetYaxis()->SetLabelSize(0.035);
  gLightColCh1->GetYaxis()->SetTitleOffset(1.4);
  gLightColCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LCresults.fResCh1, LCresults.fResCh1Err));
  
  //----- light output summed
  TCanvas *can_lcol = new TCanvas("lo_lcol","lo_lcol", 700, 500);
  can_lcol->cd();
  gPad->SetGrid(1,1);
  gLightCol->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", LCresults.fRes, LCresults.fResErr));
  
  //----- spectra
  TCanvas *can_spec_ch0 = new TCanvas("lo_spec_ch0", "lo_spec_ch0", 2000, 1200);
  can_spec_ch0->DivideSquare(npoints);
  
  TCanvas *can_spec_ch1 = new TCanvas("lo_spec_ch1", "lo_spec_ch1", 2000, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    can_spec_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->SetStats(false);
    specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh0[i]->GetYaxis()->SetTitle("counts");
    specCh0[i]->SetTitle(Form("Cherge spectrum Ch0, position %.2f mm", positions[i]));
    ctx.configureFromJson("hSpec");
    ctx.print();
    specCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    specCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    specCh0[i]->Draw();
    
    can_spec_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCh1[i]->SetStats(false);
    specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh1[i]->GetYaxis()->SetTitle("counts");
    specCh1[i]->SetTitle(Form("Cherge spectrum Ch1, position %.2f mm", positions[i]));
    specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    specCh1[i]->Draw();
  }
  
  //----- saving
  TString fname = Form("lightout_series%i.root", seriesNo);
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in lightout.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  can_lout_ch->Write();
  can_lout->Write();
  can_lcol_ch->Write();
  can_lcol->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  can->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "LIGHT_OUTPUT";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, LOUT, LOUT_ERR, LOUT_CH0, LOUT_CH0_ERR, LOUT_CH1, LOUT_CH1_ERR, LCOL, LCOL_ERR, LCOL_CH0, LCOL_CH0_ERR, LCOL_CH1, LCOL_CH1_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), LOresults.fRes, LOresults.fResErr, LOresults.fResCh0, LOresults.fResCh0Err, LOresults.fResCh1, LOresults.fResCh1Err, LCresults.fRes, LCresults.fResErr, LCresults.fResCh0, LCresults.fResCh0Err, LCresults.fResCh1, LCresults.fResCh1Err);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete lout;

  return 0;
} 
