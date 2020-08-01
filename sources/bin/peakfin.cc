// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              peakfin.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "common_options.h"
#include "SFPeakFinder.hh"
#include "SFData.hh"
#include "SFAttenuation.hh"

#include <DistributionContext.h>

#include <TCanvas.h>

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
    std::cout << "to run type: ./posres seriesNo";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }
  
  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in peakfin.cc!" << std::endl;
    return 1;
  }
  
  data->Print();
  
  double s = SFTools::GetSigmaBL(data->GetSiPM());
  std::vector <double> sigmas = {s, s};
  TString cutCh0    = SFDrawCommands::GetCut(SFCutType::SpecCh0, sigmas);
  TString cutCh1    = SFDrawCommands::GetCut(SFCutType::SpecCh1, sigmas);
  TString cutCh0Ch1 = SFDrawCommands::GetCut(SFCutType::CombCh0Ch1, sigmas);
  
  int npoints = data->GetNpoints();
  int anaGroup = data->GetAnalysisGroup();
  std::vector <double> positions = data->GetPositions();
  std::vector <int>    measurementsIDs = data->GetMeasurementsIDs();
  std::vector <TH1D*> hSpecCh0 = data->GetSpectra(0, SFSelectionType::PE, cutCh0);
  std::vector <TH1D*> hSpecCh1 = data->GetSpectra(1, SFSelectionType::PE, cutCh1);
  std::vector <TH1D*> hSpecAv  = data->GetCustomHistograms(SFSelectionType::PEAverage,cutCh0Ch1);
  
  std::vector <TH1D*> hPeakCh0;
  std::vector <TH1D*> hPeakCh1;
  std::vector <TH1D*> hPeakAv;
    
  std::vector <SFPeakFinder*> pfCh0;
  std::vector <SFPeakFinder*> pfCh1;
  std::vector <SFPeakFinder*> pfAve;
   
  
  DistributionContext ctx;
  ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

  ctx.dim = DIM1;
  ctx.x.min = 0;
  ctx.x.max = 100;
  ctx.y.min = 0;
  ctx.y.max = 1000;
  
  for(int i=0; i<npoints; i++){
    
    pfCh0.push_back(new SFPeakFinder(hSpecCh0[i], 0, 1));   // verbose = 1, tests = 1
    pfCh1.push_back(new SFPeakFinder(hSpecCh1[i], 0, 1));
    pfAve.push_back(new SFPeakFinder(hSpecAv[i], 0, 1));
    
    pfCh0[i]->FindPeakFit();
    pfCh1[i]->FindPeakFit();
    pfAve[i]->FindPeakFit();
    
    pfCh0[i]->SubtractBackground();
    pfCh1[i]->SubtractBackground();
    pfAve[i]->SubtractBackground();
    
    hPeakCh0.push_back(pfCh0[i]->GetPeak());
    hPeakCh1.push_back(pfCh1[i]->GetPeak());
    hPeakAv.push_back(pfAve[i]->GetPeak());
  }
  
  //----- results of fitting
  TCanvas *canCh0 = new TCanvas("pf_ch0", "pf_ch0", 2000, 1200);
  canCh0->DivideSquare(npoints);
  
  TCanvas *canCh1 = new TCanvas("pf_ch1", "pf_ch1", 2000, 1200);
  canCh1->DivideSquare(npoints);
  
  TCanvas *canAve = new TCanvas("pf_ave", "pf_ave", 2000, 1200);
  canAve->DivideSquare(npoints);
  
  //----- results of background subtracting
  TCanvas *canCh0_bgs = new TCanvas("pf_ch0_bgs", "pf_ch0_bgs", 2000, 1200);
  canCh0_bgs->DivideSquare(npoints);
  
  TCanvas *canCh1_bgs = new TCanvas("pf_ch1_bgs", "pf_ch1_bgs", 2000, 1200);
  canCh1_bgs->DivideSquare(npoints);
  
  TCanvas *canAve_bgs = new TCanvas("pf_ave_bgs", "pf_ave_bgs", 2000, 1200);
  canAve_bgs->DivideSquare(npoints);
  
  //----- drawing
  
  double maxval = 0;
  
  for(int i=0; i<npoints; i++){
    canCh0->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh0[i]->SetStats(0);
    hSpecCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hSpecCh0[i]->GetYaxis()->SetTitle("counts");
    hSpecCh0[i]->SetTitle(Form("PE spectrum: S%i ch0 %.1f mm", seriesNo, positions[i]));
    maxval = hSpecCh0[i]->GetBinContent(hSpecCh0[i]->GetMaximumBin());
    //hSpecCh0[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    ctx.configureFromJson("hSpec");
    ctx.print();
    hSpecCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hSpecCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hSpecCh0[i]->DrawClone();
    
    canCh0_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh0[i]->Draw();
    hSpecCh0[i]->GetFunction("fun_bg_clone")->Delete();
    hSpecCh0[i]->GetFunction("fun_gaus_clone")->Delete();
    hPeakCh0[i]->SetLineColor(kMagenta);
    hPeakCh0[i]->Draw("same");
    
    canCh1->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh1[i]->SetStats(0);
    hSpecCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hSpecCh1[i]->GetYaxis()->SetTitle("counts");
    hSpecCh1[i]->SetTitle(Form("PE spectrum: S%i ch1 %.1f mm", seriesNo, positions[i]));
    maxval = hSpecCh1[i]->GetBinContent(hSpecCh1[i]->GetMaximumBin());
    //hSpecCh1[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    hSpecCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hSpecCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hSpecCh1[i]->DrawClone();
    
    canCh1_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh1[i]->Draw();
    hSpecCh1[i]->GetFunction("fun_bg_clone")->Delete();
    hSpecCh1[i]->GetFunction("fun_gaus_clone")->Delete();
    hPeakCh1[i]->SetLineColor(kMagenta);
    hPeakCh1[i]->Draw("same");
    
    canAve->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecAv[i]->SetStats(0);
    hSpecAv[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hSpecAv[i]->GetYaxis()->SetTitle("counts");
    hSpecAv[i]->SetTitle(Form("Average charge spectrum: S%i %.1f mm", seriesNo, positions[i]));
    maxval = hSpecAv[i]->GetBinContent(hSpecAv[i]->GetMaximumBin());
    //hSpecAv[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    ctx.configureFromJson("hSpecAv");
    ctx.print();
    hSpecAv[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    hSpecAv[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    hSpecAv[i]->DrawClone();
    
    canAve_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecAv[i]->Draw();
    hSpecAv[i]->GetFunction("fun_bg_clone")->Delete();
    hSpecAv[i]->GetFunction("fun_gaus_clone")->Delete();
    hPeakAv[i]->SetLineColor(kMagenta);
    hPeakAv[i]->Draw("same");
  }
  
  //----- saving
  TString fname = Form("peakfin_series%i.root", seriesNo);
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in peakfin.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  canCh0->Write();
  canCh1->Write();
  canAve->Write();
  canCh0_bgs->Write();
  canCh1_bgs->Write();
  canAve_bgs->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "PEAK_FINDER";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES(%i, '%s')", table.Data(), seriesNo, fname_full.Data());
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  //delete att;
  
  return 0;
}
