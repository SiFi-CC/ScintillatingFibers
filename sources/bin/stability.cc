// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             stability.cc              *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFStabilityMon.hh"
#include "SFTools.hh"
#include "TCanvas.h"
#include "TLatex.h"
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"
#include <sys/stat.h> 
#include <sys/types.h>
#include "common_options.h"

int main(int argc, char **argv){
  
  if(argc<2 || argc>6){
    std::cout << "to run type: ./stability seriesNo";
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
   std::cerr << "##### Exception in stability.cc!" << std::endl;
   return 1;
  }
 
  data->Print();
  int npoints = data->GetNpoints();
  std::vector <double> positions = data->GetPositions();
    
  TString desc = data->GetDescription();
  if(!desc.Contains("Stability monitoring")){
    std::cerr << "##### Error in stability.cc! This is not stability monitoring series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cerr << "Description: " << desc << std::endl;
    return 1;
  }
  
  SFStabilityMon *stab;
  try{
    stab = new SFStabilityMon(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Error in stability.cc" << std::endl;
    return 1;
  }
  
  //----- stability of channel 0
  stab->AnalyzeStability(0);
  TGraphErrors *gCh0PeakPos  = stab->GetPeakPosGraph(0);
  TGraphErrors *gCh0Residual = stab->GetResidualsGraph(0);
  double stdDevCh0 = stab->GetStdDev(0);
  double meanCh0   = stab->GetMean(0);
  std::vector <TH1D*> specCh0 = stab->GetSpectra(0);
  
  //----- stability of channel 1
  stab->AnalyzeStability(1);
  TGraphErrors *gCh1PeakPos  = stab->GetPeakPosGraph(1);
  TGraphErrors *gCh1Residual = stab->GetResidualsGraph(1);
  double stdDevCh1 = stab->GetStdDev(1);
  double meanCh1   = stab->GetMean(1);
  std::vector <TH1D*> specCh1 = stab->GetSpectra(1);
  
  TCanvas *can = new TCanvas("can", "can", 1000, 700);
  TPad *pad_peakPos = new TPad("pad_peakPos", "pad_peakPos", 0, 0.3, 1, 1, 10, 0);
  TPad *pad_res     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);
  
  can->cd(0);
  pad_peakPos->Draw();
  pad_peakPos->cd();
  pad_peakPos->SetGrid(1,1);
  gCh0PeakPos->Draw("AP");
  gCh0PeakPos->SetTitle("511 keV peak position stability");
  gCh0PeakPos->SetMarkerColor(kBlack);
  gCh0PeakPos->SetLineColor(kBlack);
  gCh0PeakPos->GetFunction("funPol0")->SetLineColor(kGray+1);
  gCh1PeakPos->Draw("P");
  gCh1PeakPos->SetMarkerColor(kRed);
  gCh1PeakPos->SetLineColor(kRed);
  gCh1PeakPos->GetFunction("funPol0")->SetLineColor(kRed-7);
  
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
  text.SetTextColor(kBlack);
  text.DrawLatex(0.2, 0.3, Form("#bar{PP}_{ch0} = %.2f +/- %.2f", meanCh0, stdDevCh0));
  text.SetTextColor(kRed);
  text.DrawLatex(0.2, 0.2, Form("#bar{PP}_{ch1} = %.2f +/- %.2f", meanCh1, stdDevCh1));
  
  can->cd(0);
  pad_res->Draw();
  pad_res->cd();
  pad_res->SetGrid(1,1);
  gCh0Residual->Draw("AP");
  gCh0Residual->SetTitle("Residuals graph");
  gCh0Residual->SetMarkerColor(kBlack);
  gCh0Residual->SetLineColor(kBlack);
  gCh1Residual->Draw("P");
  gCh1Residual->SetMarkerColor(kRed);
  gCh1Residual->SetLineColor(kRed);
  
  TCanvas *can_ch0 = new TCanvas("can_ch0", "can_ch0", 1200, 1000);
  can_ch0->DivideSquare(npoints);
  
  TCanvas *can_ch1 = new TCanvas("can_ch1", "can_ch1", 1200, 1000);
  can_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    can_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->SetTitle(Form("PE spectrum, ch0, position %.2f", positions[i]));
    specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh0[i]->GetYaxis()->SetTitle("counts");
    specCh0[i]->Draw();
    
    can_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCh1[i]->SetTitle(Form("PE spectrum, ch1, position %.2f", positions[i]));
    specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh1[i]->GetYaxis()->SetTitle("counts");
    specCh1[i]->Draw();
  }
  
  //----- saving
  TString fname = Form("stability_series%i.root", seriesNo);
  TString outdir;
  TString dbase;

  int ret = parse_common_options(argc, argv, outdir, dbase);
  if(ret != 0) 
    exit(ret);

  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in stability.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  can->Write();
  can_ch0->Write();
  can_ch1->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "STABILITY_MON";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, CH0_MEAN, CH0_STDDEV, CH1_MEAN, CH1_STDDEV) VALUES (%i, '%s', %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), meanCh0, stdDevCh0, meanCh1, stdDevCh1);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete stab;
  
  return 0;
}
