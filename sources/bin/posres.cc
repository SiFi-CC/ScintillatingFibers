// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               posres.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFPositionRes.hh"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv){
  
  if(argc<2 || argc>6){
    std::cout << "to run type: ./posres seriesNo";
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
  
  SFPositionRes *posres;
  
  try{
    posres = new SFPositionRes(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in posres.cc!" << std::endl;
    return 0;
  }
  
  posres->AnalyzePositionResMLR();
  TGraphErrors *gMLR_MeanVsPos   = posres->GetMeanGraph("MLR");
  TGraphErrors *gMLR_PosResVsPos = posres->GetPositionResGraph("MLR");
  std::vector <TH1D*> MLR_ratios = posres->GetRatios("MLR");
  std::vector <double> MLR_results   = posres->GetPositionResSeries("MLR");
  std::vector <double> MLR_posRes    = posres->GetPositionRes("MLR");
  std::vector <double> MLR_posResErr = posres->GetPositionResError("MLR");
  
  posres->AnalyzePositionResY();
  TGraphErrors *gY_MeanVsPos   = posres->GetMeanGraph("Y");
  TGraphErrors *gY_PosResVsPos = posres->GetPositionResGraph("Y");
  std::vector <TH1D*> Y_ratios = posres->GetRatios("Y");
  std::vector <double> Y_results   = posres->GetPositionResSeries("Y");
  std::vector <double> Y_posRes    = posres->GetPositionRes("Y");
  std::vector <double> Y_posResErr = posres->GetPositionResError("Y");
  
  std::vector <TH1D*> specCh0 = posres->GetSpectra(0);
  std::vector <TH1D*> specCh1 = posres->GetSpectra(1);
  
  std::vector <SFPeakFinder*> peakFinCh0;
  std::vector <SFPeakFinder*> peakFinCh1;
  double xminCh0, xmaxCh0;
  double xminCh1, xmaxCh1;
  
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
  
  TLine line;
  line.SetLineColor(kRed);
  line.SetLineWidth(1);
  
  TCanvas *can_ratMLR = new TCanvas("can_ratMLR", "can_ratMLR", 1200, 1200);
  can_ratMLR->DivideSquare(npoints);
  
  TCanvas *can_ratY = new TCanvas("can_ratY", "can_ratY", 1200, 1200);
  can_ratY->DivideSquare(npoints);
  
  TCanvas *can_spec_ch0 = new TCanvas("can_spec_ch0", "can_spec_ch0", 1200, 1200);
  can_spec_ch0->DivideSquare(npoints);
  
  TCanvas *can_spec_ch1 = new TCanvas("can_spec_ch1", "can_spec_ch1", 1200, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    can_ratMLR->cd(i+1);
    gPad->SetGrid(1,1);
    MLR_ratios[i]->Draw();
    text.DrawLatex(0.3, 0.8, Form("#sigma = (%.2f +/- %.2f) mm", MLR_posRes[i], MLR_posResErr[i]));
    
    can_ratY->cd(i+1);
    gPad->SetGrid(1,1);
    Y_ratios[i]->Draw();
    text.DrawLatex(0.3, 0.8, Form("#sigma = (%.2f +/- %.2f) mm", Y_posRes[i], Y_posResErr[i]));
    
    peakFinCh0.push_back(new SFPeakFinder(specCh0[i], 0));
    peakFinCh0[i]->FindPeakRange(xminCh0, xmaxCh0);
    can_spec_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->Draw();
    line.DrawLine(xminCh0, 0, xminCh0, specCh0[i]->GetBinContent(specCh0[i]->GetMaximumBin()));
    line.DrawLine(xmaxCh0, 0, xmaxCh0, specCh0[i]->GetBinContent(specCh0[i]->GetMaximumBin()));
    
    peakFinCh1.push_back(new SFPeakFinder(specCh1[i], 0));
    peakFinCh1[i]->FindPeakRange(xminCh1, xmaxCh1);
    can_spec_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCh1[i]->Draw();
    line.DrawLine(xminCh1, 0, xminCh1, specCh1[i]->GetBinContent(specCh1[i]->GetMaximumBin()));
    line.DrawLine(xmaxCh1, 0, xmaxCh1, specCh1[i]->GetBinContent(specCh1[i]->GetMaximumBin()));
  }
  
  TCanvas *can_posresMLR = new TCanvas("can_posresMLR", "can_poresMLR", 1000, 800);
  can_posresMLR->Divide(2,1);
  
  can_posresMLR->cd(1);
  gPad->SetGrid(1,1);
  gMLR_MeanVsPos->Draw("AP");
  
  can_posresMLR->cd(2);
  gPad->SetGrid(1,1);
  gMLR_PosResVsPos->Draw("AP");
  text.DrawLatex(0.3, 0.8, Form("PR = (%.2f +/- %.2f) mm", MLR_results[0], MLR_results[1]));
  
  TCanvas *can_posresY = new TCanvas("pos_reco", "pos_reco", 1000, 800);
  can_posresY->Divide(2,1);
  
  can_posresY->cd(1);
  gPad->SetGrid(1,1);
  gY_MeanVsPos->Draw("AP");
  
  can_posresY->cd(2);
  gPad->SetGrid(1,1);
  gY_PosResVsPos->Draw("AP");
  
  //----- saving
  TString path = std::string(getenv("SFPATH"));
  TString fname = Form("posres_series%i.root", seriesNo);

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
      std::cerr << "##### Error in posres.cc! Unable to create new direcotry!" << std::endl;
      return 0;
    }
  }
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in posres.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 0;
  }
  
  can_ratMLR->Write();
  can_ratY->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  can_posresMLR->Write();
  can_posresY->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "POSITION_RESOLUTION";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES_MLR, POSITION_RES_MLR_ERR, POSITION_RES_Y, POSITION_RES_Y_ERR) VALUES(%i, '%s', %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), MLR_results[0], MLR_results[1], Y_results[0], Y_results[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete posres;
 
  return 1; 
}
