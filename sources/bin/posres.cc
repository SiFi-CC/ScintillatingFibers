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
  
  posres->AnalyzePositionRes();
  posres->ReconstructPos();
  TGraphErrors *gPosVsMean = posres->GetMeanGraph();
  TGraphErrors *gPosVsSigma = posres->GetPositionResGraph();
  TGraphErrors *gRecoPos = posres->GetRecoPosGraph();
  std::vector <TH1D*> ratios = posres->GetRatios();
  std::vector <TH1D*> specCh0 = posres->GetSpectra(0);
  std::vector <TH1D*> specCh1 = posres->GetSpectra(1);
  std::vector <double> results = posres->GetPositionResSeries();
  std::vector <double> posResMeas = posres->GetPositionRes();
  std::vector <double> posResMeasErr = posres->GetPositionResError();
  
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
  
  TCanvas *can_rat = new TCanvas("can_rat", "can_rat", 1200, 1200);
  can_rat->DivideSquare(npoints);
  
  TCanvas *can_spec_ch0 = new TCanvas("can_spec_ch0", "can_spec_ch0", 1200, 1200);
  can_spec_ch0->DivideSquare(npoints);
  
  TCanvas *can_spec_ch1 = new TCanvas("can_spec_ch1", "can_spec_ch1", 1200, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    can_rat->cd(i+1);
    gPad->SetGrid(1,1);
    ratios[i]->Draw();
    text.DrawLatex(0.3, 0.8, Form("#sigma = (%.2f +/- %.2f) mm", posResMeas[i], posResMeasErr[i]));
    
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
  
  TCanvas *can_posres = new TCanvas("can_posres", "can_pores", 1000, 800);
  can_posres->Divide(2,1);
  
  can_posres->cd(1);
  gPad->SetGrid(1,1);
  gPosVsMean->Draw("AP");
  
  can_posres->cd(2);
  gPad->SetGrid(1,1);
  gPosVsSigma->Draw("AP");
  text.DrawLatex(0.3, 0.8, Form("PR = (%.2f +/- %.2f) mm", results[0], results[1]));
  
  TCanvas *can_reco = new TCanvas("pos_reco", "pos_reco", 800, 800);
  gPad->SetGrid(1,1);
  gRecoPos->Draw("AP");
  
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
  
  can_rat->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  can_posres->Write();
  can_reco->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "POSITION_RESOLUTION";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES, POSITION_RES_ERR) VALUES(%i, '%s', %f, %f)", table.Data(), seriesNo, fname_full.Data(), results[0], results[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete posres;
 
  return 1; 
}
