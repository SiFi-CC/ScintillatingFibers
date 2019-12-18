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
#include <sys/stat.h> 
#include <sys/types.h> 
#include "common_options.h"

int main(int argc, char **argv){
  
  if(argc<2 || argc>6){
    std::cout << "to run type: ./posres seriesNo";
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
  
  SFPositionRes *posres;
  
  try{
    posres = new SFPositionRes(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in posres.cc!" << std::endl;
    return 1;
  }
  
  posres->AnalyzePositionRes();
  TGraphErrors *gPosRecoVsPos     = posres->GetPositionRecoGraph();
  TGraphErrors *gPosResVsPos      = posres->GetPositionResGraph();
  std::vector <TH1D*> hPosReco    = posres->GetPositionRecoDist();
  std::vector <double> results    = posres->GetPositionResSeries();
  std::vector <double> posRes     = posres->GetPositionRes();
  std::vector <double> posResErr  = posres->GetPositionResError();
  std::vector <double> posReco    = posres->GetPositionReco();
  std::vector <double> posRecoErr = posres->GetPositionRecoError();
  
  std::vector <TH1D*> spec = posres->GetSpectra();
  
  std::vector <SFPeakFinder*> peakFinCh0;
  std::vector <SFPeakFinder*> peakFinCh1;
  std::vector <double> xmin(npoints);
  std::vector <double> xmax(npoints);
  
  for(int i=0; i<npoints; i++){
    peakFinCh0.push_back(new SFPeakFinder(spec[i], 0));
    peakFinCh0[i]->FindPeakRange(xmin[i], xmax[i]);
  }
  
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.05);
  
  TLine line;
  line.SetLineColor(kRed);
  line.SetLineWidth(1);
  
  TCanvas *can_posreco = new TCanvas("can_posreco", "can_posreco", 1200, 1200);
  can_posreco->DivideSquare(npoints);
  
  TCanvas *can_spec = new TCanvas("can_spec", "can_spec", 1200, 1200);
  can_spec->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
      
    can_posreco->cd(i+1);
    gPad->SetGrid(1,1);
    hPosReco[i]->Draw();
    text.DrawLatex(0.3, 0.8, Form("#mu = (%.2f +/- %.2f) mm", posReco[i], posRecoErr[i]));
    text.DrawLatex(0.3, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", posRes[i], posResErr[i]));
    
    can_spec->cd(i+1);
    gPad->SetGrid(1,1);
    spec[i]->Draw();
    line.DrawLine(xmin[i], 0, xmin[i], spec[i]->GetBinContent(spec[i]->GetMaximumBin()));
    line.DrawLine(xmax[i], 0, xmax[i], spec[i]->GetBinContent(spec[i]->GetMaximumBin()));
  }
  
  TCanvas *can_posres = new TCanvas("can_posres", "can_pores", 1000, 800);
  can_posres->Divide(2,1);
  
  can_posres->cd(1);
  gPad->SetGrid(1,1);
  gPosRecoVsPos->Draw("AP");
  
  can_posres->cd(2);
  gPad->SetGrid(1,1);
  gPosResVsPos->Draw("AP");
  text.DrawLatex(0.3, 0.8, Form("PR = (%.2f +/- %.2f) mm", results[0], results[1]));
  
  //----- saving
  TString fname = Form("posres_series%i.root", seriesNo);
  TString outdir;
  TString dbase;

  int ret = parse_common_options(argc, argv, outdir, dbase);
  if(ret != 0) 
    exit(ret);
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in posres.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  can_posreco->Write();
  can_spec->Write();
  can_posres->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "POSITION_RESOLUTION";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES, POSITION_RES_ERR) VALUES(%i, '%s', %f, %f)", table.Data(), seriesNo, fname_full.Data(), results[0], results[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete posres;
 
  return 0; 
}
