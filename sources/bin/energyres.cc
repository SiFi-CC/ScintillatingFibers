// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             energyres.cc              *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "common_options.h"
#include "SFEnergyRes.hh"

#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>

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
    std::cout << "to run type: ./energyres seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }

  //----- accessing results of energy resolution analysis
  SFData *data; 
  try {
    data = new SFData(seriesNo);
 }
   catch(const char* message){
    std::cerr << message << std::endl;
    std::cout << "##### Exception in energyres.cc!" << std::endl;
    return 1;
  }

  data->Print();
  int npoints = data->GetNpoints();
  int anaGroup = data->GetAnalysisGroup();
  TString collimator = data->GetCollimator();
  std::vector <double> positions = data->GetPositions();
  
  SFEnergyRes* enres;
  try{
    enres= new SFEnergyRes(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in energyres.cc!" << std::endl;
    return 1;
  }
    
  enres->CalculateEnergyRes(0);
  enres->CalculateEnergyRes(1);
  enres->CalculateEnergyRes();
    
  TGraphErrors* gEnResAve = enres->GetEnergyResolutionGraph(); 
  TGraphErrors* gEnResCh0 = enres->GetEnergyResolutionGraph(0); 
  TGraphErrors* gEnResCh1 = enres->GetEnergyResolutionGraph(1); 

  std::vector <TH1D*> specAve   = enres->GetSpectra();
  std::vector <TH1D*> specCh0   = enres->GetSpectra(0);
  std::vector <TH1D*> specCh1   = enres->GetSpectra(1);

  EnergyResResults results = enres->GetResults();
  
  //----- accessing json file
  DistributionContext ctx;
  ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

  ctx.dim = DIM1;
  ctx.x.min = 0;
  ctx.x.max = 100;
  ctx.y.min = 0;
  ctx.y.max = 1000;
  
  //----- drawing energy resolution graphs
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
    
  //----- averaged channels 
  TCanvas *can_ave = new TCanvas("er_ave", "er_ave", 700, 500);
  can_ave->cd();
  gPad->SetGrid(1,1);
  gEnResAve->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                 results.fEnergyResAve, results.fEnergyResAveErr));
    
  //----- channel 0 
  TCanvas *can_ch0 = new TCanvas("er_ch0", "er_ch0", 700, 500);
  can_ch0->cd();
  gPad->SetGrid(1,1);
  gEnResCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                 results.fEnergyResCh0, results.fEnergyResCh0Err));

  //----- channel 1 
  TCanvas *can_ch1 = new TCanvas("er_ch1", "er_ch1", 700, 500);
  can_ch1->cd();
  gPad->SetGrid(1,1);
  gEnResCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                 results.fEnergyResCh1, results.fEnergyResCh1Err));

  //----- drawing spectra 
  TCanvas *can_spec_ave = new TCanvas("er_spec_ave", "er_spec_ave", 2000, 1200);
  can_spec_ave->DivideSquare(npoints);
  
  TCanvas *can_spec_ch0 = new TCanvas("er_spec_ch0", "er_spec_ch0", 2000, 1200);
  can_spec_ch0->DivideSquare(npoints);

  TCanvas *can_spec_ch1 = new TCanvas("er_spec_ch1", "er_spec_ch1", 2000, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  
  for(int i=0; i<npoints; i++){
    
    can_spec_ave->cd(i+1);
    gPad->SetGrid(1,1);
    specAve[i]->SetStats(false);
    specAve[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specAve[i]->GetYaxis()->SetTitle("counts");
    specAve[i]->SetTitle(Form("Summed PE spectrum, position %.2f mm", positions[i]));
    ctx.configureFromJson("hSpecAv");
    ctx.print();
    specAve[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    specAve[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    specAve[i]->Draw();
    
    can_spec_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->SetStats(false);
    specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh0[i]->GetYaxis()->SetTitle("counts");
    specCh0[i]->SetTitle(Form("Ch0 PE spectrum, position %.2f mm", positions[i]));
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
    specCh1[i]->SetTitle(Form("Ch1 PE spectrum, position %.2f mm", positions[i]));
    specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
    specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
    specCh1[i]->Draw();
  }
  
  //----- saving ROOT file
  TString fname = Form("enres_series%i.root", seriesNo);
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in energyres.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }

  can_ave->Write();
  can_ch0->Write();
  can_ch1->Write();
  can_spec_ave->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  file->Close();

  //----- writing results to the data base
  TString table = "ENERGY_RESOLUTION";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ENRES_AV, ENRES_AV_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), results.fEnergyResAve, results.fEnergyResAveErr, results.fEnergyResCh0, results.fEnergyResCh0Err, results.fEnergyResCh1, results.fEnergyResCh1Err);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

  delete data;
  delete enres;
  
  return 0;
} 
