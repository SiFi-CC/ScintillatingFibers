// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             energyres.cc              *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyRes.hh"
#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv){

  if(argc<2 || argc>6){
    std::cout << "to run type: ./energyres seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 0;
  }

  int seriesNo = atoi(argv[1]);

  SFData *data; 
  try {
    data = new SFData(seriesNo);
 }
   catch(const char* message){
    std::cerr << message << std::endl;
    std::cout << "##### Exception in energyres.cc!" << std::endl;
    return 0;
  }

  data->Print();
  int npoints = data->GetNpoints();
  TString collimator = data->GetCollimator();
  std::vector <double> positions = data->GetPositions();

  SFEnergyRes* enres;
  try{
    enres= new SFEnergyRes(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in energyres.cc!" << std::endl;
    return 0;
  }
    
  enres->CalculateEnergyRes(0);
  enres->CalculateEnergyRes(1);
  enres->CalculateEnergyRes();
    
  TGraphErrors* gEnResSum = enres->GetEnergyResolutionGraph(); 
  TGraphErrors* gEnResCh0 = enres->GetEnergyResolutionGraph(0); 
  TGraphErrors* gEnResCh1 = enres->GetEnergyResolutionGraph(1); 
  std::vector <double> enResSum = enres->GetEnergyResolution(); 
  std::vector <double> enResCh0 = enres->GetEnergyResolution(0); 
  std::vector <double> enResCh1 = enres->GetEnergyResolution(1); 
  std::vector <TH1D*> specSum   = enres->GetSpectraSum();
  std::vector <TH1D*> specCh0   = enres->GetSpectra(0);
  std::vector <TH1D*> specCh1   = enres->GetSpectra(1);
  std::vector <TH1D*> specCorrCh0 = enres->GetSpectraCorrected(0);
  std::vector <TH1D*> specCorrCh1 = enres->GetSpectraCorrected(1);
  std::vector <TH1D*> peaksSum;
  std::vector <TH1D*> peaksCh0;
  std::vector <TH1D*> peaksCh1;
  
  if(collimator=="Lead"){
    peaksSum = enres->GetPeaks();
    peaksCh0 = enres->GetPeaks(0);
    peaksCh1 = enres->GetPeaks(1);
  }

  //----- drawing channels
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
    
  //----- summed channels 
  TCanvas *can_sum = new TCanvas("can_sum", "can_sum", 700, 500);
  can_sum->cd();
  gPad->SetGrid(1,1);
  gEnResSum->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f)", enResSum[0], enResSum[1]));
    
  //----- Ch0 
  TCanvas *can_ch0 = new TCanvas("can_ch0", "can_ch0", 700, 500);
  can_ch0->cd();
  gPad->SetGrid(1,1);
  gEnResCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) ", enResCh0[0], enResCh0[1]));

  //----- Ch1 
  TCanvas *can_ch1 = new TCanvas("can_ch1", "can_ch1", 700, 500);
  can_ch1->cd();
  gPad->SetGrid(1,1);
  gEnResCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) ", enResCh1[0], enResCh1[1]));

  //----- Spectra 
  TCanvas *can_spec_sum = new TCanvas("can_spec_sum", "can_spec_sum" ,1200, 1200);
  can_spec_sum->DivideSquare(npoints);
  
  TCanvas *can_spec_corr_ch0 = new TCanvas("can_spec_corr_ch0", "can_spec_corr_ch0", 1200, 1200);
  can_spec_corr_ch0->DivideSquare(npoints);
  
  TCanvas *can_spec_corr_ch1 = new TCanvas("can_spec_corr_ch1", "can_spec_corr_ch1", 1200, 1200);
  can_spec_corr_ch1->DivideSquare(npoints);
  
  TCanvas *can_spec_ch0 = new TCanvas("can_spec_ch0", "can_spec_ch0", 1200, 1200);
  can_spec_ch0->DivideSquare(npoints);

  TCanvas *can_spec_ch1 = new TCanvas("can_spec_ch1", "can_spec_ch1", 1200, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  
  for(int i=0; i<npoints; i++){
    
    can_spec_sum->cd(i+1);
    gPad->SetGrid(1,1);
    specSum[i]->SetStats(false);
    specSum[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specSum[i]->GetYaxis()->SetTitle("counts");
    specSum[i]->SetTitle(Form("Summed PE spectrum, position %.2f mm", positions[i]));
    specSum[i]->Draw();
    if(collimator=="Lead"){
      peaksSum[i]->SetLineColor(kGreen+3);
      peaksSum[i]->Draw("same");    
    }

    can_spec_corr_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCorrCh0[i]->SetStats(false); 
    specCorrCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCorrCh0[i]->GetYaxis()->SetTitle("counts");
    specCorrCh0[i]->SetTitle(Form("Corrected Ch0 PE spectrum, position %.2f mm", positions[i]));
    specCorrCh0[i]->Draw();

    can_spec_corr_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCorrCh1[i]->SetStats(false);
    specCorrCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCorrCh1[i]->GetYaxis()->SetTitle("counts");
    specCorrCh1[i]->SetTitle(Form("Corrected Ch1 PE spectrum, position %.2f mm", positions[i]));
    specCorrCh1[i]->Draw();
    
    can_spec_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->SetStats(false);
    specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh0[i]->GetYaxis()->SetTitle("counts");
    specCh0[i]->SetTitle(Form("Ch0 PE spectrum, position %.2f mm", positions[i]));
    specCh0[i]->Draw();
    if(collimator=="Lead"){
      peaksCh0[i]->SetLineColor(kGreen+3);
      peaksCh0[i]->Draw("same");
    }
    
    can_spec_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCh1[i]->SetStats(false);
    specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh1[i]->GetYaxis()->SetTitle("counts");
    specCh1[i]->SetTitle(Form("Ch1 PE spectrum, position %.2f mm", positions[i]));
    specCh1[i]->Draw();
    if(collimator=="Lead"){
      peaksCh1[i]->SetLineColor(kGreen+3);
      peaksCh1[i]->Draw("same");
    }
  }
  
  //----- saving
  TString path = std::string(getenv("SFPATH"));
  TString fname = Form("enres_series%i.root", seriesNo);
  
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
      std::cerr << "##### Error in energyres.cc! Unable to create new direcotry!" << std::endl;
      return 0;
    }
  }
  
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in energyres.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 0;
  }

  can_sum->Write();
  can_ch0->Write();
  can_ch1->Write();
  can_spec_sum->Write();
  can_spec_corr_ch0->Write();
  can_spec_corr_ch1->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  file->Close();

  //----- writing results to the data base
  TString table = "ENERGY_RESOLUTION";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ENRES_SUM, ENRES_SUM_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), enResSum[0], enResSum[1], enResCh0[0], enResCh0[1], enResCh1[0], enResCh1[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

  delete data;
  delete enres;
  
  return 1;
} 
