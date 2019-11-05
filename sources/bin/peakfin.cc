// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              peakfin.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFPeakFinder.hh"
#include "SFData.hh"
#include "SFAttenuation.hh"
#include "TCanvas.h"
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
    std::cerr << "##### Exception in peakfin.cc!" << std::endl;
    return 1;
  }
  
  data->Print();
  
  SFAttenuation *att;
  try{
    att = new SFAttenuation(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in peakfin.cc!" << std::endl;
    return 1;
  }
  
  int npoints = data->GetNpoints();
  std::vector <double> positions = data->GetPositions();
  std::vector <TH1D*> hSpecCh0 = data->GetSpectra(0, SFSelectionType::PE, "ch_0.fPE>0");
  std::vector <TH1D*> hSpecCh1 = data->GetSpectra(1, SFSelectionType::PE, "ch_1.fPE>0");
  std::vector <TH1D*> hSpecSum;
    
  std::vector <SFPeakFinder*> pfCh0;
  std::vector <SFPeakFinder*> pfCh1;
  std::vector <SFPeakFinder*> pfSum;
  
  std::vector <double> customNum(4);
  customNum[1] = att->GetAttLength();
  customNum[3] = att->GetAttLength();
  double distCh0, distCh1;
  TString cut = "ch_0.fPE>0 && ch_1.fPE>0 && ch_0.fT0>0 && ch_1.fT0>0";
  
  TCanvas *canCh0 = new TCanvas("canCh0", "canCh0", 1200, 1200);
  canCh0->DivideSquare(npoints);
  
  TCanvas *canCh1 = new TCanvas("canCh1", "canCh1", 1200, 1200);
  canCh1->DivideSquare(npoints);
  
  TCanvas *canSum = new TCanvas("canSum", "canSum", 1200, 1200);
  canSum->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    distCh0 = positions[i];
    distCh1 = 100, - positions[i];
    customNum[0] = -distCh0;
    customNum[2] = -distCh1;
    hSpecSum.push_back(data->GetCustomHistogram(SFSelectionType::PEAttCorrectedSum, 
                                                cut, positions[i], customNum));
    
    pfCh0.push_back(new SFPeakFinder(hSpecCh0[i], 0, 1));   // verbose = 1, tests = 1
    pfCh1.push_back(new SFPeakFinder(hSpecCh1[i], 0, 1));
    pfSum.push_back(new SFPeakFinder(hSpecSum[i], 0, 1));
    
    //pfCh0[i]->
    //pfCh1[i]->
    //pfSum[i]->
    
    canCh0->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh0[i]->SetStats(0);
    hSpecCh0[i]->GetXaxis()->SetTitle("energy [P.E.]");
    hSpecCh0[i]->GetYaxis()->SetTitle("counts");
    hSpecCh0[i]->SetTitle(Form("PE spectrum: S%i ch0 %.1f mm", seriesNo, positions[i]));
    hSpecCh0[i]->Draw();
    
    canCh1->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh1[i]->SetStats(0);
    hSpecCh1[i]->GetXaxis()->SetTitle("energy [P.E.]");
    hSpecCh1[i]->GetYaxis()->SetTitle("counts");
    hSpecCh1[i]->SetTitle(Form("PE spectrum: S%i ch1 %.1f mm", seriesNo, positions[i]));
    hSpecCh1[i]->Draw();
    
    canSum->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecSum[i]->SetStats(0);
    hSpecSum[i]->GetXaxis()->SetTitle("energy [P.E.]");
    hSpecSum[i]->GetYaxis()->SetTitle("counts");
    hSpecSum[i]->SetTitle(Form("Summed and corrected PE spectrum: S%i %.1f mm", 
                               seriesNo, positions[i]));
    hSpecSum[i]->Draw();
  }
  
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
    std::cerr << "##### Error in peakfin.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  canCh0->Write();
  canCh1->Write();
  canSum->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "PEAK_FINDER";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES(%i, '%s')", table.Data(), seriesNo, fname_full.Data());
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete att;
  
  return 0;
}
