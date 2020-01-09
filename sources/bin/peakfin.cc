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
  
  SFAttenuation *att;
  try{
    att = new SFAttenuation(seriesNo);
  }
  catch(const char* message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in peakfin.cc!" << std::endl;
    return 1;
  }
  
  att->AttAveragedCh();
  
  int npoints = data->GetNpoints();
  std::vector <double> positions = data->GetPositions();
  std::vector <int>    measurementsIDs = data->GetMeasurementsIDs();
  std::vector <TH1D*> hSpecCh0 = data->GetSpectra(0, SFSelectionType::PE, "ch_0.fPE>0");
  std::vector <TH1D*> hSpecCh1 = data->GetSpectra(1, SFSelectionType::PE, "ch_1.fPE>0");
  std::vector <TH1D*> hSpecAv  = data->GetCustomHistograms(SFSelectionType::PEAverage, 
                                 "ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fPE>0 && ch_1.fPE>0");
  std::vector <TH1D*> hSpecSum;
  
  std::vector <TH1D*> hPeakCh0;
  std::vector <TH1D*> hPeakCh1;
  std::vector <TH1D*> hPeakAv;
  std::vector <TH1D*> hPeakSum;
    
  std::vector <SFPeakFinder*> pfCh0;
  std::vector <SFPeakFinder*> pfCh1;
  std::vector <SFPeakFinder*> pfSum;
  std::vector <SFPeakFinder*> pfAve;
  
  std::vector <double> customNum(4);
  AttenuationResults results = att->GetResults();
  customNum[1] = results.fAttCombPol1;
  customNum[3] = results.fAttCombPol1;
  double distCh0, distCh1;
  TString cut = "ch_0.fPE>0 && ch_1.fPE>0 && ch_0.fT0>0 && ch_1.fT0>0";
  
  for(int i=0; i<npoints; i++){
    distCh0 = positions[i];
    distCh1 = 100. - positions[i];
    customNum[0] = -distCh0;
    customNum[2] = -distCh1;
    hSpecSum.push_back(data->GetCustomHistogram(SFSelectionType::PEAttCorrectedSum, 
                                                cut, measurementsIDs[i], customNum));
    
    pfCh0.push_back(new SFPeakFinder(hSpecCh0[i], 1, 1));   // verbose = 1, tests = 1
    pfCh1.push_back(new SFPeakFinder(hSpecCh1[i], 1, 1));
    pfSum.push_back(new SFPeakFinder(hSpecSum[i], 1, 1));
    pfAve.push_back(new SFPeakFinder(hSpecAv[i], 1, 1));
    
    pfCh0[i]->FindPeakFit();
    pfCh1[i]->FindPeakFit();
    pfSum[i]->FindPeakFit();
    pfAve[i]->FindPeakFit();
    
    pfCh0[i]->SubtractBackground();
    pfCh1[i]->SubtractBackground();
    pfSum[i]->SubtractBackground();
    pfAve[i]->SubtractBackground();
    
    hPeakCh0.push_back(pfCh0[i]->GetPeak());
    hPeakCh1.push_back(pfCh1[i]->GetPeak());
    hPeakAv.push_back(pfAve[i]->GetPeak());
    hPeakSum.push_back(pfSum[i]->GetPeak());
  }
  
  //----- results of fitting
  TCanvas *canCh0 = new TCanvas("canCh0", "canCh0", 1200, 1200);
  canCh0->DivideSquare(npoints);
  
  TCanvas *canCh1 = new TCanvas("canCh1", "canCh1", 1200, 1200);
  canCh1->DivideSquare(npoints);
  
  TCanvas *canSum = new TCanvas("canSum", "canSum", 1200, 1200);
  canSum->DivideSquare(npoints);
  
  TCanvas *canAve = new TCanvas("canAve", "canAve", 1200, 1200);
  canAve->DivideSquare(npoints);
  
  //----- results of background subtracting
  TCanvas *canCh0_bgs = new TCanvas("canCh0_bgs", "canCh0_bgs", 1200, 1200);
  canCh0_bgs->DivideSquare(npoints);
  
  TCanvas *canCh1_bgs = new TCanvas("canCh1_bgs", "canCh1_bgs", 1200, 1200);
  canCh1_bgs->DivideSquare(npoints);
  
  TCanvas *canSum_bgs = new TCanvas("canSum_bgs", "canSum_bgs", 1200, 1200);
  canSum_bgs->DivideSquare(npoints);
  
  TCanvas *canAve_bgs = new TCanvas("canAve_bgs", "canAve_bgs", 1200, 1200);
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
    hSpecCh0[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    hSpecCh0[i]->DrawClone();
    
    canCh0_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh0[i]->Draw();
    hSpecCh0[i]->GetFunction("fun_pol1_clone")->Delete();
    hSpecCh0[i]->GetFunction("fun_expo_clone")->Delete();
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
    hSpecCh1[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    hSpecCh1[i]->DrawClone();
    
    canCh1_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecCh1[i]->Draw();
    hSpecCh1[i]->GetFunction("fun_pol1_clone")->Delete();
    hSpecCh1[i]->GetFunction("fun_expo_clone")->Delete();
    hSpecCh1[i]->GetFunction("fun_gaus_clone")->Delete();
    hPeakCh1[i]->SetLineColor(kMagenta);
    hPeakCh1[i]->Draw("same");
    
    canSum->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecSum[i]->SetStats(0);
    hSpecSum[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hSpecSum[i]->GetYaxis()->SetTitle("counts");
    hSpecSum[i]->SetTitle(Form("Summed and corrected PE spectrum: S%i %.1f mm", 
                               seriesNo, positions[i]));
    maxval = hSpecSum[i]->GetBinContent(hSpecSum[i]->GetMaximumBin());
    hSpecSum[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    hSpecSum[i]->DrawClone();
    
    canSum_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecSum[i]->Draw();
    hSpecSum[i]->GetFunction("fun_pol1_clone")->Delete();
    hSpecSum[i]->GetFunction("fun_expo_clone")->Delete();
    hSpecSum[i]->GetFunction("fun_gaus_clone")->Delete();
    hPeakSum[i]->SetLineColor(kMagenta);
    hPeakSum[i]->Draw("same");
    
    canAve->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecAv[i]->SetStats(0);
    hSpecAv[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hSpecAv[i]->GetYaxis()->SetTitle("counts");
    hSpecAv[i]->SetTitle(Form("Average charge spectrum: S%i %.1f mm", seriesNo, positions[i]));
    maxval = hSpecAv[i]->GetBinContent(hSpecAv[i]->GetMaximumBin());
    hSpecAv[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
    hSpecAv[i]->DrawClone();
    
    canAve_bgs->cd(i+1);
    gPad->SetGrid(1,1);
    hSpecAv[i]->Draw();
    hSpecAv[i]->GetFunction("fun_pol1_clone")->Delete();
    hSpecAv[i]->GetFunction("fun_expo_clone")->Delete();
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
  canSum->Write();
  canAve->Write();
  canCh0_bgs->Write();
  canCh1_bgs->Write();
  canSum_bgs->Write();
  canAve_bgs->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "PEAK_FINDER";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES(%i, '%s')", table.Data(), seriesNo, fname_full.Data());
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete att;
  
  return 0;
}
