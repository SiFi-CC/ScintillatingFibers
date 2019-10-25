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

#include "SFData.hh"
#include "SFLightOutput.hh"
#include "SFTools.hh"
#include "CmdLineConfig.hh"
#include "CmdLineOption.hh"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv){
    
  if(argc<2 || argc>6){
    std::cout << "to run type: ./lightout seriesNo";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 0;
  }
 
  int seriesNo = atoi(argv[1]);
  
  SFData *data; 
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in lightout.cc!" << std::endl;
    return 0;
  }

  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    std::cerr << "##### Error in lightout.cc! This is not regular series!" << std::endl;
    std::cerr << "Series number: " << seriesNo << std::endl;
    std::cerr << "Description: " << desc << std::endl;
    return 0;    
  }
  
  int npoints = data->GetNpoints();
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
    return 0;
  }

  //----- 
  lout->CalculateLightOut(0);
  lout->CalculateLightOut(1);
  lout->CalculateLightOut();
  
  TGraphErrors *gLightOutCh0 = lout->GetLightOutputGraph(0);
  TGraphErrors *gLightOutCh1 = lout->GetLightOutputGraph(1);
  TGraphErrors *gLightOut    = lout->GetLightOutputGraph();
  std::vector <double> lightOutCh0 = lout->GetLightOutput(0);
  std::vector <double> lightOutCh1 = lout->GetLightOutput(1);
  std::vector <double> lightOut    = lout->GetLightOutput();
  std::vector <TH1D*>  specCh0     = lout->GetSpectra(0);
  std::vector <TH1D*>  specCh1     = lout->GetSpectra(1);
  
  //-----
  lout->CalculateLightCol(0);
  lout->CalculateLightCol(1);
  lout->CalculateLightCol();
  
  TGraphErrors *gLightColCh0 = lout->GetLightColGraph(0);
  TGraphErrors *gLightColCh1 = lout->GetLightColGraph(1);
  TGraphErrors *gLightCol    = lout->GetLightColGraph();
  std::vector <double> lightColCh0 = lout->GetLightCol(0);
  std::vector <double> lightColCh1 = lout->GetLightCol(1);
  std::vector <double> lightCol    = lout->GetLightCol();
  
  //----- drawing
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.04);
  
  //----- light output channels 0 and 1
  TCanvas *can_lout_ch = new TCanvas("can_lout_ch", "can_lout_ch", 1200, 600);
  can_lout_ch->Divide(2,1);
  
  can_lout_ch->cd(1);
  gPad->SetGrid(1,1);
  gLightOutCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightOutCh0[0], lightOutCh0[1]));
  
  can_lout_ch->cd(2);
  gPad->SetGrid(1,1);
  gLightOutCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightOutCh1[0], lightOutCh1[1]));
  
  //----- light output summed
  TCanvas *can_lout = new TCanvas("can_lout","can_lout", 700, 500);
  can_lout->cd();
  gPad->SetGrid(1,1);
  gLightOut->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightOut[0], lightOut[1]));

   //----- light collection channels 0 and 1
  TCanvas *can_lcol_ch = new TCanvas("can_lcol_ch", "can_lcol_ch", 1200, 600);
  can_lcol_ch->Divide(2,1);
  
  can_lcol_ch->cd(1);
  gPad->SetGrid(1,1);
  gLightColCh0->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightColCh0[0], lightColCh0[1]));
  
  can_lcol_ch->cd(2);
  gPad->SetGrid(1,1);
  gLightColCh1->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightColCh1[0], lightColCh1[1]));
  
  //----- light output summed
  TCanvas *can_lcol = new TCanvas("can_lcol","can_lcol", 700, 500);
  can_lcol->cd();
  gPad->SetGrid(1,1);
  gLightCol->Draw("AP");
  text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) ph/MeV", lightCol[0], lightCol[1]));
  
  //----- spectra
  TCanvas *can_spec_ch0 = new TCanvas("can_spec_ch0", "can_spec_ch0", 1200, 1200);
  can_spec_ch0->DivideSquare(npoints);
  
  TCanvas *can_spec_ch1 = new TCanvas("can_spec_ch1", "can_spec_ch1", 1200, 1200);
  can_spec_ch1->DivideSquare(npoints);
  
  for(int i=0; i<npoints; i++){
    can_spec_ch0->cd(i+1);
    gPad->SetGrid(1,1);
    specCh0[i]->SetStats(false);
    specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh0[i]->GetYaxis()->SetTitle("counts");
    specCh0[i]->SetTitle(Form("Cherge spectrum Ch0, position %.2f mm", positions[i]));
    specCh0[i]->Draw();
    
    can_spec_ch1->cd(i+1);
    gPad->SetGrid(1,1);
    specCh1[i]->SetStats(false);
    specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    specCh1[i]->GetYaxis()->SetTitle("counts");
    specCh1[i]->SetTitle(Form("Cherge spectrum Ch1, position %.2f mm", positions[i]));
    specCh1[i]->Draw();
  }
  
  //----- saving
  TString path = std::string(getenv("SFPATH"));
  TString fname = Form("lightout_series%i.root", seriesNo);
  
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
      std::cerr << "##### Error in lightout.cc! Unable to create new direcotry!" << std::endl;
      return 0;
    }
  }

  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");
  
  if(!file->IsOpen()){
    std::cerr << "##### Error in lightout.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 0;
  }
  
  can_lout_ch->Write();
  can_lout->Write();
  can_lcol_ch->Write();
  can_lcol->Write();
  can_spec_ch0->Write();
  can_spec_ch1->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "LIGHT_OUTPUT";
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, LOUT, LOUT_ERR, LOUT_CH0, LOUT_CH0_ERR, LOUT_CH1, LOUT_CH1_ERR, LCOL, LCOL_ERR, LCOL_CH0, LCOL_CH0_ERR, LCOL_CH1, LCOL_CH1_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), lightOut[0], lightOut[1], lightOutCh0[0], lightOutCh0[1], lightOutCh1[0], lightOutCh1[1], lightCol[0], lightCol[1], lightColCh0[0], lightColCh0[1], lightColCh1[0], lightColCh1[1]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  delete data;
  delete lout;

  return 1;
} 
