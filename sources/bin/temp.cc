// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               temp.cc                 *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************  

#include "SFTemperature.hh"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
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
    std::cout << "to run type: ./temp seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }

  SFData *data;
  
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in temp.cc!" << std::endl;
    return 1;
  }
  
  int npoints = data->GetNpoints();
  std::vector <double> positions = data->GetPositions();
  //std::vector <int> times = data->GetStartTimes();
  data->Print();
  
  SFTemperature *temp;
  
  try{
    temp = new SFTemperature(seriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    std::cerr << "##### Exception in temp.cc!" << std::endl;
    return 1;
  }
  
  int nsensors = 4;
  
  //----- reading data base
  std::vector <TString> sensorIDs;
  
  TString dbname = std::string(getenv("SFDATA")) + "/DB/ScintFib_2.db";
  sqlite3 *database;
  int status = sqlite3_open(dbname, &database);
  SFTools::CheckDBStatus(status, database);
  
  sqlite3_stmt *statement;
  TString query = Form("SELECT OUTSIDE, REFERENCE, CH0, CH1 FROM TEMP_SENSOR WHERE SERIES_ID = %i", seriesNo);
  status = sqlite3_prepare_v2(database, query, -1, &statement, nullptr);
  SFTools::CheckDBStatus(status, database);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
      const unsigned char *tmp_out = sqlite3_column_text(statement, 0);
      const unsigned char *tmp_ref = sqlite3_column_text(statement, 1);
      const unsigned char *tmp_ch0 = sqlite3_column_text(statement, 2);
      const unsigned char *tmp_ch1 = sqlite3_column_text(statement, 3);
  
      sensorIDs.push_back(std::string(reinterpret_cast<const char*>(tmp_out)));
      sensorIDs.push_back(std::string(reinterpret_cast<const char*>(tmp_ref)));
      sensorIDs.push_back(std::string(reinterpret_cast<const char*>(tmp_ch0)));
      sensorIDs.push_back(std::string(reinterpret_cast<const char*>(tmp_ch1)));
  }
  //-----
  
  std::vector <TString> sensorLoc = {"outside", "reference", "ch0", "ch1"};
  std::vector <int> colors = {608, 806, 409, 856};
  std::vector <TGraphErrors*> gTemp;
  std::vector <TGraphErrors*> gTempAv;
  std::vector <TemperatureResults> results;
  
  TString textStr;
  TLatex text;
  text.SetTextSize(0.03);
  text.SetTextFont(42);
  text.SetNDC(true);
  
  TCanvas *can_temp = new TCanvas("temp", "temp", 800, 800);
  gPad->SetGrid(1,1);
  
  TCanvas *can_av = new TCanvas("temp_ave", "temp_ave", 800, 800);
  gPad->SetGrid(1,1);
  
  TLegend *leg_temp = new TLegend();
  TLegend *leg_av = new TLegend();
  
  for(int i=0; i<nsensors; i++){
      
    temp->CalcAverageTempSeries(sensorIDs[i]);
    results.push_back(temp->GetAverageTempSeries(sensorIDs[i]));
    
    temp->BuildTempPlot(sensorIDs[i]);
    gTemp.push_back(temp->GetTempPlot(sensorIDs[i])); 
    gTemp[i]->SetMarkerColor(colors[i]);
    gTemp[i]->SetLineColor(colors[i]);
    leg_temp->AddEntry(gTemp[i], sensorLoc[i], "PEL");
    std::cout << gTemp[i]->GetName() << "\t" << gTemp[i]->GetN() << std::endl;
    can_temp->cd();
    
    if(i==0){
      gTemp[0]->SetTitle(Form("Temperature during S%i", seriesNo));
      gTemp[0]->GetYaxis()->SetRangeUser(20,30);
      gTemp[0]->Draw("AP");
    }
    else
      gTemp[i]->Draw("P");
    
    textStr = std::string(gTemp[i]->GetName()) + 
              Form(" #rightarrow %.2f +/- %.2f deg C", 
                   results[i].fTemp, results[i].fTempErr);
    text.SetTextColor(colors[i]);
    text.DrawLatex(0.5, 0.6+(0.035*i), textStr);
    
    temp->BuildTempPlotAverage(sensorIDs[i]);
    gTempAv.push_back(temp->GetTempPlotAverage(sensorIDs[i]));
    gTempAv[i]->SetMarkerColor(colors[i]);
    gTempAv[i]->SetLineColor(colors[i]);
    leg_av->AddEntry(gTempAv[i], sensorLoc[i], "PEL");
    can_av->cd();

    if(i==0){
      gTempAv[0]->SetTitle(Form("Average temperature per mesurement in S%i", seriesNo));  
      gTempAv[0]->GetYaxis()->SetRangeUser(20,30);
      gTempAv[0]->Draw("AP");
    }
    else
      gTempAv[i]->Draw("P");
  }
  
  can_temp->cd();
  leg_temp->Draw();
  
  can_av->cd();
  leg_av->Draw();
  
  //TLine line;
  //line.SetLineColor(kGray);
  //line.SetLineWidth(1);
  //double time; 
  //can_temp->cd();
  
  //for(int i=0; i<npoints; i++){
  //  time = (times[i]-times[0])/60.;
  //  line.DrawLine(time,20,time,30);  
  //}

  //----- saving
  TString fname = Form("temp_series%i.root", seriesNo);
  TString fname_full = outdir + "/" + fname;
  TString dbname_full = outdir + "/" + dbase;
  
  TFile *file = new TFile(fname_full, "RECREATE");

  if(!file->IsOpen()){
    std::cerr << "##### Error in temp.cc!" << std::endl;
    std::cerr << "Couldn't open file: " << fname_full << std::endl;
    return 1;
  }
  
  can_av->Write();
  can_temp->Write();
  file->Close();
  
  //----- writing results to the data base
  TString table = "TEMPERATURE";
  query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, OUT_ID, OUT_TEMP, OUT_ERR, REF_ID, REF_TEMP, REF_ERR, CH0_ID, CH0_TEMP, CH0_ERR, CH1_ID, CH1_TEMP, CH1_ERR) VALUES (%i, '%s', '%s', %f, %f, '%s', %f, %f, '%s', %f, %f, '%s', %f, %f)", table.Data(), seriesNo, fname_full.Data(), sensorIDs[0].Data(), results[0].fTemp, results[0].fTempErr, sensorIDs[1].Data(), results[1].fTemp, results[1].fTempErr, sensorIDs[2].Data(), results[2].fTemp, results[2].fTempErr, sensorIDs[3].Data(), results[3].fTemp, results[3].fTempErr);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  return 0;    
}
