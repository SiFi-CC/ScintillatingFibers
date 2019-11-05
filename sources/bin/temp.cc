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
 
  if(argc<2 || argc>6){
    std::cout << "to run type: ./temp seriesNo ";
    std::cout << "-out path/to/output -db database" << std::endl;
    return 1;
  }
    
  int seriesNo = atoi(argv[1]);

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
  std::vector <int> times = data->GetStartTimes();
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
  
  std::vector <int> colors = {608, 806, 409, 856};
  std::vector <SFTempSensor> sensors = {SFTempSensor::Sensor_Ch01, SFTempSensor::Sensor_Ch10,
                                        SFTempSensor::Sensor_Ref, SFTempSensor::Sensor_Out};
  int nsensors = 4;
  std::vector <TGraphErrors*> gTemp;
  std::vector <TGraphErrors*> gTempAv;
  std::vector <double> avTemp;
  std::vector <double> avTempErr;
  std::vector <double> results;
  
  TString textStr;
  TLatex text;
  text.SetTextSize(0.03);
  text.SetTextFont(42);
  text.SetNDC(true);
  
  TCanvas *can_temp = new TCanvas("can_temp", "can_temp", 800, 800);
  gPad->SetGrid(1,1);
  
  TCanvas *can_av = new TCanvas("can_av", "can_av", 800, 800);
  gPad->SetGrid(1,1);
  
  TLegend *leg_temp = new TLegend();
  TLegend *leg_av = new TLegend();
  
  for(int i=0; i<nsensors; i++){
      
    temp->CalcAverageTempSeries(sensors[i]);
    results = temp->GetAverageTempSeries(sensors[i]);
    avTemp.push_back(results[0]);
    avTempErr.push_back(results[1]);
    
    temp->BuildTempPlot(sensors[i]);
    gTemp.push_back(temp->GetTempPlot(sensors[i])); 
    gTemp[i]->SetMarkerColor(colors[i]);
    gTemp[i]->SetLineColor(colors[i]);
    leg_temp->AddEntry(gTemp[i], gTemp[i]->GetName(), "PEL");
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
              Form(" #rightarrow %.2f +/- %.2f deg C", avTemp[i], avTempErr[i]);
    text.SetTextColor(colors[i]);
    text.DrawLatex(0.5, 0.6+(0.035*i), textStr);
    
    temp->BuildTempPlotAverage(sensors[i]);
    gTempAv.push_back(temp->GetTempPlotAverage(sensors[i]));
    gTempAv[i]->SetMarkerColor(colors[i]);
    gTempAv[i]->SetLineColor(colors[i]);
    leg_av->AddEntry(gTempAv[i], gTempAv[i]->GetName(), "PEL");
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
  
  TLine line;
  line.SetLineColor(kGray);
  line.SetLineWidth(1);
  double time; 
  can_temp->cd();
  
  for(int i=0; i<npoints; i++){
    time = (times[i]-times[0])/60.;
    line.DrawLine(time,20,time,30);  
  }

  //----- saving
  TString fname = Form("temp_series%i.root", seriesNo);
  TString outdir;
  TString dbase;

  int ret = parse_common_options(argc, argv, outdir, dbase);
  if(ret != 0) 
    exit(ret);

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
  TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, TEMP_446D45, TEMP_ERR_446D45, TEMP_044F45, TEMP_ERR_044F45, TEMP_2BAD44, TEMP_ERR_2BAD44, TEMP_8F1F46, TEMP_ERR_8F1F46) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), avTemp[3], avTempErr[3], avTemp[0], avTempErr[0], avTemp[1], avTempErr[1], avTemp[2], avTempErr[2]);
  SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
  
  return 0;    
}
