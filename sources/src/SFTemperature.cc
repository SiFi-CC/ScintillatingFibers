// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTemperature.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// ***************************************** 

#include "SFTemperature.hh"

ClassImp(SFTemperature);

//------------------------------------------------------------------
SFTemperature::SFTemperature(int seriesNo): fSeriesNo(seriesNo),
                                            fAvTemp_446D45(-1),
                                            fAvTemp_044F45(-1),
                                            fAvTemp_2BAD44(-1),
                                            fAvTemp_8F1F46(-1),
                                            fAvTempErr_446D45(-1),
                                            fAvTempErr_044F45(-1),
                                            fAvTempErr_2BAD44(-1),
                                            fAvTempErr_8F1F46(-1),
                                            fData(nullptr),
                                            fTempPlot_446D45(nullptr),
                                            fTempPlot_044F45(nullptr),
                                            fTempPlot_2BAD44(nullptr),
                                            fTempPlot_8F1F46(nullptr),
                                            fTempPlotAv_446D45(nullptr),
                                            fTempPlotAv_044F45(nullptr),
                                            fTempPlotAv_2BAD44(nullptr),
                                            fTempPlotAv_8F1F46(nullptr) {
                                                
  try{
    fData = new SFData(fSeriesNo);
  }
  catch(const char *message){
    std::cerr << message << std::endl;
    throw "##### Exception in SFTemperature constructor!";
  }
  
  TString tempFile = fData->GetTempFile();
  
  if(tempFile=="-"){
    std::cerr << "This series doesn't have temeperature measurement!" <<std::endl;
    throw "##### Exception in SFTemperature!";
  }

  OpenDataBase();
   
}
//------------------------------------------------------------------
SFTemperature::~SFTemperature(){
    
  int status = sqlite3_close(fDB);
  if(status!=0) 
    std::cerr << "In SFData destructor. Data base corrupted!" << std::endl;
  }
//------------------------------------------------------------------
bool SFTemperature::OpenDataBase(void){
    
  const char *path = getenv("SFDATA");
  TString  DBname = std::string(path) + "/DB/ScintFib_2.db";
  int status = sqlite3_open(DBname, &fDB);
  
  if(status!=0){
   std::cerr << "##### Error in SFTemperature::OpenDataBase()!" << std::endl;
   std::cerr << "Could not access data base!" << std::endl;
   return false;
 }
 
 return true;
}
//------------------------------------------------------------------
TString SFTemperature::GetSensorID(SFTempSensor sensor){
    
   TString sensorID = "";
   
   switch(sensor){
       case SFTempSensor::Sensor_Ch01:
           sensorID = "044F45";
           break;
       case SFTempSensor::Sensor_Ch10:
           sensorID = "2BAD44";
           break;
       case SFTempSensor::Sensor_Ref:
           sensorID = "8F1F46";
           break;
       case SFTempSensor::Sensor_Out:
           sensorID = "446D45";
           break;
       default:
           std::cerr << "##### Error in SFTemperature::GetSensorID()!" << std::endl;
           std::cerr << "##### Unknown temperature sensor type!" << std::endl;
           break;
   }
   
   return sensorID;
}
//------------------------------------------------------------------
bool SFTemperature::LoadFromDB(SFTempSensor sensor, TString name){
 
  int status = 0;  
  sqlite3_stmt *statement;
  TString sensorID = GetSensorID(sensor);
  
  if(!fTime.empty())
    fTime.clear();
  if(!fTemperatures.empty())
    fTemperatures.clear();
  
  TString query;
  
  if(name=="none")
    query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID = %i AND SENSOR_ID = '28.%s0B0000'", fSeriesNo, sensorID.Data());
  else
    query = query = Form("SELECT TIME, TEMPERATURE FROM TEMPERATURES WHERE SERIES_ID = %i AND SENSOR_ID = '28.%s0B0000' AND MEASUREMENT_NAME = '%s'", fSeriesNo, sensorID.Data(), name.Data());
  
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    fTime.push_back(sqlite3_column_int(statement, 0));
    fTemperatures.push_back(sqlite3_column_double(statement, 1));
  }
  
  SFTools::CheckDBStatus(status, fDB);
  sqlite3_finalize(statement);
  
  return true;    
}
//------------------------------------------------------------------
std::vector <double> SFTemperature::CalcAverageTempMeasure(SFTempSensor sensor, TString name){
 
  LoadFromDB(sensor, name);
  double mean = SFTools::GetMean(fTemperatures);
  double stdErr = SFTools::GetStandardErr(fTemperatures);
  std::cout << mean << "\t" << stdErr << std::endl;
  std::vector <double> result;
  result.push_back(mean);
  result.push_back(stdErr);
  
  return result;
}
//------------------------------------------------------------------
bool SFTemperature::CalcAverageTempSeries(SFTempSensor sensor){
    
  LoadFromDB(sensor);
  double mean = SFTools::GetMean(fTemperatures);
  double stDev = SFTools::GetStandardErr(fTemperatures); 
  
  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          fAvTemp_044F45 = mean;
          fAvTempErr_044F45 = stDev;
          break;
      case SFTempSensor::Sensor_Ch10:
          fAvTemp_2BAD44 = mean;
          fAvTempErr_2BAD44 = stDev;
          break;
      case SFTempSensor::Sensor_Ref:
          fAvTemp_8F1F46 = mean;
          fAvTempErr_8F1F46 = stDev;
          break;
      case SFTempSensor::Sensor_Out:
          fAvTemp_446D45 = mean;
          fAvTempErr_446D45 = stDev;
          break;
      default:
          std::cerr << "##### Error in SFTemperature::CalcAverageTempSeries()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
  
  return true;
}
//------------------------------------------------------------------
bool SFTemperature::BuildTempPlotAverage(SFTempSensor sensor){
 
  TString sensorID = GetSensorID(sensor);  
  int npoints = fData->GetNpoints();
  TString collimator = fData->GetCollimator();
  TString testBench = fData->GetTestBench();
  std::vector <TString> names = fData->GetNames();
  std::vector <double> pos = fData->GetPositions();
  std::vector <double> result;
  
  TGraphErrors *gr = new TGraphErrors(npoints);
  gr->SetName(sensorID);
  gr->SetMarkerStyle(4);
  gr->GetXaxis()->SetTitle("source position");
  gr->GetYaxis()->SetTitle("average temperature during measurement [deg C]");
    
  for(int i =0; i<npoints; i++){
    result = CalcAverageTempMeasure(sensor, names[i]);
    gr->SetPoint(i, pos[i], result[0]);
    gr->SetPointError(i, SFTools::GetPosError(collimator, testBench), result[1]);
  }
  
  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          fTempPlotAv_044F45 = gr;
          break;
      case SFTempSensor::Sensor_Ch10:
          fTempPlotAv_2BAD44 = gr;
          break;
      case SFTempSensor::Sensor_Ref:
          fTempPlotAv_8F1F46 = gr;
          break;
      case SFTempSensor::Sensor_Out:
          fTempPlotAv_446D45 = gr;
          break;
      default:
          std::cerr << "##### Error in SFTemperature::BuildTempPlotAverage()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
          
  return true;
}
//------------------------------------------------------------------
bool SFTemperature::BuildTempPlot(SFTempSensor sensor){
  
  LoadFromDB(sensor);
  int size = fTime.size();  
  TString sensorID = GetSensorID(sensor);  
  double time;
  std::cout << "size: " << size << std::endl;
  
  TGraphErrors *gr = new TGraphErrors(size);
  gr->SetName(sensorID);
  gr->SetMarkerStyle(4);
  gr->GetXaxis()->SetTitle("minutes since beginning of measurement");
  gr->GetYaxis()->SetTitle("temperature [deg C]");
  
  for(int i=0; i<size; i++){
      time = (fTime[i]-fTime[0])/60.;
      gr->SetPoint(i, time, fTemperatures[i]);
  }
  
  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          fTempPlot_044F45 = gr;
          break;
      case SFTempSensor::Sensor_Ch10:
          fTempPlot_2BAD44 = gr;
          break;
      case SFTempSensor::Sensor_Ref:
          fTempPlot_8F1F46 = gr;
          break;
      case SFTempSensor::Sensor_Out:
          fTempPlot_446D45 = gr;
          break;
      default:
          std::cerr << "##### Error in SFTemperature::BuildTempPlot()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
  
  return true;
}
//------------------------------------------------------------------
std::vector <double> SFTemperature::GetAverageTempSeries(SFTempSensor sensor){
    
  std::vector <double> temp;

  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          temp.push_back(fAvTemp_044F45);
          temp.push_back(fAvTempErr_044F45);
          break;
      case SFTempSensor::Sensor_Ch10:
          temp.push_back(fAvTemp_2BAD44);
          temp.push_back(fAvTempErr_2BAD44);
          break;
      case SFTempSensor::Sensor_Ref:
          temp.push_back(fAvTemp_8F1F46);
          temp.push_back(fAvTempErr_8F1F46);
          break;
      case SFTempSensor::Sensor_Out:
          temp.push_back(fAvTemp_446D45);
          temp.push_back(fAvTempErr_446D45);
          break;
      default:
          std::cerr << "##### Error in SFTemperature::GetAverageTempSeries()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
  
  return temp;
}
//------------------------------------------------------------------
TGraphErrors *SFTemperature::GetTempPlot(SFTempSensor sensor){
    
  TGraphErrors *temp;

  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          temp = fTempPlot_044F45;
          break;
      case SFTempSensor::Sensor_Ch10:
          temp = fTempPlot_2BAD44;
          break;
      case SFTempSensor::Sensor_Ref:
          temp = fTempPlot_8F1F46;
          break;
      case SFTempSensor::Sensor_Out:
          temp = fTempPlot_446D45;
          break;
      default:
          std::cerr << "##### Error in SFTemperature::GetTempPlot()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
  
  if(temp==nullptr){
    std::cerr << "##### Error in SFTemperature::GetTempPlot()!" << std::endl;
    std::cerr << "Empty graph! Please check!" << std::endl;
    std::abort();
  }
  
  return temp;
}
//------------------------------------------------------------------
TGraphErrors *SFTemperature::GetTempPlotAverage(SFTempSensor sensor){
  
  TGraphErrors *temp;
  
  switch(sensor){
      case SFTempSensor::Sensor_Ch01:
          temp = fTempPlotAv_044F45;
          break;
      case SFTempSensor::Sensor_Ch10:
          temp = fTempPlotAv_2BAD44;
          break;
      case SFTempSensor::Sensor_Ref:
          temp = fTempPlotAv_8F1F46;
          break;
      case SFTempSensor::Sensor_Out:
          temp = fTempPlotAv_446D45;
          break;
      default:
          std::cerr << "##### Error in SFTemperature::GetTempPlotAverage()!" << std::endl;
          std::cerr << "Unknown temperature sensor! Please check!" << std::endl;
          break;
  }
  
  if(temp==nullptr){
    std::cerr << "##### Error in SFTemperature::GetTempPlotAverage()!" << std::endl;
    std::cerr << "Empty graph! Please check!" << std::endl;
    std::abort();
  }
  
  return temp;    
}
//------------------------------------------------------------------
void SFTemperature::Print(void){
  std::cout << "\n-------------------------------------------" << std::endl;
  std::cout << "This is print out of SFTemperature class object" << std::endl;
  std::cout << "Experimental series number: " << fSeriesNo << std::endl;
  std::cout << "-------------------------------------------\n" << std::endl;
  return;
}
//------------------------------------------------------------------
