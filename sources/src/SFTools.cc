// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFTools.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFTools.hh"

ClassImp(SFTools);

//------------------------------------------------------------------
/// Returns index in the fNames and fPositions arrays for the 
/// measurement of requested source position in mm. 
/// If measurements in analyzed series don't have unique positions 
/// a number of measurement should be passed. Measurements counting 
/// starts at 1.
/// \param positions - vector containing all source position in the 
/// measurement series
/// \param position - source position of the requested measurement
int SFTools::GetIndex(std::vector <double> positions, double position){
 
    
    int index = -1;
    int npoints = positions.size();
    
    if(fabs(positions[0]-positions[1])<1E-10){
      index = position-1;
      return index;
    }
    
    for(int i=0; i<npoints; i++){
      if(fabs(positions[i]-position)<1){
        index = i;
        break;
      }
    }
    
    if(index==-1){
      std::cerr << "##### Error in SFTools::GetIndex()! Incorrect position!" << std::endl;
      std::abort();
    }
    
    return index;
}
//------------------------------------------------------------------
int SFTools::GetSeriesNo(TString hname_tstr){

    std::string hname = std::string(hname_tstr);
    int nletters = hname.length();
    char letters[nletters];
    strcpy(letters,hname.c_str());
    
    int iposition = -1;
    for(int i=0; i<nletters; i++){
      if(letters[i]=='_'){
          iposition = i;
          break;
      }
    }
    
    if(iposition==-1){
      std::cerr << "##### Error in SFTools::GetSeriesNo()!" << std::endl;
      std::cerr << "Cannot interpret spectrum name!" << std::endl;
      std::abort();
    }
    
    TString seriesName = std::string(&letters[1], &letters[iposition]);
    int seriesNo = atoi(seriesName);
    
    return seriesNo;
} 
//------------------------------------------------------------------
double SFTools::GetPosError(TString collimator, TString testBench){
 
  double err = -1;
  
  if(collimator.Contains("Lead"))
    err = 2.0;   // mm
  else if(collimator.Contains("Electronic") && testBench.Contains("Aachen"))
    err = 1.5;   // mm
  else if(collimator.Contains("Electronic") && testBench.Contains("Krakow"))
    err = 1.0;   // mm

  if(fabs(err+1)<1E-10){
    std::cerr << "##### Error in SFTools::GetPosError()!" << std::endl;
    std::cerr << "Incorrect position error value!" << std::endl;
    std::abort();
  }
    
  return err;
}
//------------------------------------------------------------------
void SFTools::CheckDBStatus(int status, sqlite3 *database){
  
  if(!((status==SQLITE_OK) || (status==SQLITE_DONE))){
      std::cerr << "##### SQL Error in SFTools::CheckDBStatus: " << sqlite3_errmsg(database) << std::endl;  
      std::abort();
    }
    
  return;
}
//------------------------------------------------------------------
bool SFTools::SaveResultsDB(TString database, TString table, 
                            TString query, int seriesNo){
 
  std::cout << "----- Saving results in the databse: " << database << std::endl;
  std::cout << "----- Accessing table: " << table << std::endl;
  
  sqlite3 *resultsDB;
  int status = -1;
  sqlite3_stmt *statement;
  time_t now = time(nullptr);
  ctime(&now);
  
  //--- opening data base
  status = sqlite3_open(database, &resultsDB);
  CheckDBStatus(status, resultsDB);
  
  //--- executing given query
  status = sqlite3_prepare_v2(resultsDB, query, -1, &statement, nullptr);
  CheckDBStatus(status, resultsDB);
  status = sqlite3_step(statement);
  CheckDBStatus(status, resultsDB);
  
  //--- adding time
  TString time_query = Form("UPDATE %s SET DATE = %lld WHERE SERIES_ID = %i", table.Data(), (long long) now, seriesNo);
  status = sqlite3_prepare_v2(resultsDB, time_query, -1, &statement, nullptr);
  CheckDBStatus(status, resultsDB);
  status = sqlite3_step(statement);
  CheckDBStatus(status, resultsDB);
  
  //--- closing data base
  status = sqlite3_close_v2(resultsDB);
  CheckDBStatus(status, resultsDB);
  
  return true;
}
//------------------------------------------------------------------
