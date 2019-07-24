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
