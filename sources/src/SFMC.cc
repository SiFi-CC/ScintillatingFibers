// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               SFMC.hh                 *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************


#include "SFMC.hh"

ClassImp(SFMC);

//------------------------------------------------------------------
char   *MCPath = getenv("SFDATA");
//------------------------------------------------------------------
TString SFMC::MCNames_1[9] = {"LuAG_K/K_T_P10mm","LuAG_K/K_T_P20mm","LuAG_K/K_T_P30mm",
			       "LuAG_K/K_T_P40mm","LuAG_K/K_T_P50mm","LuAG_K/K_T_P60mm",
			       "LuAG_K/K_T_P70mm","LuAG_K/K_T_P80mm","LuAG_K/K_T_P90mm"};
//------------------------------------------------------------------
double SFMC::MCPositions_1[9] = {10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};

//------------------------------------------------------------------
///Default constructor. If this constructor is used the series 
///number should be set via SetDetails(int seriesNo) function.
SFMC::SFMC(){
 Reset();
 cout << "##### Warning in SFMC constructor!" << endl;
 cout << "You are using the default constructor. Required to set the series number manually." <<endl;
}
//------------------------------------------------------------------
///Standard constructor (recommended).
///\param seriesNo is number of simulated series.
SFMC::SFMC(int seriesNo){
  SetDetails(seriesNo);
}
//------------------------------------------------------------------
///Default destructor.
SFMC::~SFMC(){
}
//------------------------------------------------------------------
///Sets all details of selected simulated series. If default constructor
///was used, this function needs to be called explicitly with the number of 
///of requested series as an argument. The following attributes are set 
///within this function:
bool SFMC::SetDetails(int seriesNo){
 
  if(seriesNo!=1 ){
    cout << "##### Error in SFMC::SetDetails()!" << endl;
    cout << "There is no simulation with the number " << seriesNo << " in the database!" << endl;
    return false;
  }
  
  Reset();
  MCSeriesNo = seriesNo;
  
  ///- number of measurements in the series
  if(MCSeriesNo==1) MCNpoints = 6;

  ///- fiber type
  if(MCSeriesNo==1) MCFiber = "LuAG:Ce";
  
  ///- measurement description
  switch(MCSeriesNo){
    case 1: MCDesc = "Simulation of LuAG:Ce";   break;
  }
  
  ///- simulated positions of the radioactive source for all measurements
  switch(MCSeriesNo){
    case 1: MCPositions = MCPositions_1; break;
  }
  
  ///- names of simulated measurements in this series 
  switch(MCSeriesNo){
    case 1: MCNames = MCNames_1; break;
  }
  //~ cout << "SFMC: Finished SetDetails()" << endl;
  return true;
}
//------------------------------------------------------------------
/// Returns proper selection for GetSpectra. Loads the names for the coresponding histogramms from the files. 
/// \param type - type of spectrum to be drawn. Possible options are:
/// Else, Compton, PhotonPeak 
TString SFMC::GetSelection(TString type,TString cut){
 
  TString selection;
  
  if(type=="Compton"){
	if(cut =="511") selection = "C511";
	else if (cut=="1275") selection = "C1275";
	else cout << "The cut on the Compton spectrum is not properly defined.Please set meaningful cut." << endl;
  }
  else if(type=="Else") selection = "Else";
  else if(type=="Photon"){
    if(cut =="511") selection = "P511";
	else if (cut=="1275") selection = "P1275";
	else cout << "The cut on the Compton spectrum is not properly defined.Please set meaningful cut." << endl;
  }
  else{
    cout << "##### Error in SFMC::GetSelection()! Incorrect type!" << endl;
    cout << "Possible options are: Compton, Photon, Else" << endl;
    return "";
  }
  
  return selection;
}
//------------------------------------------------------------------
/// Returns index in the MCNames and MCPositions arrays for the 
/// measurement of requested source position in mm. 
/// If measurements in analyzed series don't have unique positions 
/// a number of measurement should be passed. Measurements counting 
/// starts at 1.
int SFMC::GetIndex(double position){
   
  int index = -1;
  
  if(MCSeriesNo==1){
    index = int(position/10)-1;
  }

  if(index==-1){
   cout << "##### Error in SFMC::GetIndex()! Incorrect position!" << endl;
  }
  
  return index;
}
//------------------------------------------------------------------
/// Returns single spectrum of requested type.
/// \param type - type of the spectrum. Possible options are: Compton, Photon, Else
/// \param cut - energy cut : 511, 1275
/// \param position - simulated position of the source in mm.  
///
TH1D* SFMC::GetSpectrum(TString type, TString cut, double position){

  int index = GetIndex(position);
  TString fname = string(MCPath)+"MC/"+MCNames[index]+".root";
  TFile *MCfile = new TFile(fname,"READ");
  
  TString selection = GetSelection(type,cut);
  MCSpectrum = (TH1D*)((TH1D*)MCfile->FindObjectAny(selection))->Clone();
  TString hname = MCFiber+selection+Form("_pos%.1f",position);
  MCSpectrum->SetName(hname);
  MCSpectrum->SetTitle(hname);
  return MCSpectrum; 
}
//------------------------------------------------------------------
/// Returns a vector with all spectra of requested type.
/// \param type - type of the spectrum. Possible options are: Compton, Photon, Else
/// \param cut - energy cut : 511, 1275
///
vector <TH1D*> SFMC::GetSpectra(TString type, TString cut){
  MCSpectra.clear();
  for(int i=0; i<MCNpoints; i++){
	MCSpectra.push_back(GetSpectrum(type,cut,MCPositions[i]));
  }
  
  return MCSpectra;
}
//------------------------------------------------------------------
///Resets all members of the class to their default values
void SFMC::Reset(void){
 MCSeriesNo      = 0;
 MCNpoints       = 0;
 MCFiber         = "dummy";
 MCDesc          = "dummy"; 
 MCNames         = NULL;
 MCPositions     = NULL;
 MCSpectrum      = NULL;
 if(!MCSpectra.empty()) MCSpectra.clear();
 }
//------------------------------------------------------------------
///Prints details of currently analyzed experimental series
void SFMC::Print(void){
 cout << "\n\n------------------------------------------------" << endl;
 cout << "This is Print() for SFMC class object" << endl;
 cout << "Number of the experimental series: " << MCSeriesNo << endl;
 cout << MCDesc << endl;
 cout << "Number of measurements in this series: " << MCNpoints << endl;
 cout << "Fiber: " << MCFiber << endl;
 cout << "List of measurements in this series:" << endl;
 for(int i=0; i<MCNpoints; i++){
  cout << setw(20);
  cout << MCNames[i] << "\t" << Form("%.1f",MCPositions[i]) << " mm" << endl; 
 }
 cout << "\n" << endl;
}
//------------------------------------------------------------------
