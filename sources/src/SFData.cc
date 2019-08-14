// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFData.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"

ClassImp(SFData);

//------------------------------------------------------------------
// constants
static const char  *gPath = getenv("SFDATA");  // path to the experimental data
static const int    gBaselineMax = 50;         // number of samples for base line determination
static const double gmV          = 4.096;      // coefficient to calibrate ADC channels to mV
//------------------------------------------------------------------
/// Default constructor. If this constructor is used the series 
/// number should be set via SetDetails(int seriesNo) function.
SFData::SFData(): fSeriesNo(-1),
                  fNpoints(-1),
                  fFiber("dummy"),
                  fSource("dummy"),
                  fCollimator("dummy"),
                  fDesc("dummy"),
                  fTestBench("dummy"),
                  fSiPM("dummy"),
                  fSpectrum(nullptr),
                  fHist(nullptr),
                  fHist2D(nullptr),
                  fSignalProfile(nullptr),
                  fSignal(nullptr) {
                      
 std::cout << "##### Warning in SFData constructor!" << std::endl;
 std::cout << "You are using the default constructor. Set the series number & open data base!" << std::endl;
}
//------------------------------------------------------------------ 
/// Standard constructor (recommended).
/// \param seriesNo is number of experimental series to analyze.
SFData::SFData(int seriesNo): fSeriesNo(seriesNo),
                              fNpoints(-1),
                              fFiber("dummy"),
                              fSource("dummy"),
                              fCollimator("dummy"),
                              fDesc("dummy"),
                              fTestBench("dummy"),
                              fSiPM("dummy"),
                              fSpectrum(nullptr),
                              fHist(nullptr),
                              fHist2D(nullptr),
                              fSignalProfile(nullptr),
                              fSignal(nullptr) {
                                  
 bool db_stat  = OpenDataBase("ScintFib_2.db");
 bool set_stat = SetDetails(seriesNo);
 if(!db_stat || !set_stat){
   throw "##### Exception in SFData constructor!";
 }
}
//------------------------------------------------------------------
/// Default destructor.
SFData::~SFData(){
    
 int status = sqlite3_close(fDB);
 if(status!=0) 
   std::cerr << "In SFData destructor. Data base corrupted!" << std::endl;
}
//------------------------------------------------------------------
/// Opens SQLite3 data base containing details of experimental series
/// and measurements.
bool SFData::OpenDataBase(TString name){

 TString db_name = std::string(gPath) + "/DB/" + name;
 int status = sqlite3_open(db_name, &fDB);
 
 if(status!=0){
   std::cerr << "##### Error in SFData::OpenDataBase()!" << std::endl;
   std::cerr << "Could not access data base!" << std::endl;
   return false;
 }
 
 return true;
}
//------------------------------------------------------------------
/// Sets all details of selected experimental series. If default constructor
/// was used, this function needs to be called explicitly with the number of 
/// of requested series as an argument. 
/// \param seriesNo - number of experimental series
/// 
/// The following attributes are set within this function:
bool SFData::SetDetails(int seriesNo){
  
  TString query;
  sqlite3_stmt *statement;
  int status;
  
  if(fSeriesNo==-1)
    fSeriesNo = seriesNo;
  
  //-----Checking if series number is valid
  int maxSeries;
  query = "SELECT COUNT(*) FROM SERIES";
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    maxSeries = sqlite3_column_int(statement, 0);
  }
  
  SFTools::CheckDBStatus(status, fDB);
  
  sqlite3_finalize(statement);

  if(fSeriesNo<1 || fSeriesNo>maxSeries){
   std::cerr << "##### Error in SFData::SetDetails()! Series number out of range!" << std::endl;
   return false;
  }
  //-----
  
  //----- Setting series attributes
  ///- fiber type 
  ///- radioactive source type
  ///- test bench type
  ///- collimator type 
  ///- number of measurements in the series
  ///- description of the series
  query = Form("SELECT FIBER, SOURCE, TEST_BENCH, COLLIMATOR, SIPM, NO_MEASUREMENTS, DESCRIPTION FROM SERIES WHERE SERIES_ID = %i", fSeriesNo);
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *fiber = sqlite3_column_text(statement, 0);
    const unsigned char *source = sqlite3_column_text(statement, 1);
    const unsigned char *test_bench = sqlite3_column_text(statement, 2);
    const unsigned char *collimator = sqlite3_column_text(statement, 3);
    const unsigned char *sipm = sqlite3_column_text(statement, 4);
    const unsigned char *description = sqlite3_column_text(statement, 6);
    fNpoints = sqlite3_column_int(statement, 5);
    fFiber = std::string(reinterpret_cast<const char*>(fiber));
    fSource = std::string(reinterpret_cast<const char*>(source));
    fDesc = std::string(reinterpret_cast<const char*>(description));
    fCollimator = std::string(reinterpret_cast<const char*>(collimator));
    fTestBench = std::string(reinterpret_cast<const char*>(test_bench));
    fSiPM = std::string(reinterpret_cast<const char*>(sipm));
  }
  
  SFTools::CheckDBStatus(status, fDB);
  
  sqlite3_finalize(statement);
  //-----
  
  //----- Setting measurements attributes
  ///- list of measurements names
  ///- list of measurements duration times
  ///- list of source positions
  ///- list of measurements starting times
  ///- list of measurements stopping times
  query = Form("SELECT MEASUREMENT_NAME, DURATION_TIME, SOURCE_POSITION, START_TIME, STOP_TIME FROM MEASUREMENT WHERE SERIES_ID = %i", fSeriesNo);
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *name = sqlite3_column_text(statement, 0);
    fNames.push_back(std::string(reinterpret_cast<const char*>(name)));
    fTimes.push_back(sqlite3_column_int(statement, 1));
    fPositions.push_back(sqlite3_column_double(statement, 2));
    fStart.push_back(sqlite3_column_int(statement, 3));
    fStop.push_back(sqlite3_column_int(statement, 4));
  }
  
  SFTools::CheckDBStatus(status, fDB);
  
  sqlite3_finalize(statement);
   
  return true;
}
//------------------------------------------------------------------
/// Parses given cut and checks if signal fulfills conditions specified by it. 
/// \param sig - currently analyzed signal, as read from the tree
/// \param cut - a logic cut to select specific signals.
/// The following syntax of the cut is acceptable: 
/// - inequality signs: '<' and '>'
/// - single expressions, e.g. "fAmp>50", "ch_0.fPE>10", "fT0<100"
/// - double expressions with &&, e.g. "fPE>10 && fPE<100", "ch_0.fT0>0 && ch_0.fT0<500" 
///
/// If empty cut is passed returns true.
bool SFData::InterpretCut(DDSignal *sig, TString cut){
 
  bool result = true;
  if(cut=="" || cut==" ") return result;
  
  double amp = sig->GetAmplitude();
  double charge = sig->GetCharge();
  double pe = sig->GetPE();
  double t0 = sig->GetT0();
  double tot = sig->GetTOT();
  
  //convert TString into string and char[]
  std::string cut_str = std::string(cut);
  int nletters = cut_str.length();
  char letters[nletters];
  strcpy(letters,cut_str.c_str());
  
  //splitting cut into expressions
  int iposition = -1;
  int nexpressions = 0;
  std::string expression[2];
  
  for(int i=0; i<nletters; i++){
   if(letters[i]=='&' && letters[i+1]=='&'){
     iposition = i;
     break;
    }
  }
  
  if(iposition==-1){
    nexpressions = 1;
    expression[0] = cut_str;
  }
  else{
    nexpressions = 2;
    expression[0] = std::string(&letters[0], &letters[iposition]);
    expression[1] = std::string(&letters[iposition+2], &letters[nletters]);
  }
    
  //extracting doubles 
  int nletters_expr = 0;
  int istop = -1;
  char letters_expr[100];
  std::string number_str[2];
  double number[2];
  
  for(int i=0; i<nexpressions; i++){
   istop = -1;
   nletters_expr = expression[i].length();
   strcpy(letters_expr, expression[i].c_str());
   for(int ii=nletters_expr; ii>0; ii--){
    if(letters_expr[ii]=='<' || letters_expr[ii]=='>'){
      istop = ii+1;
      break;
    }
   }
   if(istop==-1){
     std::cerr << "#### Error in SFData::InterpretCut! Incorrect cut syntax!" << std::endl;
     std::cerr << "Missing '<' or '>'." << std::endl;
     return false;
   }
   number_str[i] = std::string(&letters_expr[istop], &letters_expr[nletters_expr]);
   number[i] = atof(number_str[i].c_str());
  }

  //checking logic
  bool logic[2] = {true, true};
  
  for(int i=0; i<nexpressions; i++){
    
   if(expression[i].find("fAmp")!=std::string::npos){		//cut on fAmp
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = amp<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
      logic[i] = amp>number[i]; 
     }
   }
   
   else if(expression[i].find("fPE")!=std::string::npos){		//cut on fPE
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = pe<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = pe>number[i]; 
     }
   }
   
   else if(expression[i].find("fCharge")!=std::string::npos){	//cut on fCharge
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = charge<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = charge>number[i]; 
     }
   }
   
   else if(expression[i].find("fT0")!=std::string::npos){		//cut on fT0
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = t0<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = t0>number[i]; 
     }
   } 
   
   else if(expression[i].find("fTOT")!=std::string::npos){		//cut on fTOT
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = tot<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = tot>number[i]; 
     }
   } 
   
   else{
    std::cerr << "#### Error in SFData::InterpretCut! Incorrect cut syntax!" << std::endl;
    std::cerr << "Incorrect type. Available types are: fAmp, fCharge, fPE, fT0 and fTOT." << std::endl;
    return false;
   }
  }
  
  result = logic[0] && logic[1];
  
  return result;
}
//------------------------------------------------------------------
/// Returns single spectrum of requested type.
/// \param ch - chennel number
/// \param sel_type - type of the spectrum. Possible options are: fAmp, fCharge, fPE, fT0 and fTOT
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
/// \param position - position of the source in mm. If analyzed series doesn't have
/// unique positions a number of measurement should be passed here. Numbering starts at 1. 
///
/// It is possible to have spectrum with cut or raw spectrum as recorded. In the latter case pass
/// empty string as cut.
TH1D* SFData::GetSpectrum(int ch, SFSelectionType sel_type, TString cut, double position){

  int index = SFTools::GetIndex(fPositions, position);
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TString tname = std::string("tree_ft");
  TTree *tree = (TTree*)file->Get(tname);
  fSpectrum = new TH1D();
  
  gUnique+=1;
  TString selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch);
  int stat = tree->Draw(selection, cut);
  fSpectrum = (TH1D*)gROOT->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_ch%i_pos%.1f_", fSeriesNo, ch, position)+SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  fSpectrum->SetName(hname);
  fSpectrum->SetTitle(htitle);
  
  return fSpectrum; 
}
//------------------------------------------------------------------
/// Returns a vector with all spectra of requested type.
/// \param ch - channel number
/// \param sel_type - type of spectra (fAmp, fCharge, fPE, fT0, fTOT)
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
///
/// Like with sigle spectrum, it is possible to have cut and raw spectra.
std::vector <TH1D*> SFData::GetSpectra(int ch, SFSelectionType sel_type, TString cut){

  bool empty = fSpectra.empty();
  for(int i=0; i<fNpoints; i++){
   if(empty) fSpectra.push_back(new TH1D());
   if(fDesc.Contains("Regular series"))
     fSpectra[i] = GetSpectrum(ch, sel_type, cut, fPositions[i]);
   else  
     fSpectra[i] = GetSpectrum(ch, sel_type, cut, i+1);
  }
  
  return fSpectra;
}
//------------------------------------------------------------------
/// Returns single requested custom 1D histogram.
/// \param sel_type - predefined selection type
/// \param cut - cut for drawn events. Also TTree-style syntax
/// \param position - position of source in mm. If position is not unique a 
/// number of measurement should be entered.
TH1D* SFData::GetCustomHistogram(SFSelectionType sel_type, TString cut, double position, 
                                std::vector <double> customNumbers){
  
  int index = SFTools::GetIndex(fPositions, position);
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TString tname = "tree_ft";
  TTree *tree = (TTree*)file->Get(tname);
  fHist = new TH1D();
  
  gUnique+=1;
  TString selection; 
  selection = SFDrawCommands::GetSelection(sel_type, gUnique, customNumbers);
  tree->Draw(selection, cut);
  fHist = (TH1D*)gROOT->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_pos%.1f_", fSeriesNo, position) + SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  fHist->SetName(hname);
  fHist->SetTitle(htitle);
  
  return fHist;
}
//------------------------------------------------------------------
/// Returns a vector of requested custom 1D histograms for all measurements in this series.
/// \param sel_type - predefined selection type
/// \param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
/// all events will be drawn.
std::vector <TH1D*> SFData::GetCustomHistograms(SFSelectionType sel_type, TString cut){
  
  bool empty = fHists.empty();
  for(int i=0; i<fNpoints; i++){
    if(empty) fHists.push_back(new TH1D());
    if(fDesc.Contains("Regular series"))
      fHists[i] = GetCustomHistogram(sel_type, cut, fPositions[i]);
    else 
      fHists[i] = GetCustomHistogram(sel_type, cut, i+1);
  }
  
  return fHists;
}
//------------------------------------------------------------------
TH1D* SFData::GetCustomHistogram(int ch, SFSelectionType sel_type, TString cut, double position, 
                                 std::vector <double> customNumbers){
    
  int index = SFTools::GetIndex(fPositions, position);
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TString tname = "tree_ft";
  TTree *tree = (TTree*)file->Get(tname);
  fHist = new TH1D();
  
  gUnique+=1;
  TString selection;
  if(customNumbers.empty())
    selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch);
  else 
    selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch, customNumbers);
  tree->Draw(selection, cut);
  fHist = (TH1D*)gROOT->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_pos%.1f_", fSeriesNo, position)+ SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  fHist->SetName(hname);
  fHist->SetTitle(htitle);
  
  return fHist;
    
}
//------------------------------------------------------------------
/// Returns single requested correlation 2D histogram.
/// \param sel_type - predefined selection type
/// \param cut - cut for drawn events. Also TTree-style syntax
/// \param position - position of source in mm. If position is not unique a 
/// number of measurement should be entered.
TH2D* SFData::GetCorrHistogram(SFSelectionType sel_type, TString cut, double position){
  
  int index = SFTools::GetIndex(fPositions, position);
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TString tname = std::string("tree_ft");
  TTree *tree = (TTree*)file->Get(tname);
  fHist2D = new TH2D();
  
  gUnique+=1;
  TString selection = SFDrawCommands::GetSelection(sel_type, gUnique);
  tree->Draw(selection, cut, "colz");
  fHist2D = (TH2D*)gROOT->FindObjectAny(Form("htemp%.i", gUnique));
  TString hname = Form("S%i_pos%.1f_", fSeriesNo, position) + SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  fHist2D->SetName(hname);
  fHist2D->SetTitle(htitle);
  
  return fHist2D;
}
//------------------------------------------------------------------
/// Returns a vector of requested 2D correlation histograms for all measurements in 
/// this series.
/// \param sel_type - predefined selection type
/// \param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
/// all events will be drawn.
std::vector <TH2D*> SFData::GetCorrHistograms(SFSelectionType sel_type, TString cut){
  
  bool empty = fHists2D.empty();
  for(int i=0; i<fNpoints; i++){
    if(empty) fHists2D.push_back(new TH2D());
    if(fDesc.Contains("Regular series"))
      fHists2D[i] = GetCorrHistogram(sel_type, cut, fPositions[i]);
    else 
      fHists2D[i] = GetCorrHistogram(sel_type, cut, i+1);
  }
  
  return fHists2D;
}
//------------------------------------------------------------------
/// Returns averaged signal.
/// \param ch - channel number 
/// \param position - position of the source in mm. If position is not unique pass measurement number here.
/// Measurement numbering starts at 1.
/// \param cut - logic cut to choose signals. Syntax of this cut is explained in InterpretCut() function
/// \param number - number of signals to be averaged
/// \param bl - flag for base line subtraction. If true - base line will be subtracted, if false - it won't
///
/// If no cut is required pass an empty string.
TProfile* SFData::GetSignalAverage(int ch, double position, TString cut, int number, bool bl){

  TProfile *sig = nullptr;
  
  if(fTestBench=="Krakow"){
    sig = GetSignalAverageKrakow(ch, position, cut, number, bl);
  }
  else if(fTestBench=="Aachen"){
    sig = GetSignalAverageAachen(ch, position, cut, number);    
  }
  else{
    std::cerr << "##### Error in SFData::GetSignalAverage()!" << std::endl;
    std::cerr << "Unknown data format!" << std::endl;
    std::abort();
  }
  
  return sig;
}
//------------------------------------------------------------------
TProfile* SFData::GetSignalAverageKrakow(int ch, double position, TString cut, int number, bool bl){
 
  int index = SFTools::GetIndex(fPositions, position);
  const int ipoints = 1024;
  float x;
    
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = std::string(gPath) + fNames[index] + Form("/wave_%i.dat", ch);
  std::ifstream input(iname, std::ios::binary);
  
  TString hname = "sig_profile";
  TString htitle = "sig_profile";
  fSignalProfile = new TProfile(hname, htitle, ipoints, 0, ipoints, "");
  
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int counter = 0;
  int infile = 0;
  bool condition = true;
  double firstT0 = 0.;
  
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   condition = InterpretCut(sig, cut);
   if(condition && fabs(firstT0)<1E-10) firstT0 = sig->GetT0();
   if(condition && fabs(sig->GetT0()-firstT0)<1){
     infile = sizeof(x)*ipoints*i;
     if(bl){
       input.seekg(infile);
       baseline = 0.;
       for(int ii=0; ii<gBaselineMax; ii++){
         input.read((char*)&x, sizeof(x));
         baseline += x/gmV;
       }
       baseline = baseline/gBaselineMax;
     }
     input.seekg(infile);
     for(int ii=0; ii<ipoints; ii++){
       input.read((char*)&x, sizeof(x));
       if(bl) fSignalProfile->Fill(ii, (x/gmV)-baseline);
       else   fSignalProfile->Fill(ii, (x/gmV));
     }
     if(counter<number) counter++;
     else break;
    }
  }
  
  hname = Form("S%i_ch%i_pos_%.1f_sig_num_%i", fSeriesNo, ch, position, counter);
  htitle = hname + " " + cut;
  fSignalProfile->SetName(hname);
  fSignalProfile->SetTitle(htitle);
  
  if(counter<number){ 
    std::cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
              << " out of " << number << " plotted." << std::endl;
    std::cout << "Position: " << position << "\t channel: " << ch << std::endl; 
  }
  
  input.close();
  
  return fSignalProfile;
}
//------------------------------------------------------------------
TProfile* SFData::GetSignalAverageAachen(int ch, double position, TString cut, int number){
  
  int index = SFTools::GetIndex(fPositions, position);
  const int ipoints = 1024;
  float x;
  
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = std::string(gPath) + fNames[index] + "/waves.root";
  TFile* iFile = new TFile(iname, "READ");
  TTree* iTree = (TTree*)iFile->Get("wavetree");
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  TString bname = Form("voltages_ch_%i", ch);
  iTree->SetBranchAddress(bname, &iVolt);  
  
  TString hname = "sig_profile";
  TString htitle = "sig_profile";
  fSignalProfile = new TProfile(hname, htitle, ipoints, 0, ipoints, "");
  
  int nentries = tree->GetEntries();
  int counter = 0;
  bool condition = true;
  double firstT0 = 0.;
  
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   condition = InterpretCut(sig, cut);
   if(condition && fabs(firstT0)<1E-10) firstT0 = sig->GetT0();
   if(condition && fabs(sig->GetT0()-firstT0)<1){
     iTree->GetEntry(i);
     for(int ii=0; ii<ipoints; ii++){
       fSignalProfile->Fill(ii+1, (*iVolt)[ii]);
     }
     if(counter<number) counter++;
     else break;
    }
  }
  
  hname = Form("S%i_ch%i_pos_%.1f_sig_num_%i", fSeriesNo, ch, position, counter);
  htitle = hname + " " + cut;
  fSignalProfile->SetName(hname);
  fSignalProfile->SetTitle(htitle);
  
  if(counter<number){ 
    std::cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
              << " out of " << number << " plotted." << std::endl;
    std::cout << "Position: " << position << "\t channel: " << ch << std::endl; 
  }
  
  return fSignalProfile;
}
//------------------------------------------------------------------
TH1D* SFData::GetSignal(int ch, double position, TString cut, int number, bool bl){
 
  TH1D *sig = nullptr;
  
  if(fTestBench=="Krakow"){
    sig = GetSignalKrakow(ch, position, cut, number, bl);   
  }
  else if(fTestBench=="Aachen"){
    sig = GetSignalAachen(ch, position, cut, number);
  }
  else{
    std::cerr << "##### Error in SFData::GetSignal()" << std::endl;
    std::cerr << "Unknown data format!" << std::endl;
    std::abort();
  }
  
  return sig;
}
//------------------------------------------------------------------
TH1D* SFData::GetSignalKrakow(int ch, double position, TString cut, int number, bool bl){

  int index = SFTools::GetIndex(fPositions, position);
  const int ipoints = 1024;
  float x; 
  
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = std::string(gPath) + fNames[index] + Form("/wave_%i.dat", ch);
  std::ifstream input(iname, std::ios::binary);
  
  TString hname = Form("S%i_ch%i_pos_%.1f_sig_no%i", fSeriesNo, ch, position, number);
  TString htitle = hname + " " + cut; 
  
  fSignal = new TH1D(hname, htitle, ipoints, 0, ipoints);
  
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int infile = 0;
  int counter = 0;
  bool condition = true;
  
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    condition = InterpretCut(sig, cut);
    if(condition){
      counter++;
        if(counter!=number) continue;
        infile = sizeof(x)*ipoints*i;
        if(bl){
          input.seekg(infile);
          baseline = 0;
          for(int ii=0; ii<gBaselineMax; ii++){
            input.read((char*)&x, sizeof(x));
            baseline += x/gmV;
          }
          baseline = baseline/gBaselineMax;
        }  
        input.seekg(infile);
        for(int ii=1; ii<ipoints+1; ii++){
          input.read((char*)&x, sizeof(x));
          if(bl) fSignal->SetBinContent(ii, (x/gmV)-baseline);
          else   fSignal->SetBinContent(ii, (x/gmV));
        }
    }
  }
  
  input.close();
  
  return fSignal;
}
//------------------------------------------------------------------
/// Returns single signal.
/// \param ch - channel number
/// \param position - source position in mm. If there's no unique position pass a number of measurement here.
/// \param cut - logic cut to choose signals. Syntax of this cut is explained in InterpretCut() function
/// \param number - requested number of the signal to be drawn
/// \param bl - flag for base line subtraction. See GetSignalAverage()
///
/// If no cut is needed an empty string should be passed.
TH1D* SFData::GetSignalAachen(int ch, double position, TString cut, int number){
 
  int index = SFTools::GetIndex(fPositions, position);
  const int ipoints = 1024;
  float x;
  
  TString fname = std::string(gPath) + fNames[index] + "/results.root";
  TFile *file = new TFile(fname, "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = std::string(gPath) + fNames[index] + "/waves.root";
  TFile* iFile = new TFile(iname, "READ");
  TTree* iTree = (TTree*)iFile->Get("wavetree");
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  TString bname = Form("voltages_ch_%i", ch);
  iTree->SetBranchAddress(bname, &iVolt);
  
  TString hname = Form("S%i_ch%i_pos_%.1f_sig_no%i", fSeriesNo, ch, position, number);
  TString htitle = hname + " " + cut;
  
  fSignal = new TH1D(hname, htitle, ipoints, 0, ipoints);
  
  int nentries = tree->GetEntries();
  int counter = 0;
  bool condition = true;
  
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    iTree->GetEntry(i);
    condition = InterpretCut(sig, cut);
    if(condition){
      counter++;
      if(counter!=number) continue;
      for(int ii=0; ii<ipoints; ii++){
        fSignal->SetBinContent(ii+1, (*iVolt)[ii]);
      }
    }
  }
  
  return fSignal;
}
//------------------------------------------------------------------
/// Prints details of currently analyzed experimental series
void SFData::Print(void){
 std::cout << "\n\n------------------------------------------------" << std::endl;
 std::cout << "This is Print() for SFData class object" << std::endl;
 std::cout << "Number of the experimental series: " << fSeriesNo << std::endl;
 std::cout << fDesc << std::endl;
 std::cout << "Collimator: " << fCollimator << std::endl;
 std::cout << "Test bench: " << fTestBench << std::endl;
 std::cout << "Number of measurements in this series: " << fNpoints << std::endl;
 std::cout << "Fiber: " << fFiber << std::endl;
 std::cout << "Radioactive source: " << fSource << std::endl;
 std::cout << "List of measurements in this series:" << std::endl;
 for(int i=0; i<fNpoints; i++){
  std::cout << std::setw(30);
  std::cout << fNames[i] << "\t\t" << Form("%.1f mm", fPositions[i]) 
            << "\t\t" << Form("%i s", fTimes[i]) << "\t\t" << fStart[i]
            << "\t\t" << fStop[i] << std::endl; 
 }
 std::cout << "\n" << std::endl;
}
//------------------------------------------------------------------
