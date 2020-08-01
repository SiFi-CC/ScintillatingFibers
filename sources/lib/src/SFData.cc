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
// #include <TROOT.h>

ClassImp(SFData);

//------------------------------------------------------------------
// constants
static const char  *gPath = getenv("SFDATA");  // path to the experimental data and data base
static const int    gBaselineMax = 50;         // number of samples for base line determination
static const double gmV          = 4.096;      // coefficient to calibrate ADC channels to mV
const double ampMax = 660;
//------------------------------------------------------------------
/// Default constructor. If this constructor is used the series 
/// number should be set via SetDetails(int seriesNo) function.
SFData::SFData(): fSeriesNo(-1),
                  fNpoints(-1),
                  fAnalysisGroup(-1),
                  fFiber("dummy"),
                  fFiberLength(-1),
                  fSource("dummy"),
                  fCollimator("dummy"),
                  fDesc("dummy"),
                  fTestBench("dummy"),
                  fSiPM("dummy"),
                  fOvervoltage(-1),
                  fCoupling("dummy"),
                  fTempFile("dummy") {
                      
 std::cout << "##### Warning in SFData constructor!" << std::endl;
 std::cout << "You are using the default constructor. Set the series number & open data base!" << std::endl;
}
//------------------------------------------------------------------ 
/// Standard constructor (recommended).
/// \param seriesNo is number of experimental series to analyze.
SFData::SFData(int seriesNo): fSeriesNo(seriesNo),
                              fNpoints(-1),
                              fAnalysisGroup(-1),
                              fFiber("dummy"),
                              fFiberLength(-1),
                              fSource("dummy"),
                              fCollimator("dummy"),
                              fDesc("dummy"),
                              fTestBench("dummy"),
                              fSiPM("dummy"),
                              fOvervoltage(-1),
                              fCoupling("dummy"),
                              fTempFile("dummy") {
                                  
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
/// \param name - name of the data base file. 
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
/// \param seriesNo - number of experimental series to analyze.
/// 
/// The following attributes of the experimental series are set within this 
/// function:
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
  ///- liber length [mm]
  ///- radioactive source type
  ///- test bench type
  ///- collimator type
  ///- SiPM type
  ///- overvoltage [V]
  ///- coupling type
  ///- number of measurements in the series
  ///- name of the measurement log file
  ///- name of the temperature log file
  ///- description of the series
  ///- analysis group number 
  query = Form("SELECT FIBER, FIBER_LENGTH, SOURCE, TEST_BENCH, COLLIMATOR, SIPM, OVERVOLTAGE, COUPLING, NO_MEASUREMENTS, LOG_FILE, TEMP_FILE, DESCRIPTION, ANALYSIS_GROUP FROM SERIES WHERE SERIES_ID = %i", fSeriesNo);
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *fiber = sqlite3_column_text(statement, 0);
    const unsigned char *source = sqlite3_column_text(statement, 2);
    const unsigned char *test_bench = sqlite3_column_text(statement, 3);
    const unsigned char *collimator = sqlite3_column_text(statement, 4);
    const unsigned char *sipm = sqlite3_column_text(statement, 5);
    const unsigned char *coupling = sqlite3_column_text(statement, 7);
    const unsigned char *logfile = sqlite3_column_text(statement, 9);
    const unsigned char *tempfile = sqlite3_column_text(statement, 10);
    const unsigned char *description = sqlite3_column_text(statement, 11);
    fNpoints = sqlite3_column_int(statement, 8);
    fFiberLength = sqlite3_column_double(statement, 1);
    fOvervoltage = sqlite3_column_double(statement, 6);
    fFiber = std::string(reinterpret_cast<const char*>(fiber));
    fSource = std::string(reinterpret_cast<const char*>(source));
    fDesc = std::string(reinterpret_cast<const char*>(description));
    fCollimator = std::string(reinterpret_cast<const char*>(collimator));
    fTestBench = std::string(reinterpret_cast<const char*>(test_bench));
    fSiPM = std::string(reinterpret_cast<const char*>(sipm));
    fCoupling = std::string(reinterpret_cast<const char*>(coupling));
    fLogFile = std::string(reinterpret_cast<const char*>(logfile));
    fTempFile = std::string(reinterpret_cast<const char*>(tempfile));
    fAnalysisGroup = sqlite3_column_int(statement, 12);
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
  ///- list of measurements IDs
  query = Form("SELECT MEASUREMENT_NAME, DURATION_TIME, SOURCE_POSITION, START_TIME, STOP_TIME, MEASUREMENT_ID FROM MEASUREMENT WHERE SERIES_ID = %i", fSeriesNo);
  status = sqlite3_prepare_v2(fDB, query, -1, &statement, nullptr);
  
  SFTools::CheckDBStatus(status, fDB);
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *name = sqlite3_column_text(statement, 0);
    fNames.push_back(std::string(reinterpret_cast<const char*>(name)));
    fTimes.push_back(sqlite3_column_int(statement, 1));
    fPositions.push_back(sqlite3_column_double(statement, 2));
    fStart.push_back(sqlite3_column_int(statement, 3));
    fStop.push_back(sqlite3_column_int(statement, 4));
    fMeasureID.push_back(sqlite3_column_int(statement, 5));
  }
  
  SFTools::CheckDBStatus(status, fDB);
  
  sqlite3_finalize(statement);
   
  return true;
}
//------------------------------------------------------------------
SFSignal* SFData::ConvertSignal(DDSignal *sig){
    
    SFSignal *converted = new SFSignal();
    converted->fAmp       = sig->GetAmplitude();
    converted->fCharge    = sig->GetCharge();
    converted->fPE        = sig->GetPE();
    converted->fT0        = sig->GetT0();
    converted->fTOT       = sig->GetTOT();
    converted->fSDDSignal = false;
    
    return converted;
}
//------------------------------------------------------------------
SFSignal* SFData::ConvertSignal(SDDSignal *sig){
    
    SFSignal *converted = new SFSignal();
    converted->fAmp       = sig->GetAmplitude();
    converted->fCharge    = sig->GetCharge();
    converted->fPE        = sig->GetPE();
    converted->fT0        = sig->GetT0();
    converted->fTOT       = sig->GetTOT();
    converted->fBL        = sig->GetBL();
    converted->fBLsig     = sig->GetBLSigma();
    converted->fPileUp    = sig->GetPileUp();
    converted->fVeto      = sig->GetVeto();
    converted->fSDDSignal = true;
    
    return converted;
}
//------------------------------------------------------------------
/// Parses given cut and checks if signal fulfills conditions specified by it. 
/// \param sig - currently analyzed signal, as read from the tree
/// \param cut - a logic cut to select specific signals.
///
/// For floating point variables (fAmp, fCharge, fPE, fT0, fTOT, fBL, fBL_sigma) 
/// the following syntax of the cut is acceptable: 
/// - inequality signs: '<' and '>'
/// - single expressions, e.g. "fAmp>50", "ch_0.fPE>10", "fT0<100"
/// - double expressions with &&, e.g. "fPE>10 && fPE<100", "ch_0.fT0>0 && ch_0.fT0<500" 
///
/// For intiger variables (fPileUp and fVeto) the following syntax is acceptable:
/// - equality sign '=='
/// - single expressions, eg. "fVeto==1", "fPileUp==0"
/// - double expressions, eg. "fVeto==0 && fPileUp==0"
///
/// If empty cut is passed returns true.
bool SFData::InterpretCut(SFSignal *sig, TString cut){
  
  bool result = true;
  if(cut=="" || cut==" ") return result;
 
  double amp    = sig->fAmp;
  double charge = sig->fCharge;
  double pe     = sig->fPE;
  double t0     = sig->fT0;
  double tot    = sig->fTOT;
  double bl     = sig->fBL;
  double blsig  = sig->fBLsig;
  int pileup    = sig->fPileUp;
  int veto      = sig->fVeto;
  bool SDDflag  = sig->fSDDSignal;  
  
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
    else if(letters_expr[ii]=='='){
      istop = ii+2;
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
    
   if(expression[i].find("fAmp")!=std::string::npos){   //cut on fAmp
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = amp<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
      logic[i] = amp>number[i]; 
     }
   }
   
   else if(expression[i].find("fPE")!=std::string::npos){   //cut on fPE
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = pe<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = pe>number[i]; 
     }
   }
   
   else if(expression[i].find("fCharge")!=std::string::npos){   //cut on fCharge
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = charge<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = charge>number[i]; 
     }
   }
   
   else if(expression[i].find("fT0")!=std::string::npos){   //cut on fT0
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = t0<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = t0>number[i]; 
     }
   } 
   
   else if(expression[i].find("fTOT")!=std::string::npos){   //cut on fTOT
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = tot<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = tot>number[i]; 
     }
   }  
   
   else if(SDDflag && expression[i].find("fBL")!=std::string::npos){   //cut on fBL
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = bl<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = bl>number[i]; 
     }
   }
   
   else if(SDDflag && expression[i].find("fBL_sigma")!=std::string::npos){   //cut on fBL_sigma
     if(expression[i].find("<")!=std::string::npos){
       logic[i] = blsig<number[i];
     }
     else if(expression[i].find(">")!=std::string::npos){
       logic[i] = blsig>number[i]; 
     }
   } 
   
   else if(SDDflag && expression[i].find("fPileUp")!=std::string::npos){   //cut on fPileUp
     if(expression[i].find("==")!=std::string::npos){
       logic[i] = pileup==number[i];
     }
   } 
   
   else if(SDDflag && expression[i].find("fVeto")!=std::string::npos){   //cut on fVeto
     if(expression[i].find("==")!=std::string::npos){
       logic[i] = veto==number[i];
     }
   } 
   
   else{
    std::cerr << "#### Error in SFData::InterpretCut! Incorrect cut syntax!" << std::endl;
    std::cerr << "Incorrect type. Available types are: fAmp, fCharge, fPE, fT0, fTOT, fBL, fBL_sigma, fPileUp and fVeto." << std::endl;
    return false;
   }
  }
  
  result = logic[0] && logic[1];
  
  return result;
}
//------------------------------------------------------------------
/// Accesses ROOT file and returns tree containing measured data for 
/// the requested measurement.
/// \param ID - measurement ID
SLoop* SFData::GetTree(int ID){
    
  int index = SFTools::GetIndex(fMeasureID, ID);
  std::string fname = std::string(SFTools::FindData(fNames[index]))+"/sifi_results.root";
  
  SLoop * loop = new SLoop();
  loop->addFile(fname);
  loop->setInput({});
  
  if(loop==nullptr){
    std::cerr << "#### Error in SFData::GetTree!" << std::endl;
    std::cerr << "Requested tree doesn't exist!" << std::endl;
    std::abort();
  }
  
  return loop;
}
//------------------------------------------------------------------
/// Returns single spectrum of requested type.
/// \param ch - chennel number
/// \param sel_type - type of the spectrum, as defined in SFDrawCommands class
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
/// \param ID - ID of requested measurements 
///
/// It is possible to have spectrum with cut or raw spectrum as recorded. In the latter case pass
/// empty string as cut.
TH1D* SFData::GetSpectrum(int ch, SFSelectionType sel_type, TString cut, int ID){

  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  TString fname = SFTools::FindData(fNames[index]);  
  TFile *file = new TFile(fname+"/sifi_results.root", "READ");
  TString tname = std::string("S");
  TTree *tree = (TTree*)file->Get(tname);
  
  gUnique+=1;
  TString selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch);

  tree->Draw(selection, cut);
//   TH1D *spec = (TH1D*)gROOT->FindObjectAny(Form("htemp%i", gUnique));
  TH1D *spec = (TH1D*)gDirectory->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_ch%i_pos%.1f_ID%i_", fSeriesNo, ch, position, ID) + SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  spec->SetName(hname);
  spec->SetTitle(htitle);
  
  return spec;
}
//------------------------------------------------------------------
/// Returns a vector with all spectra of requested type.
/// \param ch - channel number
/// \param sel_type - type of spectra (see SFDrawCommands)
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
///
/// Like with sigle spectrum, it is possible to have cut and raw spectra.
std::vector <TH1D*> SFData::GetSpectra(int ch, SFSelectionType sel_type, TString cut){

  std::vector <TH1D*> spectra;  
  for(int i=0; i<fNpoints; i++){
    spectra.push_back(GetSpectrum(ch, sel_type, cut, fMeasureID[i]));
  }
  
  return spectra;
}
//------------------------------------------------------------------
/// Returns single requested custom 1D histogram.
/// \param sel_type - predefined selection type (see SFDrawCommands)
/// \param cut - cut for drawn events. Also TTree-style syntax
/// \param ID - ID of requested measurement 
/// \param customNumbers - vector containing set of numbers necessary for 
/// the selection of events to be drawn on the histogram
TH1D* SFData::GetCustomHistogram(SFSelectionType sel_type, TString cut, int ID, 
                                std::vector <double> customNumbers){
  
  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  TString fname = SFTools::FindData(fNames[index]);
  TFile *file = new TFile(fname+"/sifi_results.root", "READ");
  TString tname = "S";
  TTree *tree = (TTree*)file->Get(tname);
  
  gUnique+=1;
  TString selection; 
  selection = SFDrawCommands::GetSelection(sel_type, gUnique, customNumbers);
  tree->Draw(selection, cut);
  TH1D* hist = (TH1D*)gDirectory->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_pos%.1f_ID%i_", fSeriesNo, position, ID) + SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  hist->SetName(hname);
  hist->SetTitle(htitle);
  
  return hist;
}
//------------------------------------------------------------------
/// Returns a vector of requested custom 1D histograms for all measurements in this series.
/// \param sel_type - predefined selection type (see SFDrawCommands)
/// \param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
/// all events will be drawn.
std::vector <TH1D*> SFData::GetCustomHistograms(SFSelectionType sel_type, TString cut){
  
  std::vector <TH1D*> hists;
  for(int i=0; i<fNpoints; i++){
    hists.push_back(GetCustomHistogram(sel_type, cut, fMeasureID[i]));
  }
  
  return hists;
}
//------------------------------------------------------------------
/// Returns single requested custom 1D histogram.
/// \param ch - channel number
/// \param sel_type - predefined selection type (see SFDrawCommands)
/// \param cut - cut for drawn events. Also TTree-style syntax
/// \param ID - ID of requested measurement
/// \param customNumbers - vector containing set of numbers necessary for 
/// the selection of events to be drawn on the histogram.
TH1D* SFData::GetCustomHistogram(int ch, SFSelectionType sel_type, TString cut, int ID, 
                                 std::vector <double> customNumbers){
    
  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  TString fname = SFTools::FindData(fNames[index]);
  TFile *file = new TFile(fname+"/sifi_results.root", "READ");
  TString tname = "S";
  TTree *tree = (TTree*)file->Get(tname);
  
  gUnique+=1;
  
  TString selection;
  if(customNumbers.empty())
    selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch);
  else 
    selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch, customNumbers);
  
  tree->Draw(selection, cut);
  TH1D* hist = (TH1D*)gDirectory->FindObjectAny(Form("htemp%i", gUnique));
  TString hname = Form("S%i_pos%.1f_ID%i_", fSeriesNo, position, ID)+ SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  hist->SetName(hname);
  hist->SetTitle(htitle);
  
  return hist;
    
}
//------------------------------------------------------------------
/// Returns single requested 2D correlation histogram.
/// \param sel_type - predefined selection type (see SFDrawCommands)
/// \param cut - cut for drawn events. Also TTree-style syntax
/// \param ID - ID of requested measurement
TH2D* SFData::GetCorrHistogram(SFSelectionType sel_type, TString cut, int ID, int ch){
  
  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  TString fname = SFTools::FindData(fNames[index]);
  TFile *file = new TFile(fname+"/sifi_results.root", "READ");
  TString tname = std::string("S");
  TTree *tree = (TTree*)file->Get(tname);
  
  gUnique+=1;
  
  TString selection;
  if(ch==-1)
      selection = SFDrawCommands::GetSelection(sel_type, gUnique);
  else 
      selection = SFDrawCommands::GetSelection(sel_type, gUnique, ch);
  
  tree->Draw(selection, cut, "colz");
  TH2D* hist = (TH2D*)gDirectory->FindObjectAny(Form("htemp%.i", gUnique));
  TString hname = Form("S%i_pos%.1f_ID%i_", fSeriesNo, position, ID) + SFDrawCommands::GetSelectionName(sel_type);
  TString htitle = hname + " " + cut;
  hist->SetName(hname);
  hist->SetTitle(htitle);
  
  return hist;
}
//------------------------------------------------------------------
/// Returns a vector of requested 2D correlation histograms for all measurements in 
/// this series.
/// \param sel_type - predefined selection type (see SFDrawCommands)
/// \param cut - cut for drawn events. Also TTree-style syntax.
std::vector <TH2D*> SFData::GetCorrHistograms(SFSelectionType sel_type, TString cut, int ch){
  
  std::vector <TH2D*> hists;
  for(int i=0; i<fNpoints; i++){
    hists.push_back(GetCorrHistogram(sel_type, cut, fMeasureID[i], ch));
  }
  
  return hists;
}
//------------------------------------------------------------------
TH2D *SFData::GetRefCorrHistogram(int ID, int ch){
    
    const double BL_sigma_cut = SFTools::GetSigmaBL(fSiPM);
    
    int index = SFTools::GetIndex(fMeasureID, ID);
    double position = fPositions[index];
    std::string fname = std::string(SFTools::FindData(fNames[index]))+"/sifi_results.root";
    
    SLoop * loop = new SLoop();
    loop->addFile(fname);
    loop->setInput({});
    SCategory *tSig = SCategoryManager::getCategory(SCategory::CatDDSamples);
    SCategory *tCal = SCategoryManager::getCategory(SCategory::CatFibersStackCal);
    
    TString hname = Form("S%i_pos%.1f_ID%i_PE%ivsPE2Correlation", fSeriesNo, position, ID, ch);
    TH2D *htemp = new TH2D(hname, hname, 1000, -100, 15E4, 2200, -150, 1500);
    
    int n = loop->getEntries();
    
    for(int i = 0; i < n; ++i){
        
        loop->nextEvent();
        size_t tentries = tSig->getEntries();
        
        assert(tSig->getEntries() == tCal->getEntries());

        uint coinc = 0;
        
        for (int j=0; j<tentries; ++j){

            int m, l, f;
            double mod0PE, mod1PE;
            SDDSamples *samples    = (SDDSamples *)tSig->getObject(j);
            SFibersStackCal *calib = (SFibersStackCal*)tCal->getObject(j); 
            samples->getAddress(m, l, f);
            
            if (ch == 0) {
                if(m == 0 &&
                    calib->getTimeL()>0 && calib->getQDCL()>0 &&
                    samples->getSignalL()->GetTOT()>0 &&
                    samples->getSignalL()->GetAmplitude()<ampMax &&
                    samples->getSignalL()->GetBLSigma() < BL_sigma_cut &&
                    calib->getTimeR()>0 && calib->getQDCR()>0 &&
                    samples->getSignalR()->GetTOT()>0 &&
                    samples->getSignalR()->GetAmplitude()<ampMax &&
                    samples->getSignalR()->GetBLSigma() < BL_sigma_cut)
                {
                    mod0PE = calib->getQDCL();
                    coinc |= 0x1;
                }
                else if(m == 1 &&
                    calib->getTimeL()>0 && calib->getQDCL()>0 &&
                    samples->getSignalL()->GetTOT()>0 &&
                    samples->getSignalL()->GetAmplitude()<ampMax &&
                    samples->getSignalL()->GetBLSigma() < BL_sigma_cut)
                  {
                      mod1PE = calib->getQDCL();
                      coinc |= 0x2;
                  }
                  if (coinc == 0x3)
                      htemp->Fill(mod1PE, mod0PE);
            }
            else if (ch == 1) {
                if(m == 0 &&
                    calib->getTimeL()>0 && calib->getQDCL()>0 &&
                    samples->getSignalL()->GetTOT()>0 &&
                    samples->getSignalL()->GetAmplitude()<ampMax &&
                    samples->getSignalL()->GetBLSigma() < BL_sigma_cut &&
                    calib->getTimeR()>0 && calib->getQDCR()>0 &&
                    samples->getSignalR()->GetTOT()>0 &&
                    samples->getSignalR()->GetAmplitude()<ampMax &&
                    samples->getSignalR()->GetBLSigma() < BL_sigma_cut)
                {
                    mod0PE = calib->getQDCR();
                    coinc |= 0x1;
                }
                  else if(m == 1 &&
                    calib->getTimeL()>0 && calib->getQDCL()>0 &&
                    samples->getSignalL()->GetTOT()>0 &&
                    samples->getSignalL()->GetAmplitude()<ampMax &&
                    samples->getSignalL()->GetBLSigma() < BL_sigma_cut)
                  {
                      mod1PE = calib->getQDCL();
                      coinc |= 0x2;
                  }
                if (coinc == 0x3)
                    htemp->Fill(mod1PE, mod0PE);
            }
        }
    }
    
    return htemp;
}
//------------------------------------------------------------------
std::vector <TH2D*> SFData::GetRefCorrHistograms(int ch){
    
    std::vector <TH2D*> hists;
    
    for(int i=0; i<fNpoints; i++){
        hists.push_back(GetRefCorrHistogram(fMeasureID[i], ch));
    }
    
    return hists;
}
//------------------------------------------------------------------
/// Returns averaged signal.
/// \param ch - channel number 
/// \param ID - ID of requested measurement
/// \param cut - logic cut to choose signals. Syntax of this cut is explained 
/// in InterpretCut() function
/// \param number - number of signals to be averaged
/// \param bl - flag for base line subtraction. If true - base line will be subtracted, 
/// if false - it won't
///
/// If no cut is required pass an empty string.
/// This function calls separate methods to access binary files depending on the test bench type.
TProfile* SFData::GetSignalAverage(int ch, int ID, TString cut, int number, bool bl){

  TProfile *sig = nullptr;
  
  if(fTestBench=="PL"){
    sig = GetSignalAverageKrakow(ch, ID, cut, number, bl);
  }
  else if(fTestBench=="DE"){
    sig = GetSignalAverageAachen(ch, ID, cut, number);    
  }
  else{
    std::cerr << "##### Error in SFData::GetSignalAverage()!" << std::endl;
    std::cerr << "Unknown data format!" << std::endl;
    std::abort();
  }

  return sig;
}
//------------------------------------------------------------------
/// This private function allows to access averaged signals recorded 
/// with the Krakow test bench. It opens binary file corresponding
/// to the chosen measurement and channel. Based on the digitized data
/// saved in the ROOT TTree it searches for signals which fulfil given cut
/// and performs averaging using ROOT's TProfile object. 
/// \param ch - channel number
/// \param ID - measuement ID
/// \param cut - logic cut to choose signals (syntax explained in InterpretCutt()
/// \param number - number of signals to be averaged
/// \param bl - flag for base line subtraction - if true baseline will be 
/// subtracted, if false - it will not.
TProfile* SFData::GetSignalAverageKrakow(int ch, int ID, TString cut, int number, bool bl){
 
  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  const int ipoints = 1024;
  float x;
    
  double BL_sigma_cut = SFTools::GetSigmaBL(fSiPM);
  
  std::string fname = std::string(SFTools::FindData(fNames[index]));
  SLoop * loop = new SLoop();
  loop->addFile(fname+"/sifi_results.root");
  loop->setInput({});
  SCategory *tSig = SCategoryManager::getCategory(SCategory::CatDDSamples); 
  SCategory *tCal = SCategoryManager::getCategory(SCategory::CatFibersStackCal);
  
  TString iname = fname + Form("/wave_%i.dat", ch);
  std::ifstream input(iname, std::ios::binary);

  if(! input.is_open()){ 
    std::cerr << "##### Error in SFData::GetSignalKraków()! Cannot open binary file!" << std::endl;
    std::cerr << iname << std::endl;
    std::abort();
  }
  
  TString hname = "sig_profile";
  TString htitle = "sig_profile";
  TProfile *psig = new TProfile(hname, htitle, ipoints, 0, ipoints, "");
  
  int nloop = loop->getEntries();
  float baseline = 0.;
  int counter = 0;
  int infile = 0;
  bool condition = true;
  float firstT0 = 0.;
  float PE, T0;
  
  for(int i = 0; i < nloop; ++i){
        
    //loop->nextEvent();
    loop->getEvent(i);
    size_t tentries = tSig->getEntries();
        
    assert(tSig->getEntries() == tCal->getEntries());
    
    for (int j=0; j<tentries; ++j){
        
      int m, l, f;
      SDDSamples *samples    = (SDDSamples *)tSig->getObject(j);
      SFibersStackCal *calib = (SFibersStackCal*)tCal->getObject(j);
      SDDSignal *sigL = (SDDSignal*)samples->getSignalL();
      SDDSignal *sigR = (SDDSignal*)samples->getSignalR();
      samples->getAddress(m, l, f);

      TProfile *hptr  = psig;
      SDDSignal *sptr = nullptr;
      
      if (ch == 0 && m == 0) {
          sptr = sigL;
          PE = calib->getQDCL();
          T0 = calib->getTimeL();
      }
      else if (ch == 1 && m == 0) {
          sptr = sigR;
          PE = calib->getQDCR();
          T0 = calib->getTimeR();
      }
      else if (ch == 2 && m == 1) {
          sptr = sigL;
          PE = calib->getQDCL();
          T0 = calib->getTimeL();
      }
      else {
          hptr = nullptr;
          PE = -100;
          T0 = -100;
      }

      if (hptr) {
        SFSignal *conv_sig = ConvertSignal(sptr); 
        conv_sig->fPE = PE;
        conv_sig->fT0 = T0;
        condition = InterpretCut(conv_sig, cut);
        if(condition && fabs(firstT0) < 1E-10 
           && conv_sig->fBLsig < BL_sigma_cut) 
            firstT0 = T0;
        if(condition && fabs(T0-firstT0) < 1
           && conv_sig->fBLsig< BL_sigma_cut){
          infile = sizeof(x)*ipoints*i;
          if(bl)  
            baseline = sptr->GetBL();
          input.seekg(infile);
          for(int ii=1; ii<ipoints+1; ii++){
            input.read(reinterpret_cast<char*>(&x), sizeof(float));
            if(bl) hptr->Fill(ii, (x-baseline)/gmV);
            else   hptr->Fill(ii, (x/gmV));
          }
          if(counter<number) 
            counter++;
          else 
              break;
          delete conv_sig;
        }
      }
    }
  }
  
  hname = Form("S%i_ch%i_pos_%.1f_ID%i_sig_num_%i", fSeriesNo, ch, position, ID, counter);
  htitle = hname + " " + cut;
  psig->SetName(hname);
  psig->SetTitle(htitle);
  
  if(counter<number){ 
    std::cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
              << " out of " << number << " plotted." << std::endl;
    std::cout << "Position: " << position << "\t channel: " << ch << std::endl; 
  }
  
  input.close();
  
  return psig;
}
//------------------------------------------------------------------
/// This private function allows to access averaged signals recorded 
/// with the Aachen test bench. It opens ROOT file and TTree corresponding
/// to the chosen measurement and channel. Based on the digitized data
/// saved in the ROOT TTree it searches for signals which fulfil given cut
/// and performs averaging using ROOT's TProfile object. 
/// \param ch - channel number
/// \param ID - measuement ID
/// \param cut - logic cut to choose signals (syntax explained in InterpretCutt()
/// \param number - number of signals to be averaged.
TProfile* SFData::GetSignalAverageAachen(int ch, int ID, TString cut, int number){
  
  int index = SFTools::GetIndex(fMeasureID, ID);
  double position = fPositions[index];
  const int ipoints = 1024;
  
  TString fname = SFTools::FindData(fNames[index]);
  TFile *file = new TFile(fname+"/results.root", "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = fname + "/waves.root";
  TFile* iFile = new TFile(iname, "READ");
  TTree* iTree = (TTree*)iFile->Get("wavetree");
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  TString bname = Form("voltages_ch_%i", ch);
  iTree->SetBranchAddress(bname, &iVolt);  
  
  TString hname = "sig_profile";
  TString htitle = "sig_profile";
  TProfile *psig = new TProfile(hname, htitle, ipoints, 0, ipoints, "");
  
  int nentries = tree->GetEntries();
  int counter = 0;
  bool condition = true;
  double firstT0 = 0.;
  
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   SFSignal *conv_sig = ConvertSignal(sig);  
   condition = InterpretCut(conv_sig, cut);
   if(condition && fabs(firstT0)<1E-10) firstT0 = sig->GetT0();
   if(condition && fabs(sig->GetT0()-firstT0)<1){
     iTree->GetEntry(i);
     for(int ii=0; ii<ipoints; ii++){
       psig->Fill(ii+1, (*iVolt)[ii]);
     }
     if(counter<number) counter++;
     else break;
    }
  }
  
  hname = Form("S%i_ch%i_pos_%.1f_ID%i_sig_num_%i", fSeriesNo, ch, position, ID, counter);
  htitle = hname + " " + cut;
  psig->SetName(hname);
  psig->SetTitle(htitle);
  
  if(counter<number){ 
    std::cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
              << " out of " << number << " plotted." << std::endl;
    std::cout << "Position: " << position << "\t channel: " << ch << std::endl; 
  }

  return psig;
}
//------------------------------------------------------------------
/// Returns single raw signal.
/// \param ch - channel number 
/// \param ID - ID of requested measurement
/// \param cut - logic cut to choose signals (syntax explained in InterpretCut() function)
/// \param number - number of the signal to be drawn, eg. number=1 means that first signal 
/// which fulfills given cut will be drawn
/// \param bl - flag for base line subtraction. If true - base line will be subtracted, 
/// if false - it won't
///
/// If no cut is required pass an empty string.
/// This function calls separate methods to access binary files depending on the test bench type.
TH1D* SFData::GetSignal(int ch, int ID, TString cut, int number, bool bl){
 
  TH1D *sig = nullptr;
  
  if(fTestBench=="PL"){
    sig = GetSignalKrakow(ch, ID, cut, number, bl);   
  }
  else if(fTestBench=="DE"){
    sig = GetSignalAachen(ch, ID, cut, number);
  }
  else{
    std::cerr << "##### Error in SFData::GetSignal()" << std::endl;
    std::cerr << "Unknown data format!" << std::endl;
    std::abort();
  }

  return sig;
}
//------------------------------------------------------------------
/// This private function allows to access single raw signals recorded 
/// with the Krakow test bench. It opens binary file corresponding
/// to the chosen measurement and channel. Based on the digitized data
/// saved in the ROOT TTree it searches for signals which fulfil given cut
/// and draws them. 
/// \param ch - channel number
/// \param ID - measuement ID
/// \param cut - logic cut to choose signals (syntax explained in InterpretCutt()
/// \param number - number of signals to be averaged
/// \param bl - flag for base line subtraction - if true baseline will be 
/// subtracted, if false - it will not.
TH1D* SFData::GetSignalKrakow(int ch, int ID, TString cut, int number, bool bl){

  int index = SFTools::GetIndex(fMeasureID, ID);
  const int ipoints = 1024;
  float x; 
  
  double BL_sigma_cut = SFTools::GetSigmaBL(fSiPM);
  
  std::string fname = std::string(SFTools::FindData(fNames[index]));
  double position = fPositions[index];
  
  SLoop * loop = new SLoop();
  loop->addFile(fname+"/sifi_results.root");
  loop->setInput({});
  SCategory *tSig = SCategoryManager::getCategory(SCategory::CatDDSamples);
  SCategory *tCal = SCategoryManager::getCategory(SCategory::CatFibersStackCal);
  
  TString iname = fname + Form("/wave_%i.dat", ch);
  std::ifstream input(iname, std::ios::binary);

  if(! input.is_open()){ 
    std::cerr << "##### Error in SFData::GetSignalKraków()! Cannot open binary file!" << std::endl;
    std::cerr << iname << std::endl;
    std::abort();
  }

  TString hname = Form("S%i_ch%i_pos_%.1f_ID%i_sig_no%i", fSeriesNo, ch, position, ID, number);
  TString htitle = hname + " " + cut; 
  
  TH1D *hsig = new TH1D(hname, htitle, ipoints, 0, ipoints);
  
  int nloop = loop->getEntries();
  double baseline = 0.;
  int infile = 0;
  int counter = 0;
  bool condition = true;
  double PE, T0;
  
  for(int i = 0; i < nloop; ++i){
        
    //loop->nextEvent();
    loop->getEvent(i);
    size_t tentries = tSig->getEntries();
    
    assert(tSig->getEntries() == tCal->getEntries());
    
    for (int j=0; j<tentries; ++j){
        
      int m, l, f;
      SDDSamples *samples    = (SDDSamples *)tSig->getObject(j);
      SFibersStackCal *calib = (SFibersStackCal*)tCal->getObject(j);
      SDDSignal *sigL = (SDDSignal*)samples->getSignalL();
      SDDSignal *sigR = (SDDSignal*)samples->getSignalR();
      samples->getAddress(m, l, f);

      TH1 * hptr = hsig;
      SDDSignal * sptr = nullptr;
      
      if (ch == 0 && m == 0) {
          sptr = sigL;
          PE = calib->getQDCL();
          T0 = calib->getTimeL();
      }
      else if (ch == 1 && m == 0) {
          sptr = sigR;
          PE = calib->getQDCR();
          T0 = calib->getTimeR();
      }
      else if (ch == 2 && m == 1) {
          sptr = sigL;
          PE = calib->getQDCL();
          T0 = calib->getTimeL();
      }
      else {
          hptr = nullptr;
          PE = -100;
          T0 = -100;
      }

      if (hptr) {
        SFSignal *conv_sig = ConvertSignal(sptr); 
        conv_sig->fPE = PE;
        conv_sig->fT0 = T0;
        condition = InterpretCut(conv_sig, cut);
        if(condition && conv_sig->fBLsig < BL_sigma_cut){
          counter++;
          if(counter!=number) continue;
          infile = sizeof(x)*ipoints*i;
          if(bl) baseline = sptr->GetBL();
          input.seekg(infile);
          for(int ii=1; ii<ipoints+1; ii++){
            input.read(reinterpret_cast<char*>(&x), sizeof(float));
            if(bl) hptr->SetBinContent(ii, (x-baseline)/gmV);
            else   hptr->SetBinContent(ii, (x/gmV));
          }
        }
      }
    }
  }
      
  input.close();
  
  return hsig;
}
//------------------------------------------------------------------
/// This private function allows to access raw signals recorded 
/// with the Aachen test bench. It opens ROOT file and TTree corresponding
/// to the chosen measurement and channel. Based on the digitized data
/// saved in the ROOT TTree it searches for signals which fulfil given cut
/// and draws them. 
/// Returns single signal.
/// \param ch - channel number
/// \param ID - measurement ID
/// \param cut - logic cut to choose signals. Syntax of this cut is explained in InterpretCut() function
/// \param number - requested number of the signal to be drawn.
TH1D* SFData::GetSignalAachen(int ch, int ID, TString cut, int number){
 
  int index = SFTools::GetIndex(fMeasureID, ID);
  const int ipoints = 1024;
  
  TString fname = SFTools::FindData(fNames[index]);
  double position = fPositions[index];
  TFile *file = new TFile(fname+"/results.root", "READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i", ch), &sig);
  
  TString iname = fname + "/waves.root";
  TFile* iFile = new TFile(iname, "READ");
  TTree* iTree = (TTree*)iFile->Get("wavetree");
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  TString bname = Form("voltages_ch_%i", ch);
  iTree->SetBranchAddress(bname, &iVolt);
  
  TString hname = Form("S%i_ch%i_pos_%.1f_ID%i_sig_no%i", fSeriesNo, ch, position, ID, number);
  TString htitle = hname + " " + cut;
  
  TH1D *hsig = new TH1D(hname, htitle, ipoints, 0, ipoints);
  
  int nentries = tree->GetEntries();
  int counter = 0;
  bool condition = true;
  
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    iTree->GetEntry(i);
    SFSignal *conv_sig = ConvertSignal(sig);  
    condition = InterpretCut(conv_sig, cut);
    if(condition){
      counter++;
      if(counter!=number) continue;
      for(int ii=0; ii<ipoints; ii++){
        hsig->SetBinContent(ii+1, (*iVolt)[ii]);
      }
    }
  }
  
  return hsig;
}
//------------------------------------------------------------------
/// Prints details of currently analyzed experimental series.
void SFData::Print(void){
 std::cout << "\n\n------------------------------------------------" << std::endl;
 std::cout << "This is Print() for SFData class object" << std::endl;
 std::cout << "Number of the experimental series: " << fSeriesNo << std::endl;
 std::cout << fDesc << std::endl;
 std::cout << "Collimator: " << fCollimator << std::endl;
 std::cout << "Test bench: " << fTestBench << std::endl;
 std::cout << "Number of measurements in this series: " << fNpoints << std::endl;
 std::cout << "Fiber: " << fFiber << std::endl;
 std::cout << "Fiber length: " << fFiberLength << " mm" << std::endl;
 std::cout << "Coupling: " << fCoupling << std::endl;
 std::cout << "Radioactive source: " << fSource << std::endl;
 std::cout << "SiPM: " << fSiPM <<std::endl;
 std::cout << "Overvoltage: " << fOvervoltage << " V" << std::endl;
 std::cout << "Temperature logfile: " << fTempFile << std::endl;
 std::cout << "List of measurements in this series:" << std::endl;
 for(int i=0; i<fNpoints; i++){
  std::cout << std::setw(30);
  std::cout << fNames[i] << "\t\t" << Form("%.1f mm", fPositions[i]) 
            << "\t\t" << Form("%i s", fTimes[i]) << "\t\t" << fStart[i]
            << "\t\t" << fStop[i] << "\t\t" << fMeasureID[i] << std::endl; 
 }
 std::cout << "\n" << std::endl;
}
//------------------------------------------------------------------
