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
char  *gPath        = getenv("SFDATA");
double gUnique      = 0.;
int    gBaselineMax = 50;
double gmV          = 4.096;
//------------------------------------------------------------------
///Default constructor. If this constructor is used the series 
///number should be set via SetDetails(int seriesNo) function.
SFData::SFData(){
 Reset();
 cout << "##### Warning in SFData constructor!" << endl;
 cout << "You are using the default constructor. Set the series number & open data base!" <<endl;
}
//------------------------------------------------------------------
///Standard constructor (recommended).
///\param seriesNo is number of experimental series to analyze.
///
///By deafult data analyzed with fixed threshold in DD6 is accessed.
///If you need constant fraction data use SetThreshold() function.
SFData::SFData(int seriesNo){
 bool db_stat  = OpenDataBase("ScintFib");
 bool set_stat = SetDetails(seriesNo);
 if(!db_stat || !set_stat){
   throw "##### Exception in SFData constructor!";
 }
}
//------------------------------------------------------------------
///Standard constructor.
///\param seriesNo is number of the experimental series to analyze.
///\param threshold is threshold type in DD6 preliminary data analysis.
///Possible options are: "ft" - fixed threshold and "cf" - constant fraction.
SFData::SFData(int seriesNo, TString threshold){
 bool db_stat  = OpenDataBase("ScintFib");
 bool set_stat = SetDetails(seriesNo);
 bool thr_stat = SetThreshold(threshold);
 if(!db_stat || !set_stat || !thr_stat){
   throw "##### Exception in SFData constructor!";
 }
}
//------------------------------------------------------------------
///Default destructor.
SFData::~SFData(){
 int status = sqlite3_close(fDB);
 if(status!=0) 
   cout << "In SFData destructor. Data base corrupted!" << endl;
 //else 
 //  cout << "In SFData destructor. Data base clossed succesfully!" << endl;
}
//------------------------------------------------------------------
///Opens SQLite3 data base containing details of experimental series
///and measurements.
bool SFData::OpenDataBase(TString name){

 TString db_name = string(gPath)+"/DB/"+name;
 int status = sqlite3_open(db_name,&fDB);
 
 if(status!=0){
   cout << "##### Error in SFData::OpenDataBase()!" << endl;
   cout << "Could not access data base!" << endl;
   return false;
 }
 //else
 //  cout << "Data base opened succesfully!" << endl;
 
 return true;
}
//------------------------------------------------------------------
///Sets all details of selected experimental series. If default constructor
///was used, this function needs to be called explicitly with the number of 
///of requested series as an argument. The following attributes are set 
///within this function:
bool SFData::SetDetails(int seriesNo){
  
  TString query;
  sqlite3_stmt *statement;
  int status;
  
  Reset();
  fSeriesNo = seriesNo;
  fThreshold = "ft";	//default - fixed threshold
  
  //-----Checking if series number is valid
  int maxSeries;
  query = "SELECT COUNT(*) FROM SERIES";
  status = sqlite3_prepare_v2(fDB,query,-1,&statement,NULL);
  
  if(status!=SQLITE_OK){
    cout << "##### SQL Error: " <<  sqlite3_errmsg(fDB) << endl;
  }
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    maxSeries = sqlite3_column_int(statement,0);
  }
  
  if(status!=SQLITE_DONE){
    cout << "##### SQL Error: " << sqlite3_errmsg(fDB) << endl;
  }
  
  sqlite3_finalize(statement);

  if(fSeriesNo<1 || fSeriesNo>maxSeries){
   cout << "##### Error in SFData::SetDetails()! Series number out of range!" << endl;
   return false;
  }
  //-----
  
  //----- Setting series attributes
  ///- fiber type 
  ///- measurement type 
  ///- radioactive source type
  ///- number of measurements in the series
  ///- description of the series
  query = Form("SELECT FIBER, SOURCE, NO_MEASUREMENTS, DESCRIPTION, MEASURE_TYPE FROM SERIES WHERE SERIES_NO = %i",fSeriesNo);
  status = sqlite3_prepare_v2(fDB,query,-1,&statement,NULL);
  
  if(status!=SQLITE_OK){
    cout << "##### SQL Error: " <<  sqlite3_errmsg(fDB) << endl;
    return false;
  }
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *fiber = sqlite3_column_text(statement,0);
    const unsigned char *source = sqlite3_column_text(statement,1);
    const unsigned char *description = sqlite3_column_text(statement,3);
    const unsigned char *measure_type = sqlite3_column_text(statement,4);
    fNpoints = sqlite3_column_int(statement,2);
    fFiber = string(reinterpret_cast<const char*>(fiber));
    fSource = string(reinterpret_cast<const char*>(source));
    fDesc = string(reinterpret_cast<const char*>(description));
    fType = string(reinterpret_cast<const char*>(measure_type));
  }
  
  if(status!=SQLITE_DONE){
    cout << "##### SQL Error: " << sqlite3_errmsg(fDB) << endl;
    return false;
  }
  
  sqlite3_finalize(statement);
  //-----
  
  //----- Setting measurements attributes
  ///- list of measurements names
  ///- list of source positions
  ///- list of measurements times
  query = Form("SELECT NAME, POSITION, TIME FROM MEASUREMENTS WHERE SERIES_NO = %i",fSeriesNo);
  status = sqlite3_prepare_v2(fDB,query,-1,&statement,NULL);
  
  if(status!=SQLITE_OK){
    cout << "##### SQL Error: " <<  sqlite3_errmsg(fDB) << endl;
    return false;
  }
  
  while((status=sqlite3_step(statement)) == SQLITE_ROW){
    const unsigned char *name = sqlite3_column_text(statement,0);
    fNames.push_back(string(reinterpret_cast<const char*>(name)));
    fPositions.push_back(sqlite3_column_double(statement,1));
    fTimes.push_back(sqlite3_column_double(statement,2));
  }
  
  if(status!=SQLITE_DONE){
    cout << "##### SQL Error: " << sqlite3_errmsg(fDB) << endl;
    return false;
  }
  
  sqlite3_finalize(statement);
   
  return true;
}
//------------------------------------------------------------------
/// Sets the flag to identify the tree with data which should be accessed.
/// \param threshold - type of threshold used during data analysis in DD6.
/// Possible options are: ft - fixed threshold (default) and cf - constant
/// fraction. 
bool SFData::SetThreshold(TString threshold){
  if(fType=="Lead"){
  	if(!(threshold=="ft" || threshold=="cf")){
    	cout << "##### Error in SFData::SetThreshold()! Incorrect threshold type!" << endl;
    	cout << "Possible options are: ft for fixed threshold and cf for constant fraction" << endl;
    	return false;
  	}
  	fThreshold = threshold;
  	return true;
  }
  else{
    	if(!(threshold=="ft")){
		cout << "##### Error in SFData::SetThreshold()! Set threshold not possible!" << endl;
    		cout << "For the measurements with the electric collimator is it not possible to set the threshold type" << endl;
		return false;
	}
	else return true;	
  }
}
//------------------------------------------------------------------
/// Returns proper selection for Draw() method of TTree. It sets binning, 
/// ranges and unique names for created histograms.
/// \param ch - channel number
/// \param type - type of spectrum to be drawn. Possible options are:
/// fAmp, fCharge, fPE, fT0 and fTOT
TString SFData::GetSelection(int ch, TString type){
 
  TString selection;
  gUnique = gRandom->Uniform(0,1);
  
  if(type=="fAmp")
    selection = Form("ch_%i.fAmp>>htemp%.7f(1000,0,700)",ch,gUnique);
  else if(type=="fCharge")
    selection = Form("ch_%i.fCharge>>htemp%.7f(1000,-1E4,2.5E5)",ch,gUnique);
  else if(type=="fPE")
    selection = Form("ch_%i.fPE>>htemp%.7f(1000,-150,1200)",ch,gUnique);
  else if(type=="fT0")
    selection = Form("ch_%i.fT0>>htemp%.7f(1000,-110,1100)",ch,gUnique);
  else if(type=="fTOT")
    selection = Form("ch_%i.fTOT>>htemp%.7f(1000,-110,1100)",ch,gUnique);
  else{
    cout << "##### Error in SFData::GetSelection()! Incorrect type!" << endl;
    cout << "Possible options are: fAmp, fCharge, fPE, fT0, fTOT" << endl;
    return "";
  }
    
  return selection;
}
//------------------------------------------------------------------
/// Returns selection for Draw method of TTree for custom histograms. 
/// It sets binning, ranges and unique names for created histograms. 
/// \param selection is custom selection entered by user.
TString SFData::GetSelectionCustom(TString selection){
 
  TString selectAndDraw;
  gUnique = gRandom->Uniform(0,1);
  
  if(selection=="log(sqrt(ch_1.fPE/ch_0.fPE))")
    selectAndDraw = selection+Form(">>htemp%.7f(500,-5,5)",gUnique);
  else if(selection=="ch_0.fT0-ch_1.fT0")
    selectAndDraw = selection+Form(">>htemp%.7f(500,-50,50)",gUnique);
  else if(selection=="sqrt(ch_0.fPE*ch_1.fPE)")
    selectAndDraw = selection+Form(">>htemp%.7f(1000,-150,1200)",gUnique);
  else if(selection=="sqrt(ch_0.fAmp*ch_1.fAmp)")
    selectAndDraw = selection+Form(">>htemp%.7f(1000,0,700)",gUnique);
  else if(selection=="ch_0.fPE:ch_1.fPE")
    selectAndDraw = selection+Form(">>htemp%.7f(1000,-150,1200,1000,-150,1200)",gUnique);
  else if(selection=="ch_0.fAmp:ch_1.fAmp")
    selectAndDraw = selection+Form(">>htemp%.7f(1000,0,700,1000,0,700)",gUnique);
  else if(selection=="ch_0.fT0:ch_1.fT0")
    selectAndDraw = selection+Form(">>htemp%.7f(1000,-110,1100,1000,-110,1100)",gUnique);
  else{
    cout << "##### Warning in SFData::GetSelectionCustom()!" << endl;
    cout << "Unknown selection! Deafault selection used!" << endl;
    selectAndDraw = selection+Form(">>htemp%.7f",gUnique);
  }
  
  return selectAndDraw;
}
//------------------------------------------------------------------
/// Returns index in the fNames and fPositions arrays for the 
/// measurement of requested source position in mm. 
/// If measurements in analyzed series don't have unique positions 
/// a number of measurement should be passed. Measurements counting 
/// starts at 1.
int SFData::GetIndex(double position){
   
  int index = -1;

  if(!fDesc.Contains("Regular series")){
    index = position-1;
    return index;
  }

  for(int i=0; i<fNpoints; i++){
    if(fabs(fPositions[i]-position)<3){
      index = i;
      break;
    }
  }

  if(index==-1){
   cout << "##### Error in SFData::GetIndex()! Incorrect position!" << endl;
   return index; 
  }
  
  return index;
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
  string cut_str = string(cut);
  int nletters = cut_str.length();
  char letters[nletters];
  strcpy(letters,cut_str.c_str());
  
  //splitting cut into expressions
  int iposition = -1;
  int nexpressions = 0;
  string expression[2];
  
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
    expression[0] = string(&letters[0], &letters[iposition]);
    expression[1] = string(&letters[iposition+2], &letters[nletters]);
  }
    
  //extracting doubles 
  int nletters_expr = 0;
  int istop = -1;
  char letters_expr[100];
  string number_str[2];
  double number[2];
  
  for(int i=0; i<nexpressions; i++){
   istop = -1;
   nletters_expr = expression[i].length();
   strcpy(letters_expr,expression[i].c_str());
   for(int ii=nletters_expr; ii>0; ii--){
    if(letters_expr[ii]=='<' || letters_expr[ii]=='>'){
      istop = ii+1;
      break;
    }
   }
   if(istop==-1){
     cout << "#### Error in SFData::InterpretCut! Incorrect cut syntax!" << endl;
     cout << "Missing '<' or '>'." << endl;
     return false;
   }
   number_str[i] = string(&letters_expr[istop], &letters_expr[nletters_expr]);
   number[i] = atof(number_str[i].c_str());
  }

  //checking logic
  bool logic[2] = {true,true};
  
  for(int i=0; i<nexpressions; i++){
    
   if(expression[i].find("fAmp")!=string::npos){		//cut on fAmp
     if(expression[i].find("<")!=string::npos){
       logic[i] = amp<number[i];
     }
     else if(expression[i].find(">")!=string::npos){
      logic[i] = amp>number[i]; 
     }
   }
   
   else if(expression[i].find("fPE")!=string::npos){		//cut on fPE
     if(expression[i].find("<")!=string::npos){
       logic[i] = pe<number[i];
     }
     else if(expression[i].find(">")!=string::npos){
       logic[i] = pe>number[i]; 
     }
   }
   
   else if(expression[i].find("fCharge")!=string::npos){	//cut on fCharge
     if(expression[i].find("<")!=string::npos){
       logic[i] = charge<number[i];
     }
     else if(expression[i].find(">")!=string::npos){
       logic[i] = charge>number[i]; 
     }
   }
   
   else if(expression[i].find("fT0")!=string::npos){		//cut on fT0
     if(expression[i].find("<")!=string::npos){
       logic[i] = t0<number[i];
     }
     else if(expression[i].find(">")!=string::npos){
       logic[i] = t0>number[i]; 
     }
   } 
   
   else if(expression[i].find("fTOT")!=string::npos){		//cut on fTOT
     if(expression[i].find("<")!=string::npos){
       logic[i] = tot<number[i];
     }
     else if(expression[i].find(">")!=string::npos){
       logic[i] = tot>number[i]; 
     }
   } 
   
   else{
    cout << "#### Error in SFData::InterpretCut! Incorrect cut syntax!" << endl;
    cout << "Incorrect type. Available types are: fAmp, fCharge, fPE, fT0 and fTOT." << endl;
    return false;
   }
  }
  
  result = logic[0] && logic[1];
  
  return result;
}
//------------------------------------------------------------------
/// Returns single spectrum of requested type.
/// \param ch - chennel number
/// \param type - type of the spectrum. Possible options are: fAmp, fCharge, fPE, fT0 and fTOT
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
/// \param position - position of the source in mm. If analyzed series doesn't have
/// unique positions a number of measurement should be passed here. Numbering starts at 1. 
///
/// It is possible to have spectrum with cut or raw spectrum as recorded. In the latter case pass
/// empty string as cut.
TH1D* SFData::GetSpectrum(int ch, TString type, TString cut, double position){

  int index = GetIndex(position);
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  fSpectrum = new TH1D();
  
  TString selection = GetSelection(ch,type);
  tree->Draw(selection,cut);
  fSpectrum = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
  TString hname = Form("S%i_ch%i_pos%.1f_",fSeriesNo,ch,position)+type+string("_")+fThreshold;
  TString htitle = hname+" "+cut;
  fSpectrum->SetName(hname);
  fSpectrum->SetTitle(htitle);
  
  return fSpectrum; 
}
//------------------------------------------------------------------
/// Returns a vector with all spectra of requested type.
/// \param ch - channel number
/// \param type - type of spectra (fAmp, fCharge, fPE, fT0, fTOT)
/// \param cut - logic cut for drawn events (syntax like for Draw() method of TTree)
///
/// Like with sigle spectrum, it is possible to have cut and raw spectra.
vector <TH1D*> SFData::GetSpectra(int ch, TString type, TString cut){

  bool empty = fSpectra.empty();
  for(int i=0; i<fNpoints; i++){
   if(empty) fSpectra.push_back(new TH1D());
   if(fDesc.Contains("Regular series"))
     fSpectra[i] = GetSpectrum(ch,type,cut,fPositions[i]);
   else 
     fSpectra[i] = GetSpectrum(ch,type,cut,i+1);
  }
  
  return fSpectra;
}
//------------------------------------------------------------------
///Returns single requested custom 1D histogram.
///\param selection - selection like for Draw() method of TTree.
///\param cut - cut for drawn events. Also TTree-style syntax.
///\param position - position of source in mm. If position is not unique a 
///number of measurement should be entered.
TH1D* SFData::GetCustomHistogram(TString selection, TString cut, double position){
  
  int index = GetIndex(position);
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  fHist = new TH1D();
  
  TString selectAndDraw = GetSelectionCustom(selection);
  tree->Draw(selectAndDraw,cut);
  fHist = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
  TString hname = Form("S%i_pos%.1f_",fSeriesNo,position)+selection+string("_")+fThreshold;
  TString htitle = hname+" "+cut;
  fHist->SetName(hname);
  fHist->SetTitle(htitle);
  
  return fHist;
}
//------------------------------------------------------------------
///Returns a vector of requested custom 1D histograms for all measurements in this series.
///\param selection - selection like for Draw() method of TTree, e.g. "log(ch_0.fPE/ch_1.fPE)". 
///Redirection to specific histogram should not be entered here, it is done inside the function.
///\param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
///all events will be drawn.
vector <TH1D*> SFData::GetCustomHistograms(TString selection, TString cut){
  
  bool empty = fHists.empty();
  for(int i=0; i<fNpoints; i++){
    if(empty) fHists.push_back(new TH1D());
    if(fDesc.Contains("Regular series"))
      fHists[i] = GetCustomHistogram(selection,cut,fPositions[i]);
    else 
      fHists[i] = GetCustomHistogram(selection,cut,i+1);
  }
  
  return fHists;
}
//------------------------------------------------------------------
///Returns single requested correlation 2D histogram.
///\param selection - selection like for Draw() method of TTree.
///\param cut - cut for drawn events. Also TTree-style syntax.
///\param position - position of source in mm. If position is not unique a 
///number of measurement should be entered.
TH2D* SFData::GetCorrHistogram(TString selection, TString cut, double position){
  
  int index = GetIndex(position);
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  fHist2D = new TH2D();
  
  TString selectAndDraw = GetSelectionCustom(selection);
  tree->Draw(selectAndDraw,cut,"colz");
  fHist2D = (TH2D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
  TString hname = selection+Form("S%i_pos%.1f_",fSeriesNo,position)+selection+string("_")+fThreshold;
  TString htitle = hname+" "+cut;
  fHist2D->SetName(hname);
  fHist2D->SetTitle(htitle);
  
  return fHist2D;
}
//------------------------------------------------------------------
///Returns a vector of requested 2D correlation histograms for all measurements in 
///this series.
///\param selection - selection like for Draw() method of TTree, e.g. "ch_0.fPE:ch_1.fPE". 
///Redirection to specific histogram should not be entered here, it is done inside the function.
///\param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
///all events will be drawn.
vector <TH2D*> SFData::GetCorrHistograms(TString selection, TString cut){
  
  bool empty = fHists2D.empty();
  for(int i=0; i<fNpoints; i++){
    if(empty) fHists2D.push_back(new TH2D());
    if(fDesc.Contains("Regular series"))
      fHists2D[i] = GetCorrHistogram(selection,cut,fPositions[i]);
    else 
      fHists2D[i] = GetCorrHistogram(selection,cut,i+1);
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
  
  int index = GetIndex(position);
  const int ipoints = 1024;
  float x;
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i",ch),&sig);
  
  TFile* iFile;
  TTree* iTree;
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  TString iname; 
  ifstream input;
  if(fType.Contains("Lead")){
  	iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
  	input.open(iname,ios::binary);
  }
  else if(fType.Contains("Electric")){
	iname = string(gPath)+fNames[index]+"/waves.root";
	iFile = new TFile(iname,"READ");
	iTree = (TTree*)iFile->Get("wavetree");
	iname = Form("voltages_ch_%i",ch);
	iTree->SetBranchAddress(iname,&iVolt);
  }
  TString hname = "sig_profile";
  TString htitle = "sig_profile";
  fSignalProfile = new TProfile(hname,htitle,ipoints,0,ipoints,"");
    
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int counter = 0;
  int infile = 0;
  bool condition = true;
  double firstT0 = 0.;
  
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   condition = InterpretCut(sig,cut);
   if(condition && fabs(firstT0)<1E-10) firstT0 = sig->GetT0();
   if(condition && fabs(sig->GetT0()-firstT0)<1){
     if(fType.Contains("Lead")){
     	infile = sizeof(x)*ipoints*i;
     	if(bl){
     	  input.seekg(infile);
     	  baseline = 0.;
     	  for(int ii=0; ii<gBaselineMax; ii++){
     	    input.read((char*)&x,sizeof(x));
     	    baseline += x/gmV;
     	  }
     	  baseline = baseline/gBaselineMax;
     	}
     	input.seekg(infile);
     	for(int ii=0; ii<ipoints; ii++){
     	  input.read((char*)&x,sizeof(x));
     	  if(bl) fSignalProfile->Fill(ii,(x/gmV)-baseline);
     	  else   fSignalProfile->Fill(ii,(x/gmV));
     	}
     }
     else if(fType.Contains("Electric")){
	iTree->GetEntry(i);
	for(int ii=0; ii<ipoints; ii++){
		  fSignalProfile->Fill(ii+1,(*iVolt)[ii]);
	}
     }
     if(counter<number) counter++;
     else break;
    }
  }
  
  hname = Form("S%i_ch%i_pos_%.1f_sig_num_%i_",fSeriesNo,ch,position,counter)+fThreshold;
  htitle = hname +" "+cut;
  fSignalProfile->SetName(hname);
  fSignalProfile->SetTitle(htitle);
  
  if(counter<number){ 
    cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
         << " out of " << number << " plotted." << endl;
    cout << "Position: " << position << "\t channel: " << ch << endl; 
  }
  
  input.close();
  
  return fSignalProfile;
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
TH1D* SFData::GetSignal(int ch, double position, TString cut, int number, bool bl){
 
  int index = GetIndex(position);
  const int ipoints = 1024;
  float x;
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i",ch),&sig);
  
  TString iname; 
  ifstream input;
  TFile* iFile;
  TTree* iTree;
  TVectorT<float>* iVolt= new TVectorT<float>(ipoints);
  
  if(fType.Contains("Lead")){
	iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
	input.open(iname,ios::binary);
  }
  else if(fType.Contains("Electric")){
	iname = string(gPath)+fNames[index]+"/waves.root";
	iFile = new TFile(iname,"READ");
	iTree = (TTree*)iFile->Get("wavetree");
	iname = Form("voltages_ch_%i",ch);
	iTree->SetBranchAddress(iname,&iVolt);
  }

  TString hname = Form("S%i_ch%i_pos_%.1f_sig_no%i_",fSeriesNo,ch,position,number)+fThreshold;
  TString htitle = hname+" "+cut;
  
  fSignal = new TH1D(hname,htitle,ipoints,0,ipoints);
  
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int infile = 0;
  int counter = 0;
  bool condition = true;
  if(fType=="Lead"){
	  for(int i=0; i<nentries; i++){
		  tree->GetEntry(i);
		  condition = InterpretCut(sig,cut);
		  if(condition){
			  counter++;
			  if(counter!=number) continue;
			  infile = sizeof(x)*ipoints*i;
			  if(bl){
				  input.seekg(infile);
				  baseline = 0;
				  for(int ii=0; ii<gBaselineMax; ii++){
					  input.read((char*)&x,sizeof(x));
					  baseline += x/gmV;
				  }
				  baseline = baseline/gBaselineMax;
			  }
			  input.seekg(infile);
			  for(int ii=1; ii<ipoints+1; ii++){
				  input.read((char*)&x,sizeof(x));
				  if(bl) fSignal->SetBinContent(ii,(x/gmV)-baseline);
				  else   fSignal->SetBinContent(ii,(x/gmV));
			  }
		  }
	  }
  	  input.close();
  }
  else{
	  for(int i=0; i<nentries; i++){
		  tree->GetEntry(i);
		  iTree->GetEntry(i);
		  condition = InterpretCut(sig,cut);
		  if(condition){
			  counter++;
			  if(counter!=number) continue;
			  for(int ii=0; ii<ipoints; ii++){
				  fSignal->SetBinContent(ii+1,(*iVolt)[ii]);
			  }
		  }
	  }


  }
  
  return fSignal;
}
//------------------------------------------------------------------
///Resets all members of the class to their default values
void SFData::Reset(void){
 fSeriesNo      = 0;
 fNpoints       = 0;
 fFiber         = "dummy";
 fSource        = "dummy";
 fDesc          = "dummy"; 
 fType		= "dummy"; 
 fThreshold     = "dummy"; 
 fSpectrum      = NULL;
 fHist          = NULL;
 fHist2D        = NULL;
 fSignalProfile = NULL;
 fSignal        = NULL;
 fNames.clear();
 fPositions.clear();
 fTimes.clear();
 fSpectra.clear();
 fHists.clear();
 fHists2D.clear();
}
//------------------------------------------------------------------
///Prints details of currently analyzed experimental series
void SFData::Print(void){
 cout << "\n\n------------------------------------------------" << endl;
 cout << "This is Print() for SFData class object" << endl;
 cout << "Number of the experimental series: " << fSeriesNo << endl;
 cout << fDesc << endl;
 cout << "Type of measurements: " << fType << endl;
 cout << "Type of threshold used in DD6 data analysis: " << fThreshold << endl;
 cout << "Number of measurements in this series: " << fNpoints << endl;
 cout << "Fiber: " << fFiber << endl;
 cout << "Radioactive source: " << fSource << endl;
 cout << "List of measurements in this series:" << endl;
 for(int i=0; i<fNpoints; i++){
  cout << setw(30);
  cout << fNames[i] << "\t\t" << Form("%.1f mm",fPositions[i]) 
       << "\t\t" << Form("%.1f s",fTimes[i]) << endl; 
 }
 cout << "\n" << endl;
}
//------------------------------------------------------------------
