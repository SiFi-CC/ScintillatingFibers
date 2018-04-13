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
double gUnique = 0.;
char   *gPath = getenv("SFDATA");
int    gBaselineMax = 50;
double gmV = 4.096;
//------------------------------------------------------------------
TString SFData::fNames_1[9] = {"2018_01_20_13_56","2018_01_20_14_43","2018_01_20_15_30",
			       "2018_01_22_9_29","2018_01_20_13_07","2018_01_22_10_16",
			       "2018_01_22_11_02","2018_01_22_11_51","2018_01_22_12_40"};
TString SFData::fNames_2[9] = {"2018_01_22_16_59","2018_01_22_17_54","2018_01_22_18_42",
			       "2018_01_23_9_54","2018_01_23_10_46","2018_01_23_11_38",
			       "2018_01_23_12_27","2018_01_23_14_56","2018_01_23_15_46"};
TString SFData::fNames_3[9] = {"2018_01_25_9_20","2018_01_25_9_52","2018_01_25_10_23",
			       "2018_01_25_10_54","2018_01_25_8_47","2018_01_25_11_26",
			       "2018_01_25_12_00","2018_01_25_12_32","2018_01_25_13_10"};
TString SFData::fNames_4[9] = {"2018_01_25_14_05","2018_01_25_14_36","2018_01_25_15_08",
			       "2018_01_25_15_39","2018_01_25_16_11","2018_01_25_16_47",
			       "2018_01_25_17_19","2018_01_25_17_50","2018_01_26_8_27"};
TString SFData::fNames_5[9] = {"2018_01_26_10_03","2018_01_26_10_36","2018_01_26_11_07",
			       "2018_01_26_11_39","2018_01_26_12_11","2018_01_26_12_42",
			       "2018_01_26_13_14","2018_01_26_13_46","2018_01_26_14_18"};
TString SFData::fNames_6[5] = {"2018_01_20_13_07","2018_01_24_16_27","2018_01_24_17_25",
			       "2018_01_24_18_27","2018_01_25_8_47"};
TString SFData::fNames_7[1] = {"2018_01_22_15_13"};
TString SFData::fNames_8[1] = {"2018_01_23_14_07"};
TString SFData::fNames_9[2] = {"2018_01_19_17_17","2018_01_22_15_37"};
//------------------------------------------------------------------
double SFData::fPositions_1[9] = {10.3,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
double SFData::fPositions_2[9] = {10.2,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
double SFData::fPositions_3[9] = {10.4,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
double SFData::fPositions_4[9] = {10.5,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
double SFData::fPositions_5[9] = {10.4,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
double SFData::fPositions_6[5] = {50.0,50.0,50.0,50.0,50.0};
double SFData::fPositions_7[1] = {0.00};
double SFData::fPositions_8[1] = {0.00};
double SFData::fPositions_9[2] = {0.00,0.00};

//------------------------------------------------------------------
///Default constructor. If this constructor is used the series 
///number should be set via SetDetails(int seriesNo) function.
SFData::SFData(){
 Reset();
 cout << "##### Warning in SFData constructor!" << endl;
 cout << "You are using the default constructor. Set the series number." <<endl;
}
//------------------------------------------------------------------
///Standard constructor (recommended).
///\param seriesNo is number of experimental series to analyze.
SFData::SFData(int seriesNo){
  SetDetails(seriesNo);
}
//------------------------------------------------------------------
SFData::SFData(int seriesNo, TString threshold){
 SetDetails(seriesNo);
 if(!(threshold=="ft" || threshold=="cf")){
  throw "##### Error in SFData constructor! Incorrect threshold type!"; 
 }
 fThreshold = threshold; 
}
//------------------------------------------------------------------
///Default destructor.
SFData::~SFData(){
}
//------------------------------------------------------------------
///Sets all details of selected experimental series. If default constructor
///was used, this function needs to be called explicitly with the number of 
///of requested series as an argument. The following attributes are set 
///within this function:
bool SFData::SetDetails(int seriesNo){
 
  if(seriesNo>9 || seriesNo<1){
    cout << "##### Error in SFData::SetDetails()!" << endl;
    cout << "There is no " << seriesNo << " series in the database!" << endl;
    return false;
  }
  
  Reset();
  fSeriesNo = seriesNo;
  fThreshold = "ft";	//default - fixed threshold
  
  ///- number of measurements in the series
  if(fSeriesNo<6) fNpoints = 9;
  else if(fSeriesNo==6) fNpoints = 5;
  else if(fSeriesNo==7 || fSeriesNo==8) fNpoints = 1;
  else if(fSeriesNo==9) fNpoints = 2;
  
  ///- fiber type
  if(fSeriesNo==2 || fSeriesNo==4 || fSeriesNo==8) fFiber = "LuAG:Ce (2)";
  else if(fSeriesNo==9) fFiber = "none";
  else fFiber = "LuAG:Ce (1)";
  
  ///- measurement description
  switch(fSeriesNo){
    case 1: fDesc = "First measurement";   break;
    case 2: fDesc = "First measurement";   break;
    case 3: fDesc = "Reameasuring";	   break;
    case 4: fDesc = "Reameasuring";	   break;
    case 5: fDesc = "With Al coating";	   break;
    case 6: fDesc = "Different couplings"; break;
    case 7: fDesc = "Internal activity";   break;
    case 8: fDesc = "Internal activity";   break;
    case 9: fDesc = "PE callibration";	   break;
  }
  
  ///- positions of radioactive source for all measurements
  switch(fSeriesNo){
    case 1: fPositions = fPositions_1; break;
    case 2: fPositions = fPositions_2; break;
    case 3: fPositions = fPositions_3; break;
    case 4: fPositions = fPositions_4; break;
    case 5: fPositions = fPositions_5; break;
    case 6: fPositions = fPositions_6; break;
    case 7: fPositions = fPositions_7; break;
    case 8: fPositions = fPositions_8; break;
    case 9: fPositions = fPositions_9; break;
  }
  
  ///- names of measurements in this series 
  switch(fSeriesNo){
    case 1: fNames = fNames_1; break;
    case 2: fNames = fNames_2; break;
    case 3: fNames = fNames_3; break;
    case 4: fNames = fNames_4; break;
    case 5: fNames = fNames_5; break;
    case 6: fNames = fNames_6; break;
    case 7: fNames = fNames_7; break;
    case 8: fNames = fNames_8; break;
    case 9: fNames = fNames_9; break;
  }
   
  return true;
}
//------------------------------------------------------------------
/// Sets the flag to identify the tree with data which should be accessed.
/// \param threshold - type of threshold used during data analysis in DD6.
/// Possible options are: ft - fixed threshold (default) and cf - constant
/// fraction. 
bool SFData::SetThreshold(TString threshold){
  if(!(threshold=="ft" || threshold=="cf")){
    cout << "##### Error in SFData::SetThreshold()! Incorrect threshold type!" << endl;
    cout << "Possible options are: ft for fixed threshold and cf for constant fraction" << endl;
    return false;
  }
  fThreshold = threshold;
  return true;
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
/// Returns index in the fNames and fPositions arrays for the 
/// measurement of requested source position in mm. 
/// If measurements in analyzed series don't have unique positions 
/// a number of measurement should be passed. Measurements counting 
/// starts at 1.
int SFData::GetIndex(double position){
   
  int index = -1;
  
  if(fSeriesNo>5){
    index = position-1;
    return index;
  }

  for(int i=0; i<fNpoints; i++){
    if(fabs(fPositions[i]-position)<1){
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
  TString hname = type+Form("_ch%i_pos%.1f_",ch,position)+fThreshold;
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
 
  TString fname;
  TString tname = string("tree_")+fThreshold;
  TFile *file;
  TTree *tree;
  TString selection;
  TString hname, htitle;
  bool empty = fSpectra.empty();
  
  for(int i=0; i<fNpoints; i++){
   fname = string(gPath)+fNames[i]+"/results.root";
   file = new TFile(fname,"READ");
   tree = (TTree*)file->Get(tname);
   selection = GetSelection(ch,type);
   tree->Draw(selection,cut);
   if(empty) fSpectra.push_back(new TH1D());
   fSpectra[i] = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
   hname = type+Form("_ch%i_pos%.1f_",ch,fPositions[i])+fThreshold;
   htitle = hname+" "+cut;
   fSpectra[i]->SetName(hname);
   fSpectra[i]->SetTitle(htitle);
  }
  
  return fSpectra;
}
//------------------------------------------------------------------
///Returns a vector of requested ratio histograms for all measurements in this series.
///\param selection - selection like for Draw() method of TTree, e.g. "log(ch_0.fPE/ch_1.fPE)". 
///Redirection to specific histogram should not be entered here, it is done inside the function.
///\param cut - cut for drawn events. Also TTree-style syntax. If empty string is passed here 
///all events will be drawn.
vector <TH1D*> SFData::GetRatios(TString selection, TString cut){
  
  TString fname;
  TString tname = string("tree_")+fThreshold;
  TFile *file;
  TTree *tree;
  TString selectAndDraw;
  TString hname, htitle;
  bool empty = fRatios.empty();
  
  for(int i=0; i<fNpoints; i++){
   fname = string(gPath)+fNames[i]+"/results.root";
   file = new TFile(fname,"READ");
   tree = (TTree*)file->Get(tname);
   gUnique = gRandom->Uniform(0,1);
   selectAndDraw = selection+Form(">>htemp%.7f(500,-5,5)",gUnique);
   tree->Draw(selectAndDraw,cut);
   if(empty) fRatios.push_back(new TH1D());
   fRatios[i] = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
   hname = Form("ratio_pos%.1f_",fPositions[i])+fThreshold;
   htitle = selection+" "+cut;
   fRatios[i]->SetName(hname);
   fRatios[i]->SetTitle(htitle);
  }
  
  return fRatios;
}
//------------------------------------------------------------------
///Returns single requested ratio histogram.
///\param selection - selection like for Draw() method of TTree.
///\param cut - cut for drawn events. Also TTree-style syntax.
///\param position - position of source in mm. If position is not unique a 
///number of measurement should be entered.
TH1D* SFData::GetRatio(TString selection, TString cut, double position){
  
  int index = GetIndex(position);
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TString tname = string("tree_")+fThreshold;
  TTree *tree = (TTree*)file->Get(tname);
  gUnique = gRandom->Uniform(0,1);
  TString selectAndDraw = selection+Form(">>htemp%.7f",gUnique);
  tree->Draw(selectAndDraw,cut);
  fRatio = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
  TString hname = Form("ratio_pos%.1f_",position)+fThreshold;
  TString htitle = selection+" "+cut;
  fRatio->SetName(hname);
  fRatio->SetTitle(htitle);
  
  return fRatio;
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
  
  TString iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
  ifstream input(iname,ios::binary);
 
  TString hname = Form("sig_ch%i_pos_%.1f_num_%i_",ch,position,number)+fThreshold;
  TString htitle = hname +" "+cut;
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
     if(counter<number) counter++;
     else break;
    }
  }
  
  if(counter<number) 
    cout << "##### Warning in SFData::GetSignalAverage()! " << counter 
         << " out of " << number << " plotted." << endl;
  
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
  
  TString iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
  ifstream input(iname,ios::binary);
  
  TString hname = Form("sig_ch%i_pos_%.1f_no%i_",ch,position,number)+fThreshold;
  TString htitle = hname+" "+cut;
  
  fSignal = new TH1D(hname,htitle,ipoints,0,ipoints);
  
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int infile = 0;
  int counter = 0;
  bool condition = true;
  
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
  
  return fSignal;
}
//------------------------------------------------------------------
///Resets all members of the class to their default values
void SFData::Reset(void){
 fSeriesNo      = 0;
 fNpoints       = 0;
 fFiber         = "dummy";
 fDesc          = "dummy"; 
 fThreshold     = "dummy"; 
 fNames         = NULL;
 fPositions     = NULL;
 fSpectrum      = NULL;
 fSignalProfile = NULL;
 fSignal        = NULL;
}
//------------------------------------------------------------------
///Prints details of currently analyzed experimental series
void SFData::Print(void){
 cout << "\n\n------------------------------------------------" << endl;
 cout << "This is Print() for SFData class object" << endl;
 cout << "Number of the experimental series: " << fSeriesNo << endl;
 cout << fDesc << endl;
 cout << "Type of threshold used in DD6 data analysis: " << fThreshold << endl;
 cout << "Number of measurements in this series: " << fNpoints << endl;
 cout << "Fiber: " << fFiber << endl;
 cout << "List of measurements in this series:" << endl;
 for(int i=0; i<fNpoints; i++){
  cout << setw(20);
  cout << fNames[i] << "\t" << Form("%.1f",fPositions[i]) << " mm" << endl; 
 }
 cout << "\n" << endl;
}
//------------------------------------------------------------------
