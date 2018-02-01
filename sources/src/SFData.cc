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
using namespace std;

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
//------------------------------------------------------------------
SFData::SFData(){
 Reset();
 cout << "##### Warning in SFData constructor!" << endl;
 cout << "You are using the default constructor. Set the series number." <<endl;
}
//------------------------------------------------------------------
SFData::SFData(int seriesNo){
  SetDetails(seriesNo);
}
//------------------------------------------------------------------
SFData::~SFData(){
}
//------------------------------------------------------------------
bool SFData::SetDetails(int seriesNo){
 
  if(seriesNo>9 || seriesNo<1){
    cout << "##### Error in SFData::SetDetails()!" << endl;
    cout << "There is no " << seriesNo << " series in the database!" << endl;
    return false;
  }
  
  Reset();
  fSeriesNo = seriesNo;
  
  //number of measurements in the series
  if(fSeriesNo<6) fNpoints = 9;
  else if(fSeriesNo==6) fNpoints = 5;
  else if(fSeriesNo==7 || fSeriesNo==8) fNpoints = 1;
  else if(fSeriesNo==9) fNpoints = 2;
  
  //fiber type
  if(fSeriesNo==2 || fSeriesNo==4 || fSeriesNo==8) fFiber = "LuAG:Ce (2)";
  else if(fSeriesNo==9) fFiber = "none";
  else fFiber = "LuAG:Ce (1)";
  
  //measurement description
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
  
  //positions for all measurements
  switch(fSeriesNo){
    case 1: fPositions = fPositions_1; break;
    case 2: fPositions = fPositions_2; break;
    case 3: fPositions = fPositions_3; break;
    case 4: fPositions = fPositions_4; break;
    case 5: fPositions = fPositions_5; break;
    case 6: fPositions = fPositions_6; break;
  }
  
  //names of measurements in this series 
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
TString SFData::GetSelection(int ch, TString type){
 
  TString selection;
  gUnique = gRandom->Uniform(0,1);
  
  if(type=="fAmp")
    selection = Form("ch_%i.fAmp>>htemp%.7f(1000,0,700)",ch,gUnique);
  else if(type=="fCharge")
    selection = Form("ch_%i.fCharge>>htemp%.7f(1000,-1E4,2E5)",ch,gUnique);
  else if(type=="fPE")
    selection = Form("ch_%i.fPE>>htemp%.7f(1000,-25,1000)",ch,gUnique);
  else if(type=="fT0")
    selection = Form("ch_%i.fT0>>htemp%.7f(1000,-110,500)",ch,gUnique);
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
int SFData::GetIndex(double position){
 
  int index = -1;
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
TH1D* SFData::GetSpectrum(int ch, TString type, TString cut, double position){

  int index = GetIndex(position);
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  
  TString selection = GetSelection(ch,type);
  tree->Draw(selection,cut);
  fSpectrum = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
  TString hname = type+Form("_ch%i_pos%.1f",ch,position);
  fSpectrum->SetName(hname);
  fSpectrum->SetTitle(hname);
  
  return fSpectrum; 
}
//------------------------------------------------------------------
TH1D** SFData::GetSpectra(int ch, TString type, TString cut){
 
  TString fname;
  TFile *file;
  TTree *tree;
  TString selection;
  TString hname;
  
  for(int i=0; i<fNpoints; i++){
   fname = string(gPath)+fNames[i]+"/results.root";
   file = new TFile(fname,"READ");
   TTree *tree = (TTree*)file->Get("tree_ft");
   selection = GetSelection(ch,type);
   tree->Draw(selection,cut);
   fSpectra[i] = (TH1D*)gROOT->FindObjectAny(Form("htemp%.7f",gUnique));
   hname = type+Form("_ch%i_pos%.1f",ch,fPositions[i]);
   fSpectra[i]->SetName(hname);
   fSpectra[i]->SetTitle(hname);
  }
  
  return fSpectra;
}
//------------------------------------------------------------------
TProfile* SFData::GetSignalAverage(int ch, double position, int number, bool bl){
  
  int index = GetIndex(position);
  const int ipoints = 1024;
  float x;
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i",ch),&sig);
  
  TString iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
  ifstream input(iname,ios::binary);
 
  TString hname = Form("sigav_ch%i_pos_%.1f",ch,position);
  fSignalProfile = new TProfile(hname,hname,number*ipoints,0,ipoints,"");
    
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int counter = 0;
  int infile = 0;
  bool condition = true;
  
  for(int i=0; i<nentries; i++){
   tree->GetEntry(i);
   if(condition){
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
     if(counter<number-1) counter++;
     else break;
    }
  }
  
  input.close();
  
  return fSignalProfile;
}
//------------------------------------------------------------------
TH1D* SFData::GetSignal(int ch, double position, int number, bool bl){
 
  int index = GetIndex(position);
  const int ipoints = 1024;
  float x;
  
  TString fname = string(gPath)+fNames[index]+"/results.root";
  TFile *file = new TFile(fname,"READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  DDSignal *sig = new DDSignal();
  tree->SetBranchAddress(Form("ch_%i",ch),&sig);
  
  TString iname = string(gPath)+fNames[index]+Form("/wave_%i.dat",ch);
  ifstream input(iname,ios::binary);
  
  TString hname = Form("sig_ch%i_pos_%.1f_no%i",ch,position,number);
  
  fSignal = new TH1D(hname,hname,ipoints,0,ipoints);
  
  int nentries = tree->GetEntries();
  double baseline = 0.;
  int infile = 0;
  int counter = 0;
  bool condition = true;
  
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
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
void SFData::Reset(void){
 fSeriesNo  = 0;
 fNpoints   = 0;
 fFiber     = "dummy";
 fDesc      = "dummy"; 
 fNames     = NULL;
 fPositions = NULL;
 for(int i=0; i<9; i++){
   fSpectra[i] = NULL;
 }
 fSpectrum = NULL;
 fSignalProfile = NULL;
 fSignal = NULL;
}
//------------------------------------------------------------------
void SFData::Print(void){
 cout << "\n\n------------------------------------------------" << endl;
 cout << "This is Print() for SFData class object" << endl;
 cout << "Number of the experimental series: " << fSeriesNo << endl;
 cout << fDesc << endl;
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