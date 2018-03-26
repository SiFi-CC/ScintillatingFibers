// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFAttenuation.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFAttenuation.hh"
#include "TFile.h"

ClassImp(SFAttenuation);

//------------------------------------------------------------------
///Default constructor.
SFAttenuation::SFAttenuation(){
  cout << "#### Warning in SFAttenuation constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  fSeriesNo = -1;
  fData = NULL;
}
//------------------------------------------------------------------
///Standard constructor (reccommended)
///\param seriesNo is number of experimental series to be analyzed. 
SFAttenuation::SFAttenuation(int seriesNo){
  fSeriesNo = seriesNo;
  fData = new SFData(fSeriesNo);
}
//------------------------------------------------------------------
///Default destructor.
SFAttenuation::~SFAttenuation(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------
///Method to determine attenuation length used in Pauwels et al., JINST 8 (2013) P09019.
///For both ends of the fiber one value is calculated, since averaged signal from both channels
///is taken into account.
bool SFAttenuation::AttAveragedCh(void){
 
  int npoints = fData->GetNpoints();
  double *positions = fData->GetPositions();
  TString selection = "log(sqrt(ch_1.fPE/ch_0.fPE))";
  TString cut = "ch_0.fT0>0 && ch_0.fT0<590 && ch_1.fT0>0 && ch_1.fT0<590";
  vector <TH1D*> ratio = fData->GetRatios(selection,cut);
  
  double mean, sigma;
  vector <TF1*> fun;
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("log #(){#sqrt{#frac{ch0}{ch1}}}");
  graph->SetMarkerStyle(4);
  
  for(int i=0; i<npoints; i++){
    mean = ratio[i]->GetMean();
    sigma = ratio[i]->GetRMS();
    fun.push_back(new TF1("fun","gaus",mean-sigma,mean+sigma));
    ratio[i]->Fit(fun[i],"QR");
    graph->SetPoint(i,positions[i],fun[i]->GetParameter(1));
    graph->SetPointError(i,0,fun[i]->GetParError(1));
  }
  
  TF1 *fpol1 = new TF1("fpol1","pol1",positions[0],positions[npoints-1]);
  graph->Fit(fpol1,"QR");
  double attenuation = fabs(1./fpol1->GetParameter(1));
  double att_error = fpol1->GetParError(1)/pow(fpol1->GetParameter(1),2);
  
  cout << "Attenuation lenght is: " << attenuation << " +/- " << att_error << " mm" << endl;

  TFile *file = new TFile("attenuation_no_coating_fib2.root","UPDATE");
  graph->Write();
  for(int ii=0; ii<npoints; ii++){
   ratio[ii]->Write(); 
  }
  file->Close();
  
  return true;
}
//------------------------------------------------------------------
///Method to determine attenuation length for both channels independently. 
///\param ch - channel number
bool SFAttenuation::AttSeparateCh(int ch){
 
  int npoints = fData->GetNpoints();
  double *positions = fData->GetPositions();
  TString cut = Form("ch_%i.fT0>0 && ch_%i.fT0<590 && ch_%i.fPE>0",ch,ch,ch);
  vector <TH1D*> spectra = fData->GetSpectra(ch,"fPE",cut);
  
  TGraphErrors *graph = new TGraphErrors(npoints);
  graph->GetXaxis()->SetTitle("source position [mm]");
  graph->GetYaxis()->SetTitle("P.E.");
  graph->SetMarkerStyle(4);
  
  vector <SFPeakFinder*> peakfin;
  vector <TH1D*> peaks;
  vector <TF1*> fun;
  double mean, sigma;
  
  for(int i=0; i<npoints; i++){
    peakfin.push_back(new SFPeakFinder(spectra[i],"511"));
    peaks.push_back(peakfin[i]->GetPeak());
    mean = peaks[i]->GetMean();
    sigma = peaks[i]->GetRMS();
    fun.push_back(new TF1("fun","gaus",mean-3*sigma,mean+3*sigma));
    peaks[i]->Fit(fun[i],"QR");
    graph->SetPoint(i,positions[i],fun[i]->GetParameter(1));
    graph->SetPointError(i,0,fun[i]->GetParError(1));
  }
  
  TF1 *fexp = new TF1("fexp","expo",positions[0],positions[npoints-1]);
  graph->Fit(fexp,"QR");
  double attenuation = fabs(1./fexp->GetParameter(1));
  double att_error = fexp->GetParError(1)/pow(fexp->GetParameter(1),2);
  
  cout << "Attenuation for channel "<< ch << ": " << attenuation << " +/- " << att_error << " mm" << endl;
  
  return true;
}
//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFAttenuation::Print(void){
 cout << "This is print out of SFAttenuation class object" << endl;
 cout << "Experimental series number " << fSeriesNo << endl;
}
//------------------------------------------------------------------