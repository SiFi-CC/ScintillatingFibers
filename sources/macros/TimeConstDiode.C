R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFData.hh"
#include <fstream>
#include <istream>
using namespace std;

bool GetSignals(TString path="/home/kasia/data/2018_01_16_13_04/", int ch=1){
 
  TFile *file = new TFile(Form("signals_ch%i.root",ch),"RECREATE");
  TString fname = path+Form("wave_%i.dat",ch);
  ifstream input(fname, ios::binary);
  
  int nsig    = 200; 	//number of signals to be written in the file
  int ipoints = 1024;	//number of points in one signal
  int counter = 1;
  float x;
  
  TH1D *hsig; 
  TString hname;
  
  for(int i=0; i<nsig; i++){
   if(i%10==0) cout << i << "/" << nsig << "signals loaded..." << endl;
   hname = Form("sig_%i",i);
   hsig = new TH1D(hname,hname,ipoints,0,ipoints);
   for(int ii=0; ii<ipoints; ii++){
     input.read((char*)&x,sizeof(x));
     hsig->SetBinContent(ii,x/4.096);
   }
   hsig->Write();
  }
  
  input.close();
  file->Close();
  
  return true;
}
 
 bool GetFitRange(TH1D *sig, double &xmin, double &xmax){
 
  int nbins   = sig->GetNbinsX();
  double cont = 0;
  int istart  = 0; 
  
  for(int i=0; i<nbins; i++){
    cont = sig->GetBinContent(i);
    if(cont>50){
      istart = i;
      break;
    }
  }
  
  double xstart = sig->GetBinCenter(istart);
  sig->GetXaxis()->SetRange(xstart,xstart+400);
  int imin = sig->GetMaximumBin();
  xmin = sig->GetBinCenter(imin)+5;
  xmax = xmin + 300;
  sig->GetXaxis()->UnZoom();
  
  return true;
}

bool TimeConstDiode(int ch=1){

  GetSignals();
  
  TFile *file = new TFile(Form("signals_ch%i.root",ch),"READ");
  if(!file->IsOpen()){
   cout << "Could not open requested file!" << endl;
   return false;
  }
  
  int nsig    = 50;	//number of signals for time constant calculation
  int counter = 0;	//counter of signals with succesful fit
  int status  = -1;	//fit status;
  int i       = -1;	//iterator for accessing signals
  double xmin, xmax;	//fit ranges
  TString hname;	//histogram name
  
  TH1D *temp;
  vector <TH1D*>  signals;
  vector <double> dectime;
  
  TF1 *fun = new TF1("fun","[0]*TMath::Exp(-x/[1])+[2]",0,1024);
  
  while(counter<nsig){
   i++;
   hname = Form("sig_%i",i);
   temp = (TH1D*)file->Get(hname);
   cout << "----- Fitting " << hname << endl;
   fun->SetParameters(200,50,10);
   GetFitRange(temp,xmin,xmax);
   status = temp->Fit(fun,"Q","",xmin,xmax);
   cout << "\t Fit status: " << status << endl << endl;
   if(status==0 && 
      fun->GetParameter(1)>1E-5 && 
      fun->GetParameter(1)<100  &&
      xmin<700){
      signals.push_back(temp);
      dectime.push_back(fun->GetParameter(1));
      counter++;
    }
  }
  
  //----- Drawing sample of 12 signals
  //----- & calculating average decay time
  
  TString cname = Form("signals_diode_ch%i",ch);
  TCanvas *can = new TCanvas(cname,cname,1200,900);
  can->Divide(4,3);
  
  double sum = 0;
   
  for(int i=0; i<counter; i++){
    if(i<13){
     can->cd(i+1);
     gPad->SetGrid(1,1);
     signals[i]->SetStats(0);
     signals[i]->GetXaxis()->SetTitle("time [ns]");
     signals[i]->GetYaxis()->SetTitle("amplitude [mV]");
     signals[i]->Draw();
    }
    cout << Form("%.3f ns",dectime[i]) << endl;
    sum += dectime[i];
  }
  
  double average = sum/counter;
  cout << "\n\n" << endl;
  cout << "sum: " << sum << endl;
  cout << "counter:" << counter << endl;
  cout << "average decay time: " << average << " ns" << endl;
  cout << "stdev: " << TMath::StdDev(&dectime[0],&dectime[counter-1])<< " ns" << endl; 
  
  return true;
}