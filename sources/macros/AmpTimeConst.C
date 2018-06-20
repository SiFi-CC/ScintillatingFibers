R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFData.hh"

//Finds fit range for TH1D
bool GetFitRange(TH1D *sig, double &xmin, double &xmax){
 
  int nbins   = sig->GetNbinsX();
  double cont = 0;
  int istart  = 0; 
  
  for(int i=0; i<nbins; i++){
    cont = sig->GetBinContent(i);
    if(cont>354){
      istart = i;
      break;
    }
  }
  
  double xstart = sig->GetBinCenter(istart);
  sig->GetXaxis()->SetRange(xstart,xstart+300);
  int imin = sig->GetMaximumBin();
  xmin = sig->GetBinCenter(imin);
  xmax = xmin + 200;
  sig->GetXaxis()->UnZoom();
  
  return true;
}

//Finds fit range for TProfile
bool GetFitRange(TProfile *sig, double &xmin, double &xmax){
 
  int nbins   = sig->GetNbinsX();
  double cont = 0;
  int istart  = 0; 
  
  for(int i=0; i<nbins; i++){
    cont = sig->GetBinContent(i);
    if(cont>354){
      istart = i;
      break;
    }
  }
  
  double xstart = sig->GetBinCenter(istart);
  sig->GetXaxis()->SetRange(xstart,xstart+300);
  int imin = sig->GetMaximumBin();
  xmin = sig->GetBinCenter(imin);
  xmax = xmin + 200;
  sig->GetXaxis()->UnZoom();
  
  return true;
}

//Calculates amplifier time constant as average of 50 signals
bool AmpTimeConst(int ch=0){
  
  //----- Accessing data & fitting
  
  SFData *data = new SFData(9);
  data->Print();
  
  int nsig    = 50;	//number of signals for time const calculation
  int status  = -1;	//fit status
  int counter = 0;	//counter of signals with succesful fit
  int i       = 0;	//iterator for accessing signals
  double xmin, xmax;	//fit ranges
  
  int position = 0;		//source position
  if(ch==0) position = 1;
  else if(ch==1) position = 2;
  
  vector <TH1D*>  signals;
  vector <double> dectime;
  TH1D *temp;
  
  TString cut = Form("ch_%i.fPE>0.5 && ch_%i.fPE<1.5",ch,ch);
  TF1 *fun = new TF1("fun","[0]*TMath::Exp(-x/[1])+[2]",0,1024);
  
  while(counter<nsig){
    i++;
    temp = data->GetSignal(ch,position,cut,i,false);
    cout << "----- Fitting " << temp->GetName() << endl;
    fun->SetParameters(100,10,350);
    GetFitRange(temp,xmin,xmax);
    status = temp->Fit(fun,"Q","",xmin,xmax);
    cout << "\t Fit status: " << status << endl << endl;
    if(status==0 && 
      fun->GetParameter(1)>1E-5 && 
      fun->GetParameter(1)<80   &&
      xmin<700){
      signals.push_back(temp);
      dectime.push_back(fun->GetParameter(1));
      counter++;
    }
  }
  
  //----- Drawing sample of 12 signals
  //----- & calculating average decay time
  
  TString cname = Form("signals_ch%i",ch);
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
  cout << "stdev: " << TMath::StdDev(&dectime[0],&dectime[counter-1])/sqrt(counter)<< " ns" << endl; 
  
  return true;
}

//Calculates amplifier time constant from TProfile of 30 signals
bool TimeConstAv(int ch=0){
  
  SFData *data = new SFData(9);
  data->Print();
  
  TString cut  = Form("ch_%i.fPE>0.5 && ch_%i.fPE<1.5",ch,ch);
  int position = 0;
  if(ch==0) position = 1;
  else if(ch==1) position =2;
  double xmin, xmax;
  
  TProfile *signal = data->GetSignalAverage(ch,position,cut,30,false);
  signal->GetXaxis()->SetTitle("time [ns]");
  signal->GetYaxis()->SetTitle("amplitude [mV]");
  signal->SetStats(0);
  GetFitRange(signal,xmin,xmax);
  TF1 *fun = new TF1("fun","[0]*TMath::Exp(-x/[1])+[2]",xmin,xmax);
  fun->SetParameters(100,10,350);
  signal->Fit(fun,"R");
  gPad->SetGrid(1,1);
  
  cout << "\n\n" << endl;
  cout << "Decay time: " << Form("%.4f ns",fun->GetParameter(1)) << endl;
  
  return true;
}