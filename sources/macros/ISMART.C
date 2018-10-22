R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFAttenuation.hh"
#include "../include/SFTimingRes.hh"
#include "../include/SFData.hh"
#include "../include/SFTimeConst.hh"
#include "../include/SFFitResults.hh"
#include "../include/SFEnergyResolution.hh"
#include "/home/kasia/DesktopDigitizer6/sources/include/DDSignal.hh"

bool ISMART(void){
 
  const int series = 14;
  
  SFData *data = new SFData(series);
  data->Print();
  int n = data->GetNpoints();
  vector <double> pos = data->GetPositions();
  
  //----- attenuation - ratio histograms
  
  SFAttenuation *att;
  
  try{
    att = new SFAttenuation(series);
  }
  catch(const char *message){
   cout << message << endl;
   return false;
  }
  
  att->AttAveragedCh();
  vector <TH1D*> attRatios = att->GetRatios();
  
  TCanvas *canRat = new TCanvas("canRat","canRat",1000,1000);
  gPad->SetGrid(1,1);
  
  TLegend *legR = new TLegend(0.65,0.66,0.94,0.92);
  
  Color_t col[5] = {kBlack,kRed,kBlue,kMagenta,kGreen-2};
  int counter = 0;
  
  for(int i=0; i<n; i++){
    if(i==0 || i%2==0){
      attRatios[i]->SetLineColor(col[counter]);
      attRatios[i]->SetLineWidth(1);
      attRatios[i]->GetFunction("fun")->SetLineColor(col[counter]);
      attRatios[i]->SetStats(0);
      if(i==0){ 
        attRatios[i]->Draw();
        attRatios[i]->GetXaxis()->SetRangeUser(-0.8,0.6);
        attRatios[i]->GetYaxis()->SetRangeUser(0,35000);
        attRatios[i]->SetTitle("");
        attRatios[i]->GetXaxis()->SetTitle("ln(#sqrt{Q1/Q0})");
        attRatios[i]->GetYaxis()->SetMaxDigits(1);
      }
      else 
        attRatios[i]->Draw("same");
      legR->AddEntry(attRatios[i],Form("position %.1f mm",pos[i]),"PL");
      counter++;
    }
  }
  
  legR->Draw();
  
  canRat->SaveAs("canRat.png");
  
  delete att;
  
  //----- timing resolution - T0Diff histograms
 
  SFTimingRes *tim;
  
  try{
    tim = new SFTimingRes(series,"ft","with cut");
  }
  catch(const char *message){
   cout << message << endl;
   return false;
  }
  
  vector <TH1D*> T0Diff = tim->GetT0Diff();
  
  TCanvas *canT0 = new TCanvas("canT0","canT0",1000,1000);
  gPad->SetGrid(1,1);
  
  TLegend *legT0 = new TLegend(0.65,0.66,0.94,0.92);
  
  counter=0;
  
  for(int i=0; i<n; i++){
    if(i==0 || i%2==0){
     T0Diff[i]->SetLineColor(col[counter]);
     T0Diff[i]->GetFunction("fun")->SetLineColor(col[counter]);
     T0Diff[i]->SetStats(0);
     if(i==0){
      T0Diff[i]->Draw(); 
      T0Diff[i]->GetXaxis()->SetRangeUser(-4,4);
      T0Diff[i]->GetYaxis()->SetRangeUser(0,650);
      T0Diff[i]->SetTitle("");
      T0Diff[i]->GetXaxis()->SetTitle("time difference [ns]");
     }
     else 
       T0Diff[i]->Draw("same");
     legT0->AddEntry(T0Diff[i],Form("position %.1f mm",pos[i]),"PL");
     counter++;
    }
  }
  
  legT0->Draw();
  
  canT0->SaveAs("canT0.png");
  
  delete tim;
 
  return true;
}

bool Signals(void){
 
  SFData *dataS3 = new SFData(3);
  SFData *dataS14 = new SFData(14);
  
  SFTimeConst *tconstS3;
  SFTimeConst *tconstS14;
  
  TString descS3[3] = {"base line","fast component","slow component"};
  TString descS14[2] = {"base line","decay"};
  
  try{
   tconstS3 = new SFTimeConst(3,200,false); 
   tconstS14 = new SFTimeConst(14,200,false);
  }
  catch(const char *message){
   cout << message << endl;
   return false;
  }
  
  tconstS3->FitSingleSignal(0,50);
  TProfile *sigS3 = tconstS3->GetSingleSignal(0,50);
  SFFitResults *resS3 = tconstS3->GetSingleResult(0,50);
  sigS3->SetStats(0);
  sigS3->GetXaxis()->SetTitle("time [ns]");
  sigS3->GetYaxis()->SetTitle("amplitude [mV]");
  sigS3->GetXaxis()->SetTitleSize(0.045);
  sigS3->GetYaxis()->SetTitleSize(0.045);
  sigS3->GetXaxis()->SetLabelSize(0.04);
  sigS3->GetYaxis()->SetLabelSize(0.04);
  sigS3->SetTitle("");
  
  tconstS14->FitSingleSignal(0,50);
  TProfile *sigS14 = tconstS14->GetSingleSignal(0,50);
  SFFitResults *resS14 = tconstS14->GetSingleResult(0,50);
  sigS14->SetStats(0);
  sigS14->GetXaxis()->SetTitle("time [ns]");
  sigS14->GetYaxis()->SetTitle("amplitude [mv]");
  sigS14->GetXaxis()->SetTitleSize(0.045);
  sigS14->GetYaxis()->SetTitleSize(0.045);
  sigS14->GetXaxis()->SetLabelSize(0.04);
  sigS14->GetYaxis()->SetLabelSize(0.04);
  sigS14->SetTitle("");
  
  TCanvas *canS3 = new TCanvas("canS3","canS3",700,400);
  gPad->SetGrid(1,1);
  sigS3->Draw("HIST");
  resS3->GetFunction()->Draw("same");
  vector <TF1*> compS3 = resS3->GetCompFunctions();
  TLegend *legS3 = new TLegend(0.56,0.64,0.87,0.87);
  
  for(int i=0; i<3; i++){
   compS3[i]->Draw("same"); 
   legS3->AddEntry(compS3[i],descS3[i],"L");
  }
  
  legS3->AddEntry(resS3->GetFunction(),"sum","L");
  legS3->Draw();
  
  TCanvas *canS14 = new TCanvas("canS14","canS14",700,400);
  gPad->SetGrid(1,1);
  sigS14->Draw("HIST");
  resS14->GetFunction()->Draw("same");
  vector <TF1*> compS14 = resS14->GetCompFunctions();
  TLegend *legS14 = new TLegend(0.65,0.75,0.87,0.87);
  
  for(int i=0; i<2; i++){
    compS14[i]->Draw("same"); 
    legS14->AddEntry(compS14[i],descS14[i],"L");
  }
  
  legS14->AddEntry(resS14->GetFunction(),"sum","L");
  //legS14->Draw();
  
  canS3->SaveAs("sig_LuAG.png");
  canS14->SaveAs("sig_LYSO.png");
  
  delete tconstS3, tconstS14;
  
  return true;
}

double funDecayDouble(double *x, double *par){
 double fast_dec = par[0]*TMath::Exp(-(x[0]-par[1])/par[2]); 
 double slow_dec = par[3]*TMath::Exp(-(x[0]-par[1])/par[4]); 
 double constant = par[5];
 return fast_dec+slow_dec+constant;
}

bool GAGG(void){

  const int ipoints = 1024;
  float x;
  int position = 0;
  double baseline = 0;
  double PE = 200;
  int nsignals = 50;
  double T0first = 0;
  DDSignal *sig = new DDSignal();
  
  TString path = "/home/kasia/data/2018_01_23_17_59/";
  ifstream input(path+string("wave_0.dat"),ios::binary);
  TFile *file = new TFile(path+string("results.root"),"READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  tree->SetBranchAddress("ch_0",&sig);
  int nentries = tree->GetEntries();
  
  TProfile *hSig = new TProfile("sig","sig",ipoints,0,ipoints,"");
  
  for(int i=0; i<nentries; i++){
    tree->GetEntry(i);
    baseline = 0;
    if(sig->GetPE()>(PE-0.5) && sig->GetPE()<(PE+0.5)){
      if(T0first==0)
	T0first=sig->GetT0();
      else if(abs(sig->GetT0()-T0first)<1){
	position = sizeof(x)*ipoints*i;
	input.seekg(position);
	for(int ii=0; ii<50; ii++){
	  input.read((char*)&x,sizeof(x));
	  baseline+=x/4.096;
	}
	baseline=baseline/50;
	input.seekg(position);
	for(int ii=0; ii<ipoints; ii++){
	 input.read((char*)&x,sizeof(x));
	 hSig->Fill(ii,(x/4.096)-baseline);
	}
      }
    }
  }
  
  TCanvas *can = new TCanvas("can","can",700,400);
  gPad->SetGrid(1,1);
  hSig->SetStats(0);
  hSig->GetXaxis()->SetTitle("time [ns]");
  hSig->GetYaxis()->SetTitle("amplitude [mV]");
  hSig->GetXaxis()->SetTitleSize(0.045);
  hSig->GetYaxis()->SetTitleSize(0.045);
  hSig->GetXaxis()->SetLabelSize(0.04);
  hSig->GetYaxis()->SetLabelSize(0.04);
  hSig->SetTitle("");
  hSig->Draw("HIST");
  
  double xmin = hSig->GetBinCenter(hSig->GetMaximumBin())+20.;
  double xmax = hSig->GetBinCenter(hSig->GetNbinsX());
  
  TF1 *fun_BL = new TF1("fun_BL","pol0",0,50); 
  hSig->Fit(fun_BL,"RQ");
  
  TF1* fun_fast = new TF1("fun_fast","[0]*exp(-(x-[1])/[2])",xmin,xmin+20);
  fun_fast->SetParameters(100.,100.,10.);
  fun_fast->SetParNames("A","t0","tau");
  hSig->Fit(fun_fast,"RQ");
  
  TF1* fun_slow = new TF1("fun_slow","[0]*exp(-(x-[1])/[2])",300,400);
  fun_slow->SetParameters(100.,100.,100.);
  fun_slow->SetParNames("A","t0","tau");
  fun_slow->FixParameter(1,fun_fast->GetParameter(1));
  hSig->Fit(fun_slow,"RQ");
  
  TF1* fun_all = new TF1("fun_all",funDecayDouble,xmin,xmax,6);
  fun_all->SetParNames("A_fast","t0","tau_fast","A_slow","tau_slow","const");
  fun_all->SetParameter(0,fun_fast->GetParameter(0));
  fun_all->FixParameter(1,fun_fast->GetParameter(1));
  fun_all->SetParameter(2,fun_fast->GetParameter(2));
  fun_all->SetParameter(3,fun_slow->GetParameter(0));
  fun_all->SetParameter(4,fun_slow->GetParameter(2));
  fun_all->FixParameter(5,fun_BL->GetParameter(0));
  int fitStat = hSig->Fit(fun_all,"RQ");
  
  SFFitResults *res = new SFFitResults();
  res->SetFromFunction(fun_all);
  res->Print();
  
  res->GetFunction()->Draw("same");
  vector <TF1*> comp = res->GetCompFunctions();
  TString desc[3] = {"base line","fast component","slow component"};
  
  TLegend *leg = new TLegend(0.59,0.66,0.87,0.87);
  
  for(int i=0; i<3; i++){
    comp[i]->Draw("same");
    leg->AddEntry(comp[i],desc[i],"L");
  }
  
  leg->AddEntry(res->GetFunction(),"sum","L");
  //leg->Draw();
  
  can->SaveAs("sig_GAGG.png");
  
  return true;
}

bool PESpectra(void){
 
  SFData *dataS3 = new SFData(3);
  SFData *dataS14 = new SFData(14);
  
  TH1D *histS3 = dataS3->GetSpectrum(0,"fPE","ch_0.fPE>0",50);
  TH1D *histS14 = dataS14->GetSpectrum(0,"fPE","ch_0.fPE>0",50);
  
  delete dataS3;
  delete dataS14;
  
  TString fname = "/home/kasia/data/2018_01_23_17_59/results.root";
  TFile *file = new TFile(fname,"READ");
  TTree *tree = (TTree*)file->Get("tree_ft");
  TH1D *histGAGG = new TH1D();
  tree->Draw("ch_0.fPE>>htemp(1000,-150,1200)","ch_0.fPE>0");
  histGAGG = (TH1D*)gROOT->FindObjectAny("htemp");
  
  double integral = 0;
  double max = 0;
  
  TCanvas *can = new TCanvas("canSpec","canSpec",1800,500);
  can->Divide(3,1);
  
  can->cd(1);
  gPad->SetGrid(1,1);
  histS3->SetStats(0);
  histS3->GetXaxis()->SetTitle("charge [P.E.]");
  histS3->GetYaxis()->SetTitle("counts normalized");
  histS3->GetXaxis()->SetTitleSize(0.045);
  histS3->GetYaxis()->SetTitleSize(0.045);
  histS3->GetXaxis()->SetLabelSize(0.045);
  histS3->GetYaxis()->SetLabelSize(0.045);
  histS3->GetYaxis()->SetMaxDigits(2);
  histS3->SetTitle("");
  integral = histS3->Integral();
  histS3->Scale(1./integral);
  max = histS3->GetBinContent(histS3->GetMaximumBin());
  histS3->GetXaxis()->SetRangeUser(0,800);
  histS3->GetYaxis()->SetRangeUser(0,max+0.1*max);
  histS3->Draw();
  
  can->cd(2);
  gPad->SetGrid(1,1);
  histS14->SetStats(0);
  histS14->GetXaxis()->SetTitle("charge [P.E.]");
  histS14->GetYaxis()->SetTitle("counts normalized");
  histS14->GetXaxis()->SetTitleSize(0.045);
  histS14->GetYaxis()->SetTitleSize(0.045);
  histS14->GetXaxis()->SetLabelSize(0.045);
  histS14->GetYaxis()->SetLabelSize(0.045);
  histS14->GetYaxis()->SetMaxDigits(2);
  histS14->SetTitle("");
  integral = histS14->Integral();
  histS14->Scale(1./integral);
  max = histS14->GetBinContent(histS14->GetMaximumBin());
  histS14->GetXaxis()->SetRangeUser(0,1200);
  histS14->GetYaxis()->SetRangeUser(0,max+0.1*max);
  histS14->Draw();
  
  can->cd(3);
  gPad->SetGrid(1,1);
  histGAGG->SetStats(0);
  histGAGG->GetXaxis()->SetTitle("charge [P.E.]");
  histGAGG->GetYaxis()->SetTitle("counts normalized");
  histGAGG->GetXaxis()->SetTitleSize(0.045);
  histGAGG->GetYaxis()->SetTitleSize(0.045);
  histGAGG->GetXaxis()->SetLabelSize(0.045);
  histGAGG->GetYaxis()->SetLabelSize(0.045);
  histGAGG->GetYaxis()->SetMaxDigits(2);
  histGAGG->SetTitle("");
  integral = histGAGG->Integral();
  histGAGG->Scale(1./integral);
  max = histGAGG->GetBinContent(histGAGG->GetMaximumBin());
  histGAGG->GetXaxis()->SetRangeUser(0,1000);
  histGAGG->GetYaxis()->SetRangeUser(0,max+0.1*max);
  histGAGG->Draw();
  
  can->SaveAs("canSpec.png");
  
  return true;
}


bool EnergyRes(void){
 
  SFEnergyResolution *enS3;
  SFEnergyResolution *enS4;
  SFEnergyResolution *enS5;
  SFEnergyResolution *enS14;
  
  try{
    enS3 = new SFEnergyResolution(3);
    enS4 = new SFEnergyResolution(4);
    enS5 = new SFEnergyResolution(5);
    enS14 = new SFEnergyResolution(14);
  }
  catch(const char *message){
    cout << message << endl;
    return false;
  }
  
  TGraphErrors *grS3 = enS3->GetEnergyResolutionGraph(0);
  TGraphErrors *grS4 = enS4->GetEnergyResolutionGraph(0);
  TGraphErrors *grS5 = enS5->GetEnergyResolutionGraph(0);
  TGraphErrors *grS14 = enS14->GetEnergyResolutionGraph(0);
 
  vector <double> erS3 = enS3->GetEnergyResolution(0);
  vector <double> erS4 = enS4->GetEnergyResolution(0);
  vector <double> erS5 = enS5->GetEnergyResolution(0);
  vector <double> erS14 = enS14->GetEnergyResolution(0);
  
  delete enS3, enS4, enS5, enS14;
  
  TCanvas *can = new TCanvas("can","can",1000,700);
  gPad->SetGrid(1,1);
  
  grS3->SetMarkerStyle(20);
  grS3->SetMarkerSize(1);
  grS3->SetMarkerColor(kRed);
  grS3->SetLineColor(kRed);
  grS3->GetYaxis()->SetRangeUser(0,9);
  grS3->GetYaxis()->SetTitle("Energy Resolution [%]");
  grS3->SetTitle("");
  
  grS4->SetMarkerStyle(21);
  grS4->SetMarkerSize(1);
  grS4->SetMarkerColor(kBlue);
  grS4->SetLineColor(kBlue);
  
  grS5->SetMarkerStyle(22);
  grS5->SetMarkerSize(1.2);
  grS5->SetMarkerColor(kMagenta);
  grS5->SetLineColor(kMagenta);
  
  grS14->SetMarkerStyle(23);
  grS14->SetMarkerSize(1.2);
  grS14->SetMarkerColor(kGreen+3);
  grS14->SetLineColor(kGreen+3);
  
  for(int i=0; i<9; i++){
   grS3->SetPointError(i,2,grS3->GetErrorY(i));
   grS4->SetPointError(i,2,grS4->GetErrorY(i));
   grS5->SetPointError(i,2,grS5->GetErrorY(i));
   grS14->SetPointError(i,2,grS14->GetErrorY(i)); 
  }
  
  TLegend *leg = new TLegend(0.59,0.22,0.89,0.36);
  leg->AddEntry(grS3,"LuAG:Ce (1)","PE");
  leg->AddEntry(grS4,"LuAG:Ce (2)","PE");
  leg->AddEntry(grS5,"LuAG:Ce (1) + coating","PE");
  leg->AddEntry(grS14,"LYSO:Ce","PE");
  
  grS3->Draw("AP");
  grS4->Draw("P");
  grS5->Draw("P");
  grS14->Draw("P");
  leg->Draw();
 
  TLatex text;
  text.SetTextSize(0.030);
  text.SetNDC(true);
  
  text.SetTextColor(kRed);
  text.DrawLatex(0.2,0.3,Form("ER = (%.3f +/- %.3f)%%",erS3[0],erS3[1]));
  text.SetTextColor(kBlue);
  text.DrawLatex(0.2,0.25,Form("ER = (%.3f +/- %.3f)%%",erS4[0],erS4[1]));
  text.SetTextColor(kMagenta);
  text.DrawLatex(0.2,0.20,Form("ER = (%.3f +/- %.3f)%%",erS5[0],erS5[1]));
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.2,0.15,Form("ER = (%.3f +/- %.3f)%%",erS14[0],erS14[1]));
  
  can->SaveAs("canEnRes.png");
  //delete can;
  
  return true;
}