// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               scifi.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "TimeConst.hh"
#include "TFit.hh"
#include "SFPeakFinder.hh"
#include "TGraphErrors.h"
#include "TCanvas.h"

int main(){
  
  TFile *f = new TFile("test.root","RECREATE");
  
  SFData *data = new SFData(1);
  data->Print();
  
  int n = data->GetNpoints();
  
  TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fT0>0",10);
  //~ TH1D *h2 = data1->GetSpectrum(0,"fPE","ch_0.fT0>0",1);

  SFPeakFinder *peakfin = new SFPeakFinder(h1,"511");
  peakfin->Print();
  TH1D *peak = peakfin->GetPeak();
  f->cd();
  h1->Write();
  peak->Write();
  
  //vector <TH1D*> hh1 = data->GetSpectra(0,"fAmp","");
  //vector <TH1D*> hh2 = data->GetSpectra(1,"fT0","ch_1.fT0!=-100");
  //vector <TH1D*> rr1 = data->GetRatios("log(ch_0.fPE/ch_1.fPE)","ch_0.fT0<590 && ch_0.fPE>0 && ch_1.fT0<590 && ch_1.fPE>0");
  
  //f->cd();
  
  /*for(int i=0; i<n; i++){
   hh1[i]->Write();
   hh2[i]->Write();
   rr1[i]->Write();
  }*/
  
  //~ TFit* firstfit = new TFit(h1,h2);
  //~ 
  //~ TH1D* spectrum = firstfit->GetSpectrum();
  //~ TH1D* background = firstfit->GetFittedBackground();
  //~ 
  //~ TH1D* p511 = firstfit->GetFittedPhotoPeak(511);
  //~ TH1D* p1275 = firstfit->GetFittedPhotoPeak(1275);
  //~ 
  //~ TH1D* cs511 = firstfit->GetFittedCompton(511);
  //~ TH1D* cs1275 = firstfit->GetFittedCompton(1275);
  //~ 
  //~ f->cd();
  //~ spectrum->Write();
  //~ background->Write();
  //~ p511->Write();
  //~ p1275->Write();
  //~ cs511->Write();
  //~ cs1275->Write();
  
  //~ TProfile *p1 = data->GetSignalAverage(1,30,"",20,true);
  //~ TProfile *p2 = data->GetSignalAverage(1,40,"",20,false);
  //~ vector <TProfile*> AverSig_60(9);
  //~ vector <TProfile*> AverSig_120(9);
  //~ vector <TProfile*> AverSig_200(9);
  //~ 
  //~ vector <TimeConst*> TC_60(9);
  //~ vector <TimeConst*> TC_120(9);
  //~ vector <TimeConst*> TC_200(9);
  //~ 
  //~ for(int i=0; i<9; i++){
	//~ AverSig_60[i]= data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>59.5 && ch_0.fPE<60.5",100,true);
	//~ AverSig_120[i]= data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>119.5 && ch_0.fPE<120.5",100,true);
	//~ AverSig_200[i]= data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>199.5 && ch_0.fPE<200.5",100,true);
	//~ }
  //~ 
  //~ vector<TGraphErrors*> FirstDecay(5);
  //~ vector<TGraphErrors*> SecondDecay(5); 
  //~ 
  //~ for(int j=0;j<5;j++){
	//~ FirstDecay[j] = new TGraphErrors(27);
	//~ FirstDecay[j]->SetTitle(";point;Decay Constatn in ns");
	//~ FirstDecay[j]->SetName(Form("FirstDecay_%i_%.1f",(j+1)*10,10.));
	//~ 
	//~ SecondDecay[j] = new TGraphErrors(27);
	//~ SecondDecay[j]->SetTitle(";point;Decay Constatn in ns");
	//~ SecondDecay[j]->SetName(Form("SecondDecay_%i_%.1f",(j+1)*10,10.));
	//~ 
	//~ for(int i=0; i<9; i++){
		//~ TC_60[i]= new TimeConst(AverSig_60[i],"double",(j+1)*10,10);
		//~ TC_120[i]= new TimeConst(AverSig_120[i],"double",(j+1)*10,10);
		//~ TC_200[i]= new TimeConst(AverSig_200[i],"double",(j+1)*10,10);
		//~ FirstDecay[j]->SetPoint(i,i,TC_60[i]->GetFitData()[1]);
		//~ FirstDecay[j]->SetPointError(i,0,TC_60[i]->GetFitDataError()[1]);
		//~ FirstDecay[j]->SetPoint(i+9,i+9,TC_120[i]->GetFitData()[1]);
		//~ FirstDecay[j]->SetPointError(i+9,0,TC_120[i]->GetFitDataError()[1]);
		//~ FirstDecay[j]->SetPoint(i+18,i+18,TC_200[i]->GetFitData()[1]);
		//~ FirstDecay[j]->SetPointError(i+18,0,TC_200[i]->GetFitDataError()[1]);
		//~ SecondDecay[j]->SetPoint(i,i,TC_60[i]->GetFitData()[3]);
		//~ SecondDecay[j]->SetPointError(i,0,TC_60[i]->GetFitDataError()[3]);
		//~ SecondDecay[j]->SetPoint(i+9,i+9,TC_120[i]->GetFitData()[3]);
		//~ SecondDecay[j]->SetPointError(i+9,0,TC_120[i]->GetFitDataError()[3]);
		//~ SecondDecay[j]->SetPoint(i+18,i+18,TC_200[i]->GetFitData()[3]);
		//~ SecondDecay[j]->SetPointError(i+18,0,TC_200[i]->GetFitDataError()[3]);
	//~ }
  //~ }
  //~ 
  //~ 
  //~ f->cd();
//~ 
  //~ for(int j=0; j<5; j++){
	//~ FirstDecay[j]->Write();
	//~ SecondDecay[j]->Write();
//~ }
//~ 
  //~ TH1D *s1 = data->GetSignal(0,60,"",10,true);
  //~ TH1D *s2 = data->GetSignal(0,70,"",14,false);
  //~ TH1D *s3 = data->GetSignal(0,80,"fAmp<100",1,true);
  //~ 
    
  //~ TimeConst* decay_p3_s = new TimeConst(p3,"single");
  //~ decay_p3_s->Print();
  //~ 
  //~ TimeConst* decay_p3_d = new TimeConst(p3,"double");
  //~ decay_p3_d->Print();
  
  
  //----- this part is to see if signals with low 
  //----- cut on PE are negative after bl subtraction
  //~ TCanvas *can = new TCanvas("can","can",1000,1000);
  //~ can->Divide(3,3);
  //~ TH1D *sig[9];
  //~ 
  //~ for(int i=0; i<9; i++){
    //~ sig[i] = data->GetSignal(0,50,"ch_0.fPE>19.99 && ch_0.fPE<20.01",i+1,true);
    //~ can->cd(i+1);
    //~ sig[i]->Draw();
  //~ }
  //~ 
  //~ f->cd();
  //~ can->Write();
  //~ //-----
  //~ 
  //~ f->cd();
  //~ h1->Write();
  //~ h2->Write();
  //~ p1->Write();
  //~ p2->Write();
  //~ p3->Write();
  //~ s1->Write();
  //~ s2->Write();
  //~ s3->Write();
  
  f->Close();

  delete data;
  
  return 1;
}
