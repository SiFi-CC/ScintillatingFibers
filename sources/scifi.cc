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
#include "SFAttenuation.hh"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"

int main(){
  
  //TString mode="TotalSplit";
  TString mode="NoSplit";

  TFile *f = new TFile("../results/NoSplit_Series4_50_20_20.root","RECREATE");
  
  //SFData *data = new SFData(2);
  //data->Print();
  
  //int n = data->GetNpoints();
  
  //~ SFAttenuation *att = new SFAttenuation(1);
  //~ att->AttSeparateCh(0);
  //~ att->AttSeparateCh(1);
  //~ att->AttAveragedCh();
  
   TFit* test = new TFit(4,mode);
    
   f->cd();
   for(int i=0;i<9;i++){//test->GetSpectra().size();i++){
		test->GetSpectra()[i]->Write();
		test->GetFittedTemplates()[i]->Write();
		test->GetResiduals()[i]->Write();
	}
   test->GetChi()->Write();
   if(mode=="TotalSplit"){
	test->GetGraphicWeights()->Write();
   	test->GetEnergyConstants()->Write();
   	test->GetEC_a()->Write();
   	test->GetEC_b()->Write();
   	test->GetEC_c()->Write();
   	test->GetIBG_W()->Write();
   	test->GetEBG_W()->Write();
   	test->GetPP511_W()->Write();
   	test->GetPP1275_W()->Write();
   	test->GetC511_W()->Write();
   	test->GetC1275_W()->Write();
   }
   else if(mode=="NoSplit"){

   	test->GetEC_a()->Write();
   	test->GetEC_b()->Write();
   	test->GetEC_c()->Write();
   	test->GetIBG_W()->Write();
   	test->GetSum_W()->Write();

   }
  
  //~ TH1D *h1 = data->GetSpectrum(0,"fPE","ch_0.fT0>0",10);
  //~ TH1D *h2 = data1->GetSpectrum(0,"fPE","ch_0.fT0>0",1);
  //~ 
  //~ vector <TH1D*> hh1 = data->GetSpectra(0,"fAmp","");
  //~ vector <TH1D*> hh2 = data->GetSpectra(1,"fT0","ch_1.fT0!=-100");

  //~ for(int i=0; i<n; i++){
    //~ hh1[i]->Write();
    //~ hh2[i]->Write();
  //~ }
  /*std::vector <TH1D*> h1;
  std::vector <TH1D*> h2;
  std::vector <SFPeakFinder*> peakfin;
  
  TCanvas *can = new TCanvas("can","can",1200,1000);
  can->Divide(4,3);
  
  for(int i=0; i<n; i++){
        can->cd(i+1);
	h1.push_back(data->GetSpectrum(0,"fPE","ch_0.fT0>0 && ch_0.fT0<590",(i+1)*10));
	h1[i]->Draw();
	peakfin.push_back(new SFPeakFinder(h1[i],"511"));
	peakfin[i]->Print();
	h2.push_back(peakfin[i]->GetPeak());
	h2[i]->Draw("same");
  }
  
  can->SaveAs("canv.root");
  */
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

  //delete data;
  
  return 1;
}
