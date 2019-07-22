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
#include "TFit.hh"
#include "SFPeakFinder.hh"
#include "SFAttenuation.hh"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"

int main(){
  
  TString mode="TotalSplit";
  //TString mode="NoSplit";

  TFile *f = new TFile("../results/RootFiles/TotalSplit_50_20_20_R4.root","RECREATE");
  
   TFit* TFitResults = new TFit(1,mode,4);
    
   f->cd();
   for(int i=0;i<9;i++){//TFitResults->GetSpectra().size();i++){
		TFitResults->GetSpectra()[i]->Write();
		TFitResults->GetFittedTemplates()[i]->Write();
		TFitResults->GetResiduals()[i]->Write();
	}
   TFitResults->GetChi()->Write();
   if(mode=="TotalSplit"){
	TFitResults->GetGraphicWeights()->Write();
   	TFitResults->GetEC_a()->Write();
   	TFitResults->GetEC_b()->Write();
   	TFitResults->GetEC_c()->Write();
   	TFitResults->GetIBG_W()->Write();
   	TFitResults->GetEBG_W()->Write();
   	TFitResults->GetPP511_W()->Write();
   	TFitResults->GetPP1275_W()->Write();
   	TFitResults->GetC511_W()->Write();
   	TFitResults->GetC1275_W()->Write();
   }
   else if(mode=="NoSplit"){

   	TFitResults->GetEC_a()->Write();
   	TFitResults->GetEC_b()->Write();
   	TFitResults->GetEC_c()->Write();
   	TFitResults->GetIBG_W()->Write();
   	TFitResults->GetSum_W()->Write();

   }
  
  f->Close();

  return 1;
}
