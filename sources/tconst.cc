// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              tconst.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFTimeConst.hh"
#include "SFFitResults.hh"
#include "TCanvas.h"
#include "TLatex.h"

int main(int argc, char **argv){
 
  if(argc!=2){
    cout << "to run type: ./tconst seriesNo" << endl;
    return 0;
  }
  
  int seriesNo = atoi(argv[1]);
  
  SFData *data = new SFData(seriesNo);
  int npoints = data->GetNpoints();
  double *positions = data->GetPositions();
  data->Print();
  
  SFTimeConst *tconst;
  try{
    tconst = new SFTimeConst(seriesNo,100,true); 
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  tconst->Print();
  tconst->FitAllSignals("double");
  vector <TProfile*> signals = tconst->GetAllSignals();
  vector <SFFitResults*> results = tconst->GetAllResults("double");
  vector <TF1*> compFun;
  
  TCanvas *can = new TCanvas("can","can",1200,1200);
  can->Divide(3,3);
  
  for(int i=0; i<npoints; i++){
    can->cd(i+1);
    gPad->SetGrid(1,1);
    signals[i]->Draw();
    results[i]->GetFunction()->Draw("same");
    compFun = results[i]->GetCompFunctions();
    for(int ii=0; ii<compFun.size(); ii++){
      compFun[ii]->Draw("same");
    }
  }
  
  TFile *file = new TFile(Form("../results/timeconst_series%i.root",seriesNo),"RECREATE");
  can->Write();
  file->Close();
  
  delete tconst;
  
  return 1;
}