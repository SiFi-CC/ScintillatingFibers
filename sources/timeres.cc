// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              timeres.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************  


#include "SFTimingRes.hh"
#include "TCanvas.h"
#include "TLatex.h"

int main(int argc, char **argv){
 
  if(argc!=2){
    cout << "to run type: ./timeres seriesNo" << endl;
    return 0;
  }
  
  int seriesNo = atoi(argv[1]);
  
  SFData *data = new SFData(seriesNo);
  int npoints = data->GetNpoints();
  double *positions = data->GetPositions();
  data->Print();
  
  SFTimingRes *tim_ft;
  SFTimingRes *tim_ft_cut;
  SFTimingRes *tim_cf;
  SFTimingRes *tim_cf_cut;
  
  try{
    tim_ft     = new SFTimingRes(seriesNo,"ft","no cut");
    tim_ft->Print();
    tim_ft_cut = new SFTimingRes(seriesNo,"ft","with cut");
    tim_ft_cut->Print();
    tim_cf     = new SFTimingRes(seriesNo,"cf","no cut");
    tim_cf->Print();
    tim_cf_cut = new SFTimingRes(seriesNo,"cf","with cut");
    tim_cf_cut->Print();
  }
  catch(const char* message){
    cout << message << endl;
    return 0;
  }
  
  //----- timing resolution FT, no enery cut
  vector <TH1D*>  T0diff_ft  = tim_ft->GetT0Diff();
  vector <double> tres_ft    = tim_ft->GetTimingResolutions();
  vector <double> treserr_ft = tim_ft->GetTimingResErrors();
  
  //----- timing resolution FT, with energy cut
  vector <TH1D*>  T0diff_ft_cut  = tim_ft_cut->GetT0Diff();
  vector <double> tres_ft_cut    = tim_ft_cut->GetTimingResolutions();
  vector <double> treserr_ft_cut = tim_ft_cut->GetTimingResErrors();
  
  //----- timing resolution CF, no energy cut
  vector <TH1D*>  T0diff_cf  = tim_cf->GetT0Diff();
  vector <double> tres_cf    = tim_cf->GetTimingResolutions();
  vector <double> treserr_cf = tim_cf->GetTimingResErrors();
  
  //----- timing resolution, CF with energy cut
  vector <TH1D*>  T0diff_cf_cut  = tim_cf_cut->GetT0Diff();
  vector <double> tres_cf_cut    = tim_cf_cut->GetTimingResolutions();
  vector <double> treserr_cf_cut = tim_cf_cut->GetTimingResErrors();
  
  //----- drawing
  TCanvas *can_ft = new TCanvas("can_ft","can_ft",1200,1200);
  can_ft->Divide(3,3);
  
  TCanvas *can_ft_cut = new TCanvas("can_ft_cut","cut_ft_cut",1200,1200);
  can_ft_cut->Divide(3,3);
  
  TCanvas *can_cf = new TCanvas("can_cf","can_cf",1200,1200);
  can_cf->Divide(3,3);
  
  TCanvas *can_cf_cut = new TCanvas("can_cf_cut","can_cf_cut",1200,1200);
  can_cf_cut->Divide(3,3);
  
  TLatex  text;
  TString string;
  TString title;
  text.SetNDC(true);
  text.SetTextSize(0.045);
  
  for(int i=0; i<npoints; i++){
    title = Form("ch_0.fT0 - ch_1.fT0, series %i, source position %.2f mm",seriesNo,positions[i]);
    
    can_ft->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_ft[i],treserr_ft[i]);
    T0diff_ft[i]->SetTitle(title);
    T0diff_ft[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_ft[i]->GetXaxis()->SetRangeUser(-50,50);
    T0diff_ft[i]->Draw();
    text.DrawLatex(0.15,0.8,string);
    
    can_ft_cut->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_ft_cut[i],treserr_ft_cut[i]);
    T0diff_ft_cut[i]->SetTitle(title);
    T0diff_ft_cut[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_ft_cut[i]->GetXaxis()->SetRangeUser(-20,20);
    T0diff_ft_cut[i]->Draw();
    text.DrawLatex(0.15,0.8,string);
    
    can_cf->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_cf[i],treserr_cf[i]);
    T0diff_cf[i]->SetTitle(title);
    T0diff_cf[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_cf[i]->GetXaxis()->SetRangeUser(-50,50);
    T0diff_cf[i]->Draw();
    text.DrawLatex(0.15,0.8,string);
    
    can_cf_cut->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_cf_cut[i],treserr_cf_cut[i]);
    T0diff_cf_cut[i]->SetTitle(title);
    T0diff_cf_cut[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_cf_cut[i]->GetXaxis()->SetRangeUser(-20,20);
    T0diff_cf_cut[i]->Draw();
    text.DrawLatex(0.2,0.8,string);
  }
  
  //----- saving
  TString fname = Form("../results/timingres_series%i.root",seriesNo);
  TFile *file = new TFile(fname,"RECREATE");
  can_ft->Write();
  can_ft_cut->Write();
  can_cf->Write();
  can_cf_cut->Write();
  file->Close();
  
  delete tim_ft;
  delete tim_ft_cut;
  delete tim_cf;
  delete tim_cf_cut;
  
  return 1;
}