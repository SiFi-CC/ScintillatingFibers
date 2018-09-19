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
#include "TLine.h"

int main(int argc, char **argv){
 
  if(argc!=2){
    cout << "to run type: ./timeres seriesNo" << endl;
    return 0;
  }
  
  int seriesNo = atoi(argv[1]);
  
  SFData *data;
  try{
    data = new SFData(seriesNo);
  }
  catch(const char *message){
    cout << message << endl;
    cout << "##### Exception in timeres.cc!" << endl;
    return 0;
  }
  TString desc = data->GetDescription();
  if(!desc.Contains("Regular series")){
    cout << "##### Error in timeres.cc! This is not regular series!" << endl;
    cout << "Series number: " << seriesNo << endl;
    cout << "GetDescription: " << desc << endl;
    return 0;
  }
  
  int npoints = data->GetNpoints();
  vector <double> positions = data->GetPositions();
  data->Print();
  
  TString type = data->GetMeasureType(); 
 
  SFTimingRes *tim_ft;
  SFTimingRes *tim_ft_cut;
  SFTimingRes *tim_cf;
  SFTimingRes *tim_cf_cut;
  
  try{
    tim_ft     = new SFTimingRes(seriesNo,"ft","no cut");
    tim_ft->Print();
    tim_ft_cut = new SFTimingRes(seriesNo,"ft","with cut");
    tim_ft_cut->Print();
    if((type=="Lead")){
	    tim_cf     = new SFTimingRes(seriesNo,"cf","no cut");
	    tim_cf->Print();
	    tim_cf_cut = new SFTimingRes(seriesNo,"cf","with cut");
	    tim_cf_cut->Print();
    }
  }
  catch(const char* message){
    cout << message << endl;
    cout << "##### Exception in timeres.cc!" << endl;
    return 0;
  }
  
  //----- timing resolution FT, no enery cut
  vector <TH1D*>  T0diff_ft  = tim_ft->GetT0Diff();
  vector <double> tres_ft    = tim_ft->GetTimingResolutions();
  vector <double> treserr_ft = tim_ft->GetTimingResErrors();
  vector <TH1D*>  ratio_ft   = tim_ft->GetRatios();
  TGraphErrors    *gr_ft     = tim_ft->GetT0Graph();
  gr_ft->SetName("gr_ft");
  
  //----- timing resolution FT, with energy cut
  vector <TH1D*>  T0diff_ft_cut  = tim_ft_cut->GetT0Diff();
  vector <double> tres_ft_cut    = tim_ft_cut->GetTimingResolutions();
  vector <double> treserr_ft_cut = tim_ft_cut->GetTimingResErrors();
  vector <TH1D*>  spec_ft_cut_0  = tim_ft_cut->GetSpectra(0);
  vector <TH1D*>  spec_ft_cut_1  = tim_ft_cut->GetSpectra(1);
  TGraphErrors    *gr_ft_cut     = tim_ft_cut->GetT0Graph();
  gr_ft_cut->SetName("gr_ft_cut");
  
  //----- drawing
  TCanvas *can_ft = new TCanvas("can_ft","can_ft",1200,1200);
  can_ft->Divide(3,3);
  
  TCanvas *can_ft_cut = new TCanvas("can_ft_cut","cut_ft_cut",1200,1200);
  can_ft_cut->Divide(3,3);
  
  TCanvas *can_cf = new TCanvas("can_cf","can_cf",1200,1200);
  can_cf->Divide(3,3);
  
  TCanvas *can_cf_cut = new TCanvas("can_cf_cut","can_cf_cut",1200,1200);
  can_cf_cut->Divide(3,3);
 
  TCanvas *can_ft_ratio = new TCanvas("can_ft_ratio","can_ft_ratio",1200,1200);
  can_ft_ratio->Divide(3,3);
  
  TCanvas *can_ft_spec_0 = new TCanvas("can_ft_spec_0","can_ft_spec_0",1200,1200);
  can_ft_spec_0->Divide(3,3);
  
  TCanvas *can_ft_spec_1 = new TCanvas("can_ft_spec_1","can_ft_spec_1",1200,1200);
  can_ft_spec_1->Divide(3,3);
  
  TCanvas *can_cf_ratio = new TCanvas("can_cf_ratio","can_cf_ratio",1200,1200);
  can_cf_ratio->Divide(3,3);

  TCanvas *can_cf_spec_0 = new TCanvas("can_cf_spec_0","can_cf_spec_0",1200,1200);
  can_cf_spec_0->Divide(3,3);

  TCanvas *can_cf_spec_1 = new TCanvas("can_cf_spec_1","can_cf_spec_1",1200,1200);
  can_cf_spec_1->Divide(3,3);
  
  TLatex  text;
  TString string;
  TString title;
  text.SetNDC(true);
  text.SetTextSize(0.045);
  
  TLine line;
  line.SetLineColor(kRed);
  line.SetLineStyle(9);
  
  SFPeakFinder *peakfin = new SFPeakFinder();
  
  double mean, sigma;
  double max;
  double xmin, xmax;
  double center, delta;		//changed here for smaller cut
  
  double sum_ft = 0;
  double sum_ft_cut = 0;
  
  for(int i=0; i<npoints; i++){
    title = Form("ch_0.fT0 - ch_1.fT0, series %i, source position %.2f mm",seriesNo,positions[i]);
    
    cout << "testtt" << endl;
    can_ft->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_ft[i],treserr_ft[i]);
    T0diff_ft[i]->SetTitle(title);
    T0diff_ft[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_ft[i]->GetXaxis()->SetRangeUser(-50,50);
    T0diff_ft[i]->Draw();
    text.DrawLatex(0.15,0.8,string);
    sum_ft += tres_ft[i];
    
    can_ft_ratio->cd(i+1);
    gPad->SetGrid(1,1);
    ratio_ft[i]->GetXaxis()->SetRangeUser(-1,1);
    ratio_ft[i]->GetXaxis()->SetTitle("ln(#sqrt{ch1/ch0})");
    ratio_ft[i]->SetTitle(Form("ln(#sqrt{ch1/ch0}), source position %.2f mm",positions[i]));
    max = ratio_ft[i]->GetBinContent(ratio_ft[i]->GetMaximumBin());
    mean = ratio_ft[i]->GetFunction("fun")->GetParameter(1);
    sigma = ratio_ft[i]->GetFunction("fun")->GetParameter(2);
    ratio_ft[i]->Draw();
    line.DrawLine(mean-0.5*sigma,0,mean-0.5*sigma,max);
    line.DrawLine(mean+0.5*sigma,0,mean+0.5*sigma,max);
    
    can_ft_cut->cd(i+1);
    gPad->SetGrid(1,1);
    string = Form("%.2f +/- %.2f ns",tres_ft_cut[i],treserr_ft_cut[i]);
    T0diff_ft_cut[i]->SetTitle(title);
    T0diff_ft_cut[i]->GetXaxis()->SetTitle("time difference [ns]");
    T0diff_ft_cut[i]->GetXaxis()->SetRangeUser(-20,20);
    T0diff_ft_cut[i]->Draw();
    text.DrawLatex(0.15,0.8,string);
    sum_ft_cut += tres_ft_cut[i]; 
    
    can_ft_spec_0->cd(i+1);
    gPad->SetGrid(1,1);
    spec_ft_cut_0[i]->GetXaxis()->SetTitle("charge P.E.");
    spec_ft_cut_0[i]->GetYaxis()->SetTitle("counts");
    spec_ft_cut_0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm",seriesNo,positions[i]));
    spec_ft_cut_0[i]->SetStats(false);
    spec_ft_cut_0[i]->GetYaxis()->SetMaxDigits(2);
    peakfin->SetSpectrum(spec_ft_cut_0[i],"511");
    peakfin->FindPeakRange(xmin,xmax);
    center = xmin+(xmax-xmin)/2.;	//changed here for smaller cut
    delta  = (xmax-xmin)/6; 		//
    max = spec_ft_cut_0[i]->GetBinContent(spec_ft_cut_0[i]->GetMaximumBin());
    spec_ft_cut_0[i]->Draw();
    line.DrawLine(center-delta,0,center-delta,max);	//changed here for smaller cut
    line.DrawLine(center+delta,0,center+delta,max);	//
    spec_ft_cut_0[i]->GetXaxis()->SetRangeUser(0,600);
    
    can_ft_spec_1->cd(i+1);
    gPad->SetGrid(1,1);
    spec_ft_cut_1[i]->GetXaxis()->SetTitle("charge P.E.");
    spec_ft_cut_1[i]->GetYaxis()->SetTitle("counts");
    spec_ft_cut_1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm",seriesNo,positions[i]));
    spec_ft_cut_1[i]->SetStats(false);
    spec_ft_cut_1[i]->GetYaxis()->SetMaxDigits(2);
    peakfin->SetSpectrum(spec_ft_cut_1[i],"511");
    peakfin->FindPeakRange(xmin,xmax);
    center = xmin+(xmax-xmin)/2.;	//changed here for smaller cut
    delta  = (xmax-xmin)/6; 		//
    max = spec_ft_cut_1[i]->GetBinContent(spec_ft_cut_1[i]->GetMaximumBin());
    spec_ft_cut_1[i]->Draw();
    line.DrawLine(center-delta,0,center-delta,max);		//changed here for smaller cut
    line.DrawLine(center+delta,0,center+delta,max);		//
    spec_ft_cut_1[i]->GetXaxis()->SetRangeUser(0,600);
    cout << "testtt" << endl;
  }
  
  TCanvas *can_gr_ft = new TCanvas("can_gr_ft","can_gr_ft",1000,500);
  can_gr_ft->Divide(2,1);
  can_gr_ft->cd(1);
  gPad->SetGrid(1,1);
  gr_ft->SetTitle(Form("Series %i, fixed threshold, cut on scattered events",seriesNo));
  gr_ft->Draw("AP");
  can_gr_ft->cd(2);
  gPad->SetGrid(1,1);
  gr_ft_cut->SetTitle(Form("Series %i, fixed threshold, cut on 511 keV peak",seriesNo));
  gr_ft_cut->Draw("AP");
  //-----printing
  cout << "\n\n-------------------------------" << endl;
  cout << "Average timing resolution time for:" << endl;
  cout << "\t Fixed threshold, scattered events rejected: " << sum_ft/npoints << " ns +/- "
       << TMath::StdDev(&tres_ft[0],&tres_ft[npoints-1])/sqrt(npoints) << " ns" << endl;
  cout << "\t Fixed threshold, cut on 511 keV peak: " << sum_ft_cut/npoints << " ns +/- " 
       << TMath::StdDev(&tres_ft_cut[0],&tres_ft_cut[npoints-1])/sqrt(npoints) << " ns" << endl;
  //----- saving
  TString fname = Form("../results/timingres_series%i.root",seriesNo);
  TFile *file = new TFile(fname,"RECREATE");
  can_ft->Write();
  can_ft_ratio->Write();
  can_ft_cut->Write();
  can_ft_spec_0->Write();
  can_ft_spec_1->Write();
  can_gr_ft->Write();
 
  if((type=="Lead")){
  	//----- timing resolution CF, no energy cut
	vector <TH1D*>  T0diff_cf  = tim_cf->GetT0Diff();
	vector <double> tres_cf    = tim_cf->GetTimingResolutions();
	vector <double> treserr_cf = tim_cf->GetTimingResErrors();
	vector <TH1D*>  ratio_cf   = tim_cf->GetRatios();
	TGraphErrors    *gr_cf     = tim_cf->GetT0Graph();
	gr_cf->SetName("gr_cf");

	//----- timing resolution, CF with energy cut
	vector <TH1D*>  T0diff_cf_cut  = tim_cf_cut->GetT0Diff();
	vector <double> tres_cf_cut    = tim_cf_cut->GetTimingResolutions();
	vector <double> treserr_cf_cut = tim_cf_cut->GetTimingResErrors();
	vector <TH1D*>  spec_cf_cut_0  = tim_cf_cut->GetSpectra(0);
	vector <TH1D*>  spec_cf_cut_1  = tim_cf_cut->GetSpectra(1);
	TGraphErrors    *gr_cf_cut     = tim_cf_cut->GetT0Graph(); 
	gr_cf_cut->SetName("gr_cf_cut");

  	double sum_cf = 0;
  	double sum_cf_cut = 0;
	
  	for(int i=0; i<npoints; i++){
		title = Form("ch_0.fT0 - ch_1.fT0, series %i, source position %.2f mm",seriesNo,positions[i]);
		can_cf->cd(i+1);
		gPad->SetGrid(1,1);
		string = Form("%.2f +/- %.2f ns",tres_cf[i],treserr_cf[i]);
		T0diff_cf[i]->SetTitle(title);
		T0diff_cf[i]->GetXaxis()->SetTitle("time difference [ns]");
		T0diff_cf[i]->GetXaxis()->SetRangeUser(-50,50);
		T0diff_cf[i]->Draw();
		text.DrawLatex(0.15,0.8,string);
		sum_cf += tres_cf[i];

		can_cf_ratio->cd(i+1);
		gPad->SetGrid(1,1);
		ratio_cf[i]->GetXaxis()->SetRangeUser(-1,1);
		ratio_cf[i]->GetXaxis()->SetTitle("ln(#sqrt{ch1/ch0})");
		ratio_cf[i]->SetTitle(Form("ln(#sqrt{ch1/ch0}), source position %.2f mm",positions[i]));
		max = ratio_cf[i]->GetBinContent(ratio_ft[i]->GetMaximumBin());
		mean = ratio_cf[i]->GetFunction("fun")->GetParameter(1);
		sigma = ratio_cf[i]->GetFunction("fun")->GetParameter(2);
		ratio_cf[i]->Draw();
		line.DrawLine(mean-0.5*sigma,0,mean-0.5*sigma,max);
		line.DrawLine(mean+0.5*sigma,0,mean+0.5*sigma,max);

		can_cf_cut->cd(i+1);
		gPad->SetGrid(1,1);
		string = Form("%.2f +/- %.2f ns",tres_cf_cut[i],treserr_cf_cut[i]);
		T0diff_cf_cut[i]->SetTitle(title);
		T0diff_cf_cut[i]->GetXaxis()->SetTitle("time difference [ns]");
		T0diff_cf_cut[i]->GetXaxis()->SetRangeUser(-20,20);
		T0diff_cf_cut[i]->Draw();
		text.DrawLatex(0.2,0.8,string);
		sum_cf_cut += tres_cf_cut[i];

		can_cf_spec_0->cd(i+1);
		gPad->SetGrid(1,1);
		spec_cf_cut_0[i]->GetXaxis()->SetTitle("charge P.E.");
		spec_cf_cut_0[i]->GetYaxis()->SetTitle("counts");
		spec_cf_cut_0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm",seriesNo,positions[i]));
		spec_cf_cut_0[i]->SetStats(false);
		spec_cf_cut_0[i]->GetYaxis()->SetMaxDigits(2);
		peakfin->SetSpectrum(spec_cf_cut_0[i],"511");
		peakfin->FindPeakRange(xmin,xmax);
		center = xmin+(xmax-xmin)/2.;	//changed here for smaller cut
		delta  = (xmax-xmin)/6; 		//
		max = spec_cf_cut_0[i]->GetBinContent(spec_cf_cut_0[i]->GetMaximumBin());
		spec_cf_cut_0[i]->Draw();
		line.DrawLine(center-delta,0,center-delta,max);	//changed here for smaller cut
		line.DrawLine(center+delta,0,center+delta,max);	//
		spec_cf_cut_0[i]->GetXaxis()->SetRangeUser(0,600);

		can_cf_spec_1->cd(i+1);
		gPad->SetGrid(1,1);
		spec_cf_cut_1[i]->GetXaxis()->SetTitle("charge P.E.");
		spec_cf_cut_1[i]->GetYaxis()->SetTitle("counts");
		spec_cf_cut_1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm",seriesNo,positions[i]));
		spec_cf_cut_1[i]->SetStats(false);
		spec_cf_cut_1[i]->GetYaxis()->SetMaxDigits(2);
		peakfin->SetSpectrum(spec_cf_cut_1[i],"511");
		peakfin->FindPeakRange(xmin,xmax);
		center = xmin+(xmax-xmin)/2.;	//changed here for smaller cut
		delta  = (xmax-xmin)/6; 		//
		max = spec_cf_cut_1[i]->GetBinContent(spec_cf_cut_1[i]->GetMaximumBin());
		spec_cf_cut_1[i]->Draw();
		line.DrawLine(center-delta,0,center-delta,max);	//changed here for smaller cut
		line.DrawLine(center+delta,0,center+delta,max);	//
		spec_cf_cut_1[i]->GetXaxis()->SetRangeUser(0,600);
        }

	TCanvas *can_gr_cf = new TCanvas("can_gr_cf","can_gr_cf",1000,500);
	can_gr_cf->Divide(2,1);
	can_gr_cf->cd(1);
	gPad->SetGrid(1,1);
	gr_cf->SetTitle(Form("Series %i, constant fraction, cut on scattered events",seriesNo));
	gr_cf->Draw("AP");
	can_gr_cf->cd(2);
	gPad->SetGrid(1,1);
	gr_cf_cut->SetTitle(Form("Series %i, constant fraction, cut on 511 keV peak",seriesNo));
	gr_cf_cut->Draw("AP");
  	
     	cout << "\t Constant fraction, scattered events rejected: " << sum_cf/npoints << " ns +/- "
  	     << TMath::StdDev(&tres_cf[0],&tres_cf[npoints-1])/sqrt(npoints) << " ns" << endl; 
  	cout << "\t Constant fraction, cut on 511 keV peak: " << sum_cf_cut/npoints << " ns +/- " 
  	     << TMath::StdDev(&tres_cf_cut[0],&tres_cf_cut[npoints-1])/sqrt(npoints) << " ns" << endl;
  	cout << "\n\n" << endl;
	can_cf->Write();
	can_cf_ratio->Write();
	can_cf_cut->Write();
	can_cf_spec_0->Write();
	can_cf_spec_1->Write();
	can_gr_cf->Write();
  }
  file->Close();
  
  delete tim_ft;
  delete tim_ft_cut;
  if((type=="Lead")){ 
  	delete tim_cf;
  	delete tim_cf_cut;
  } 
  return 1;
}
