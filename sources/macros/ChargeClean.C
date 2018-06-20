R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFData.hh"
#include "../include/SFPeakFinder.hh"

bool ChargeClean(int series=1){
  
  SFData *data = new SFData(series);
  data->Print();
  int npoints = data->GetNpoints();
  double *positions = data->GetPositions();
  
  //----- loading ratio histograms
  TString selection = "log(sqrt(ch_1.fPE/ch_0.fPE))";
  TString cut = "ch_0.fT0<590 && ch_0.fPE>0 && ch_0.fT0>0 && ch_1.fT0<590 && ch_1.fT0>0 && ch_1.fPE>0";
  vector <TH1D*> hRatio = data->GetCustomHistograms(selection,cut);
  
  //----- fitting Gauss to ratio histograms
  TF1 *gaus = new TF1("gaus","gaus",-100,100);
  double hMean  = 0;
  double hSigma = 0;
  
  for(int i=0; i<npoints; i++){
    hMean  = hRatio[i]->GetMean();
    hSigma = hRatio[i]->GetRMS();
    hRatio[i]->Fit(gaus,"Q","",hMean-hSigma,hMean+hSigma); 
  }
  
  //----- loading charge spectra
  vector <TH1D*> hPECh0cut;
  vector <TH1D*> hPECh1cut;
  vector <TH1D*> hPECh0;
  vector <TH1D*> hPECh1;
  TString PEcutCh0 = "ch_0.fPE>0 && ch_0.fT0>0 && ch_0.fT0<590";
  TString PEcutCh1 = "ch_1.fPE>0 && ch_1.fT0>0 && ch_1.fT0<590";
  TString PEcut;
  double fSigma = 0;
  double fMean  = 0;
  
  for(int i=0; i<npoints; i++){
    fMean = hRatio[i]->GetFunction("gaus")->GetParameter(1);
    fSigma = hRatio[i]->GetFunction("gaus")->GetParameter(2);
    PEcut = Form("log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",fMean-fSigma,fMean+fSigma);
    
    hPECh0cut.push_back(data->GetSpectrum(0,"fPE",PEcut+string("&&")+PEcutCh0+"&&"+PEcutCh1,positions[i]));
    hPECh1cut.push_back(data->GetSpectrum(1,"fPE",PEcut+string("&&")+PEcutCh0+"&&"+PEcutCh1,positions[i]));
    
    hPECh0.push_back(data->GetSpectrum(0,"fPE",PEcutCh0,positions[i]));
    hPECh1.push_back(data->GetSpectrum(1,"fPE",PEcutCh1,positions[i]));
  }
  
  //----- drawing spectra & calculating supression
  TCanvas *canCh0 = new TCanvas("canCh0","canCh0",1200,1200);
  canCh0->Divide(3,3);
  
  TCanvas *canCh1 = new TCanvas("canCh1","chaCh1",1200,1200);
  canCh1->Divide(3,3);
  
  TString string;
  TLatex text;
  text.SetTextSize(0.040);
  text.SetTextFont(42);
  text.SetTextColor(kBlack);
  text.SetNDC(true);
  
  TLine line;
  line.SetLineColor(kRed);
  line.SetLineStyle(9);
  
  SFPeakFinder *peakfinCh0 = new SFPeakFinder();
  SFPeakFinder *peakfinCh1 = new SFPeakFinder();
  double xmin, xmax;
  double max;
  
  double intCh0Peak    = 0;
  double intCh0        = 0;
  double intCh0PeakCut = 0;
  double intCh0Cut     = 0;
  double deltaCh0Peak  = 0;
  double deltaCh0      = 0;
  
  double intCh1PeakCut = 0;
  double intCh1Peak    = 0;
  double intCh1Cut     = 0;
  double intCh1        = 0;
  double deltaCh1Peak  = 0;
  double deltaCh1      = 0;
  
  for(int i=0; i<npoints; i++){
    canCh0->cd(i+1);
    gPad->SetGrid(1,1);
    hPECh0[i]->SetTitle(Form("Charge spectrum, source position %.2f mm",positions[i]));
    hPECh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hPECh0[i]->GetYaxis()->SetTitle("counts");
    hPECh0[i]->GetYaxis()->SetMaxDigits(2);
    hPECh0[i]->SetStats(false);
    peakfinCh0->SetSpectrum(hPECh0[i],"511");
    peakfinCh0->FindPeakRange(xmin,xmax);
    max = hPECh0[i]->GetBinContent(hPECh0[i]->GetMaximumBin());
    hPECh0[i]->Draw();
    hPECh0cut[i]->SetLineColor(kRed);
    hPECh0cut[i]->Draw("same");
    line.DrawLine(xmin,0,xmin,max);
    line.DrawLine(xmax,0,xmax,max);
    intCh0Peak    = hPECh0[i]->Integral(hPECh0[i]->FindBin(xmin),hPECh0[i]->FindBin(xmax));
    intCh0        = hPECh0[i]->Integral(hPECh0[i]->FindBin(0),hPECh0[i]->FindBin(xmin));
    intCh0PeakCut = hPECh0cut[i]->Integral(hPECh0cut[i]->FindBin(xmin),hPECh0cut[i]->FindBin(xmax));
    intCh0Cut     = hPECh0cut[i]->Integral(hPECh0cut[i]->FindBin(0),hPECh0cut[i]->FindBin(xmin));
    deltaCh0      = ((intCh0-intCh0Cut)/intCh0)*100;
    deltaCh0Peak  = ((intCh0Peak-intCh0PeakCut)/intCh0Peak)*100;
    string = Form("#delta_{511keV} = %.2f perc.",deltaCh0Peak);
    text.DrawLatex(0.5,0.60,string);
    string = Form("#delta_{other} = %.2f perc.",deltaCh0);
    text.DrawLatex(0.5,0.55,string);
    hPECh0[i]->GetXaxis()->SetRangeUser(0,700);
    
    canCh1->cd(i+1);
    gPad->SetGrid(1,1);
    hPECh1[i]->SetTitle(Form("Charge spectrum, source position %.2f mm",positions[i]));
    hPECh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
    hPECh1[i]->GetYaxis()->SetTitle("counts");
    hPECh1[i]->GetYaxis()->SetMaxDigits(2);
    hPECh1[i]->SetStats(false);
    peakfinCh1->SetSpectrum(hPECh1[i],"511");
    peakfinCh1->FindPeakRange(xmin,xmax);
    max = hPECh1[i]->GetBinContent(hPECh1[i]->GetMaximumBin());
    hPECh1[i]->Draw();
    hPECh1cut[i]->SetLineColor(kRed);
    hPECh1cut[i]->Draw("same");
    line.DrawLine(xmin,0,xmin,max);
    line.DrawLine(xmax,0,xmax,max);
    intCh1Peak    = hPECh1[i]->Integral(hPECh1[i]->FindBin(xmin),hPECh1[i]->FindBin(xmax));
    intCh1        = hPECh1[i]->Integral(hPECh1[i]->FindBin(0),hPECh1[i]->FindBin(xmin));
    intCh1PeakCut = hPECh1cut[i]->Integral(hPECh1cut[i]->FindBin(xmin),hPECh1cut[i]->FindBin(xmax));
    intCh1Cut     = hPECh1cut[i]->Integral(hPECh1cut[i]->FindBin(0),hPECh1cut[i]->FindBin(xmin));
    deltaCh1      = ((intCh1-intCh1Cut)/intCh1)*100;
    deltaCh1Peak  = ((intCh1Peak-intCh1PeakCut)/intCh1Peak)*100;
    string = Form("#delta_{511keV} = %.2f perc.",deltaCh1Peak);
    text.DrawLatex(0.5,0.60,string);
    string = Form("#delta_{other} = %.2f perc.",deltaCh1);
    text.DrawLatex(0.5,0.55,string);
    hPECh1[i]->GetXaxis()->SetRangeUser(0,700);
  }
  
  TLegend *legCh0 = new TLegend(0.285,0.725,0.877,0.871);
  legCh0->AddEntry(hPECh0[0],"all events","PEL");
  legCh0->AddEntry(hPECh0cut[0],"scattered events rejected","PEL");
  canCh0->cd(1);
  legCh0->Draw();
  
  TLegend *legCh1 = new TLegend(0.285,0.725,0.877,0.871);
  legCh1->AddEntry(hPECh1[0],"all events","PEL");
  legCh1->AddEntry(hPECh1cut[0],"scattered events rejected","PEL");
  canCh1->cd(1);
  legCh1->Draw();
  
  return true;
}