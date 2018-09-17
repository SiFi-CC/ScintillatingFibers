R__LOAD_LIBRARY(../../build/libScintillatingFibers.so)
#include "../include/SFAttenuation.hh"
#include "../include/SFLightOutput.hh"
#include "../include/SFTimingRes.hh"

bool JPET2018(void){
  
  //----- attenuation
  
  SFAttenuation *attS3 = new SFAttenuation(3);	//LuAG:Ce (1)
  attS3->AttAveragedCh();
  TGraphErrors *gAttS3 = attS3->GetAttGraph();
  vector <double> attlenS3 = attS3->GetAttenuation();
 
  SFAttenuation *attS4 = new SFAttenuation(4); 	//LuAG:Ce (2)
  attS4->AttAveragedCh();
  TGraphErrors *gAttS4 = attS4->GetAttGraph();
  vector <double> attlenS4 = attS4->GetAttenuation();
    
  SFAttenuation *attS5 = new SFAttenuation(5);	//with coating (1)
  attS5->AttAveragedCh();
  TGraphErrors *gAttS5 = attS5->GetAttGraph();
  vector <double> attlenS5 = attS5->GetAttenuation();
  
  SFAttenuation *attS14 = new SFAttenuation(14); 	//LYSO
  attS14->AttAveragedCh();
  TGraphErrors *gAttS14 = attS14->GetAttGraph();
  vector <double> attlenS14 = attS14->GetAttenuation();
  
  gAttS3->SetMarkerStyle(20);
  gAttS3->SetMarkerSize(1);
  gAttS3->SetMarkerColor(kRed);
  gAttS3->SetLineColor(kRed);
  
  gAttS4->SetMarkerColor(kBlue);
  gAttS4->SetMarkerStyle(21);
  gAttS4->SetMarkerSize(1);
  gAttS4->SetLineColor(kBlue);
  gAttS4->GetFunction("fpol1")->SetLineColor(kBlue);
  
  gAttS5->SetMarkerStyle(22);
  gAttS5->SetMarkerSize(1.2);
  gAttS5->SetMarkerColor(kMagenta);
  gAttS5->SetLineColor(kMagenta);
  gAttS5->SetTitle("");
  gAttS5->GetFunction("fpol1")->SetLineColor(kMagenta);
  gAttS5->GetYaxis()->SetRangeUser(-0.25,0.30);
  
  gAttS14->SetMarkerColor(kGreen+3);
  gAttS14->SetMarkerStyle(23);
  gAttS14->SetMarkerSize(1.2);
  gAttS14->SetLineColor(kGreen+3);
  gAttS14->GetFunction("fpol1")->SetLineColor(kGreen+3);
  
  TLegend *leg = new TLegend(0.133,0.732,0.434,0.882);
  leg->AddEntry(gAttS3,"LuAG:Ce (1) ","PE");
  leg->AddEntry(gAttS4,"LuAG:Ce (2)", "PE");
  leg->AddEntry(gAttS5,"LuAG:Ce (1) + coating","PE");
  leg->AddEntry(gAttS14,"LYSO:Ce","PE");
  
  TLatex text;
  text.SetTextSize(0.030);
  text.SetNDC(true);
  
  delete attS3;
  delete attS4;
  delete attS5;
  delete attS14;
  
  //----- light output
  
  SFLightOutput *loS3 = new SFLightOutput(3);	//LuAG:Ce (1)
  TGraphErrors *gLOS3 = loS3->GetLightOutputGraph("Ave");
  vector <double> vLOS3 = loS3->GetLightOutput("Ave");
  
  SFLightOutput *loS4 = new SFLightOutput(4);	//LuAG:Ce (2)
  TGraphErrors *gLOS4 = loS4->GetLightOutputGraph("Ave");
  vector <double> vLOS4 = loS4->GetLightOutput("Ave");
  
  SFLightOutput *loS5 = new SFLightOutput(5);	//with coating (1)
  TGraphErrors *gLOS5 = loS5->GetLightOutputGraph("Ave");
  vector <double> vLOS5 = loS5->GetLightOutput("Ave");
  
  SFLightOutput *loS14 = new SFLightOutput(14);	//LYSO
  TGraphErrors *gLOS14 = loS14->GetLightOutputGraph("Ave");
  vector <double> vLOS14 = loS14->GetLightOutput("Ave");
  
  gLOS3->SetMarkerStyle(20);
  gLOS3->SetMarkerSize(1);
  gLOS3->SetMarkerColor(kRed);
  gLOS3->SetLineColor(kRed);
  gLOS3->GetYaxis()->SetRangeUser(3000,8100);
  gLOS3->SetTitle("");
    
  gLOS4->SetMarkerStyle(21);
  gLOS4->SetMarkerSize(1);
  gLOS4->SetMarkerColor(kBlue);
  gLOS4->SetLineColor(kBlue);
  
  gLOS5->SetMarkerStyle(22);
  gLOS5->SetMarkerSize(1.2);
  gLOS5->SetMarkerColor(kMagenta);
  gLOS5->SetLineColor(kMagenta);
  
  gLOS14->SetMarkerStyle(23);
  gLOS14->SetMarkerSize(1.2);
  gLOS14->SetMarkerColor(kGreen+3);
  gLOS14->SetLineColor(kGreen+3);
  gLOS14->GetYaxis()->SetTitle("");
  
  delete loS3;
  delete loS5;
  delete loS4;
  delete loS14;
  
  //----- timing resolution
  
  SFTimingRes *trS3 = new SFTimingRes(3,"ft","with cut");
  vector <double> vTRS3 = trS3->GetTimingResolutions();
  TGraphErrors *gTRS3 = trS3->GetT0Graph();
  
  SFTimingRes *trS4 = new SFTimingRes(4,"ft","with cut");
  vector <double> vTRS4 = trS4->GetTimingResolutions();
  TGraphErrors *gTRS4 = trS4->GetT0Graph();
  
  SFTimingRes *trS5 = new SFTimingRes(5,"ft","with cut");
  vector <double> vTRS5 = trS5->GetTimingResolutions();
  TGraphErrors *gTRS5 = trS5->GetT0Graph();
  
  SFTimingRes *trS14 = new SFTimingRes(14,"ft","with cut");
  vector <double> vTRS14 = trS14->GetTimingResolutions();
  TGraphErrors *gTRS14 = trS14->GetT0Graph();
  
  gTRS3->SetMarkerStyle(20);
  gTRS3->SetMarkerSize(1);
  gTRS3->SetMarkerColor(kRed);
  gTRS3->SetLineColor(kRed);
  gTRS3->SetTitle("");
  gTRS3->GetYaxis()->SetRangeUser(-15,15);
  
  gTRS4->SetMarkerStyle(21);
  gTRS4->SetMarkerSize(1);
  gTRS4->SetMarkerColor(kBlue);
  gTRS4->SetLineColor(kBlue);
  
  gTRS5->SetMarkerStyle(22);
  gTRS5->SetMarkerSize(1.2);
  gTRS5->SetMarkerColor(kMagenta);
  gTRS5->SetLineColor(kMagenta);
  
  gTRS14->SetMarkerStyle(23);
  gTRS14->SetMarkerSize(1.2);
  gTRS14->SetMarkerColor(kGreen+3);
  gTRS14->SetLineColor(kGreen+3);
  
  double sumS3 = 0;
  double sumS5 = 0;
  double sumS4 = 0;
  double sumS14 = 0;
  int ipoints = vTRS3.size();
  
  for(int i=0; i<ipoints; i++){
    sumS3+=vTRS3[i];
    sumS5+=vTRS5[i];
    sumS4+=vTRS4[i];
    sumS14+=vTRS14[i];
  }
  
  double avS3 = sumS3/ipoints;
  double avS5 = sumS5/ipoints;
  double avS4 = sumS4/ipoints;
  double avS14 = sumS14/ipoints;
  double errS3 = TMath::StdDev(&vTRS3[0],&vTRS3[ipoints-1])/sqrt(ipoints);
  double errS5 = TMath::StdDev(&vTRS5[0],&vTRS5[ipoints-1])/sqrt(ipoints);
  double errS4 = TMath::StdDev(&vTRS4[0],&vTRS4[ipoints-1])/sqrt(ipoints);
  double errS14 = TMath::StdDev(&vTRS14[0],&vTRS14[ipoints-1])/sqrt(ipoints);
  
  delete trS3;
  delete trS5;
  delete trS4;
  delete trS14;
  
  //----- drawing 
  
  TCanvas *canAtt = new TCanvas("canAtt","canAtt",1000,700);
 
  gPad->SetGrid(1,1);
  gAttS5->Draw("AP");
  gAttS3->Draw("P");
  gAttS4->Draw("P");
  gAttS14->Draw("P");
  leg->Draw();
  text.SetTextColor(kRed);
  text.DrawLatex(0.14,0.65,Form("L_{att} = (%.2f +/- %.2f) mm",attlenS3[0],attlenS3[1]));
  text.SetTextColor(kBlue);
  text.DrawLatex(0.14,0.60,Form("L_{att} = (%.2f +/- %.2f) mm",attlenS4[0],attlenS4[1]));  
  text.SetTextColor(kMagenta);
  text.DrawLatex(0.14,0.55,Form("L_{att} = (%.2f +/- %.2f) mm",attlenS5[0],attlenS5[1]));
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.14,0.50,Form("L_{att} = (%.2f +/- %.2f) mm",attlenS14[0],attlenS14[1]));
  
 
  TCanvas *canLO = new TCanvas("canLO","canLO",1000,700);
   
  gPad->SetGrid(1,1);
  gLOS3->Draw("AP");
  gLOS5->Draw("P");
  gLOS4->Draw("P");
  gLOS14->Draw("P");
  text.SetTextColor(kRed);
  text.DrawLatex(0.14,0.65,Form("LO = (%.0f +/- %.0f) P.E./MeV",vLOS3[0],vLOS3[1]));
  text.SetTextColor(kBlue);
  text.DrawLatex(0.14,0.60,Form("LO = (%.0f +/- %.0f) P.E./MeV",vLOS4[0],vLOS4[1]));
  text.SetTextColor(kMagenta);
  text.DrawLatex(0.14,0.55,Form("LO = (%.0f +/- %.0f) P.E/MeV",vLOS5[0],vLOS5[1]));
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.14,0.50,Form("LO = (%.0f +/- %.0f) P.E./MeV",vLOS14[0],vLOS14[1]));
  
  
  TCanvas *canTR = new TCanvas("canTR","canTR",1000,700);
  
  gPad->SetGrid(1,1);
  gTRS3->Draw("AP");
  gTRS5->Draw("P");
  gTRS4->Draw("P");
  gTRS14->Draw("P");
  text.SetTextColor(kRed);
  text.DrawLatex(0.14,0.85,Form("TR = (%.3f +/- %.3f) ns",avS3,errS3));
  text.SetTextColor(kBlue);
  text.DrawLatex(0.14,0.80,Form("TR = (%.2f +/- %.2f) ns",avS4,errS4));
  text.SetTextColor(kMagenta);
  text.DrawLatex(0.14,0.75,Form("TR = (%.3f +/- %.3f) ns",avS5,errS5));
  text.SetTextColor(kGreen+3);
  text.DrawLatex(0.14,0.70,Form("TR = (%.4f +/- %.4f) ns",avS14,errS14));
  
  canAtt->SaveAs("canAtt.png");
  canLO->SaveAs("canLO.png");
  canTR->SaveAs("canTR.png");
  
  return true;
}