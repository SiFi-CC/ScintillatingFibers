// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             lightout.cc               *
// *             Jonas Kasper		   *
// * 	 kasper@physik.rwth-aachen.de 	   *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFLightOutput.hh"
#include "TCanvas.h"
#include "TLatex.h"

int main(int argc, char **argv){
  bool all=false;
  if(argc!=2){
    cout << "to run type: ./lightout seriesNo" << endl;
    return 0;
  }
 
  int seriesNo = atoi(argv[1]);
  if (seriesNo==0) all=true;
 
  if(!all){
	SFData *data; 
	try{
	        data = new SFData(seriesNo);
	}
	catch(const char* message){
	        cout << message << endl;
	        cout << "##### Exception in lightout.cc!" << endl;
	        return 0;
	}

	data->Print();

	SFLightOutput* lout;
	try{
	        lout= new SFLightOutput(seriesNo);
	}
	catch(const char* message){
	        cout << message << endl;
	        cout << "##### Exception in lightout.cc!" << endl;
	        return 0;
	}

	TGraphErrors* loutGraph_ave= lout->GetLightOutputGraph("Ave"); 
	TGraphErrors* loutGraphCh0_ave= lout->GetLightOutputGraph("Ave",0); 
	TGraphErrors* loutGraphCh1_ave= lout->GetLightOutputGraph("Ave",1); 
	vector<double> lightout_ave = lout->GetLightOutput("Ave"); 
	vector<double> lightoutCh0_ave = lout->GetLightOutput("Ave",0); 
	vector<double> lightoutCh1_ave = lout->GetLightOutput("Ave",1); 

	TGraphErrors* loutGraph_sep= lout->GetLightOutputGraph("Sep"); 
	TGraphErrors* loutGraphCh0_sep= lout->GetLightOutputGraph("Sep",0); 
	TGraphErrors* loutGraphCh1_sep= lout->GetLightOutputGraph("Sep",1); 
	vector<double> lightout_sep = lout->GetLightOutput("Sep"); 
	vector<double> lightoutCh0_sep = lout->GetLightOutput("Sep",0); 
	vector<double> lightoutCh1_sep = lout->GetLightOutput("Sep",1); 
	//-----drawing averaged channels
	TLatex text;
	text.SetNDC(true);
	// ----------- averaged attenutation length
	TCanvas *can_averaged_ch_ave = new TCanvas("can_averaged_ch_ave","can_averaged_ch_ave",700,500);
	can_averaged_ch_ave->cd();
	gPad->SetGrid(1,1);
	//loutGraph_ave->SetTitle(Form("Series %i, Ave, lightoutput curve",seriesNo));
	loutGraph_ave->SetTitle("");
	loutGraph_ave->Draw();
	text.SetTextSize(0.04);
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) PH./MeV",lightout_ave[0],lightout_ave[1]));

	TCanvas *can_single_ch0_ave= new TCanvas("can_single_ch0_ave","can_single_ch0_ave",700,500);
	can_single_ch0_ave->cd();
	gPad->SetGrid(1,1);
	//loutGraphCh0_ave->SetTitle(Form("Series %i Ch0, Ave, lightoutput curve",seriesNo));
	loutGraphCh0_ave->SetTitle("");
	loutGraphCh0_ave->Draw();
	text.SetTextSize(0.04);
	//text.DrawLatex(0.2,0.8,Form("LO_Ch0 = (%.2f +/- %.2f) mm",lightoutCh0_ave[0],lightoutCh0_ave[1]));
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightoutCh0_ave[0],lightoutCh0_ave[1]));

	TCanvas *can_single_ch1_ave= new TCanvas("can_single_ch1_ave","can_single_ch1_ave",700,500);
	can_single_ch1_ave->cd();
	gPad->SetGrid(1,1);
	//loutGraphCh1_ave->SetTitle(Form("Series %i Ch1, Ave, lightoutput curve",seriesNo));
	loutGraphCh1_ave->SetTitle("");
	loutGraphCh1_ave->Draw();
	text.SetTextSize(0.04);
	//text.DrawLatex(0.2,0.8,Form("LO_Ch1 = (%.2f +/- %.2f) mm",lightoutCh1_ave[0],lightoutCh1_ave[1]));
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightoutCh1_ave[0],lightoutCh1_ave[1]));
	//----------seperate attenuation length
	TCanvas *can_averaged_ch_sep = new TCanvas("can_averaged_ch_sep","can_averaged_ch_sep",700,500);
	can_averaged_ch_sep->cd();
	gPad->SetGrid(1,1);
	//loutGraph_sep->SetTitle(Form("Series %i, Sep, lightoutput curve",seriesNo));
	loutGraph_sep->SetTitle("");
	loutGraph_sep->Draw();
	text.SetTextSize(0.04);
	//text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) mm",lightout_sep[0],lightout_sep[1]));
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightout_sep[0],lightout_sep[1]));

	TCanvas *can_single_ch0_sep= new TCanvas("can_single_ch0_sep","can_single_ch0_sep",700,500);
	can_single_ch0_sep->cd();
	gPad->SetGrid(1,1);
	//loutGraphCh0_sep->SetTitle(Form("Series %i Ch0, Sep, lightoutput curve",seriesNo));
	loutGraphCh0_sep->SetTitle("");
	loutGraphCh0_sep->Draw();
	text.SetTextSize(0.04);
	//text.DrawLatex(0.2,0.8,Form("LO_Ch0 = (%.2f +/- %.2f) mm",lightoutCh0_sep[0],lightoutCh0_sep[1]));
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightoutCh0_sep[0],lightoutCh0_sep[1]));

	TCanvas *can_single_ch1_sep= new TCanvas("can_single_ch1_sep","can_single_ch1_sep",700,500);
	can_single_ch1_sep->cd();
	gPad->SetGrid(1,1);
	//loutGraphCh1_sep->SetTitle(Form("Series %i Ch1, Sep, lightoutput curve",seriesNo));
	loutGraphCh1_sep->SetTitle("");
	loutGraphCh1_sep->Draw();
	text.SetTextSize(0.04);
	//text.DrawLatex(0.2,0.8,Form("LO_Ch1 = (%.2f +/- %.2f) mm",lightoutCh1_sep[0],lightoutCh1_sep[1]));
	text.DrawLatex(0.2,0.8,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightoutCh1_sep[0],lightoutCh1_sep[1]));
	//----- saving
	TString fname = Form("../results/lighoutput_series_%i.root",seriesNo);
	TFile *file = new TFile(fname,"RECREATE");
	can_averaged_ch_ave->Write();
	can_single_ch0_ave->Write();
	can_single_ch1_ave->Write();
	can_averaged_ch_sep->Write();
	can_single_ch0_sep->Write();
	can_single_ch1_sep->Write();
	file->Close();


	delete data;
	delete lout;
  }
  else {
	vector<SFData*> data; 
	try{
		for(int i=1;i<5;i++){
	      		data.push_back(new SFData(i));
	      	}
	}
	catch(const char* message){
	        cout << message << endl;
	        cout << "##### Exception in lightout.cc!" << endl;
	        return 0;
	}

	for(int i=0;i<4;i++){
		data.at(i)->Print();
	}
	vector<SFLightOutput*> lout;
	try{
		for(int i=1;i<5;i++){
	      		lout.push_back(new SFLightOutput(i));
		}
	}
	catch(const char* message){
	        cout << message << endl;
	        cout << "##### Exception in lightout.cc!" << endl;
	        return 0;
	}

	vector<TGraphErrors*> loutGraph_ave;
	vector<TGraphErrors*> loutGraphCh0_ave;
	vector<TGraphErrors*> loutGraphCh1_ave;
	vector<vector<double>> lightout_ave ;
	vector<vector<double>> lightoutCh0_ave;
	vector<vector<double>> lightoutCh1_ave;

	vector<TGraphErrors*> loutGraph_sep;
	vector<TGraphErrors*> loutGraphCh0_sep;
	vector<TGraphErrors*> loutGraphCh1_sep;
	vector<vector<double>> lightout_sep; 
	vector<vector<double>> lightoutCh0_sep;
	vector<vector<double>> lightoutCh1_sep;

	for(int i=0;i<4;i++){
		loutGraph_ave.push_back(lout.at(i)->GetLightOutputGraph("Ave")); 
		loutGraphCh0_ave.push_back(lout.at(i)->GetLightOutputGraph("Ave",0));
		loutGraphCh1_ave.push_back(lout.at(i)->GetLightOutputGraph("Ave",1));
		lightout_ave.push_back(lout.at(i)->GetLightOutput("Ave")); 
		lightoutCh0_ave.push_back(lout.at(i)->GetLightOutput("Ave",0)); 
		lightoutCh1_ave.push_back(lout.at(i)->GetLightOutput("Ave",1)); 

		loutGraph_sep.push_back(lout.at(i)->GetLightOutputGraph("Sep")); 
		loutGraphCh0_sep.push_back(lout.at(i)->GetLightOutputGraph("Sep",0));
		loutGraphCh1_sep.push_back(lout.at(i)->GetLightOutputGraph("Sep",1));
		lightout_sep.push_back(lout.at(i)->GetLightOutput("Sep")); 
		lightoutCh0_sep.push_back(lout.at(i)->GetLightOutput("Sep",0)); 
		lightoutCh1_sep.push_back(lout.at(i)->GetLightOutput("Sep",1));
	}
	 
	TLatex text;
	text.SetNDC(true);
	// ----------- averaged attenutation length
	TCanvas *can_averaged_ch_ave = new TCanvas("can_averaged_ch_ave","can_averaged_ch_ave",700,500);
	can_averaged_ch_ave->cd();
	gPad->SetGrid(1,1);

	TCanvas *can_single_ch0_ave= new TCanvas("can_single_ch0_ave","can_single_ch0_ave",700,500);
	can_single_ch0_ave->cd();
	gPad->SetGrid(1,1);

	TCanvas *can_single_ch1_ave= new TCanvas("can_single_ch1_ave","can_single_ch1_ave",700,500);
	can_single_ch1_ave->cd();
	gPad->SetGrid(1,1);

	for(int i=0;i<4;i++){	
		can_averaged_ch_ave->cd();
		loutGraph_ave.at(i)->SetTitle("");
		//loutGraph_ave.at(i)->SetTitle("Averaged, lightoutput curve");
		loutGraph_ave.at(i)->SetMarkerStyle(20+i);
		loutGraph_ave.at(i)->SetMarkerColor(i+1);
		loutGraph_ave.at(i)->SetLineColor(i+1);
		if(i==0)loutGraph_ave.at(i)->Draw("AP");
		else loutGraph_ave.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		//text.DrawLatex(0.2,0.85-i*0.05,Form("LO%i = (%.2f +/- %.2f) Ph./MeV",i,lightout_ave.at(i)[0],lightout_ave.at(i)[1]));
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightout_ave.at(i)[0],lightout_ave.at(i)[1]));

		can_single_ch0_ave->cd();
		loutGraphCh0_ave.at(i)->SetTitle("");
		//loutGraphCh0_ave.at(i)->SetTitle("Ch0, Averaged, lightoutput curve");
		loutGraphCh0_ave.at(i)->SetMarkerStyle(20+i);
		loutGraphCh0_ave.at(i)->SetMarkerColor(i+1);
		loutGraphCh0_ave.at(i)->SetLineColor(i+1);
		if(i==0)loutGraphCh0_ave.at(i)->Draw("AP");
		else loutGraphCh0_ave.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO = (%.2f +/- %.2f) Ph./MeV",lightoutCh0_ave.at(i)[0],lightoutCh0_ave.at(i)[1]));
		//text.DrawLatex(0.2,0.85-i*0.05,Form("LO_S%i_Ch0 = (%.2f +/- %.2f) Ph./MeV",i,lightoutCh0_ave.at(i)[0],lightoutCh0_ave.at(i)[1]));

		can_single_ch1_ave->cd();
		//loutGraphCh1_ave.at(i)->SetTitle("Ch1, Averaged, lightoutput curve");
		loutGraphCh1_ave.at(i)->SetTitle("");
		loutGraphCh1_ave.at(i)->SetMarkerStyle(20+i);
		loutGraphCh1_ave.at(i)->SetMarkerColor(i+1);
		loutGraphCh1_ave.at(i)->SetLineColor(i+1);
		if(i==0)loutGraphCh1_ave.at(i)->Draw("AP");
		else loutGraphCh1_ave.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		//text.DrawLatex(0.2,0.85-i*0.05,Form("LO_Ch1 = (%.2f +/- %.2f) mm",lightoutCh1_ave.at(i)[0],lightoutCh1_ave.at(i)[1]));
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO = (%.2f +/- %.2f) ",lightoutCh1_ave.at(i)[0],lightoutCh1_ave.at(i)[1]));
	}
	//----------seperate attenuation length
	TCanvas *can_averaged_ch_sep = new TCanvas("can_averaged_ch_sep","can_averaged_ch_sep",700,500);
	can_averaged_ch_sep->cd();
	gPad->SetGrid(1,1);
	
	TCanvas *can_single_ch0_sep= new TCanvas("can_single_ch0_sep","can_single_ch0_sep",700,500);
	can_single_ch0_sep->cd();
	gPad->SetGrid(1,1);
	
	TCanvas *can_single_ch1_sep= new TCanvas("can_single_ch1_sep","can_single_ch1_sep",700,500);
	can_single_ch1_sep->cd();
	gPad->SetGrid(1,1);

	for(int i=0;i<4;i++){	
		can_averaged_ch_sep->cd();
		//loutGraph_sep.at(i)->SetTitle("Seperate, lightoutput curve");
		loutGraph_sep.at(i)->SetTitle("");
		loutGraph_sep.at(i)->SetMarkerStyle(20+i);
		loutGraph_sep.at(i)->SetMarkerColor(i+1);
		loutGraph_sep.at(i)->SetLineColor(i+1);
		if(i==0)loutGraph_sep.at(i)->Draw("AP");
		else loutGraph_sep.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO_S%i = (%.2f +/- %.2f) mm",i,lightout_sep.at(i)[0],lightout_sep.at(i)[1]));

		can_single_ch0_sep->cd();
		loutGraphCh0_sep.at(i)->SetTitle("Ch0, Seperate, lightoutput curve");
		loutGraphCh0_sep.at(i)->SetMarkerStyle(20+i);
		loutGraphCh0_sep.at(i)->SetMarkerColor(i+1);
		loutGraphCh0_sep.at(i)->SetLineColor(i+1);
		if(i==0)loutGraphCh0_sep.at(i)->Draw("AP");
		else loutGraphCh0_sep.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO_S%i_Ch0 = (%.2f +/- %.2f) mm",i,lightoutCh0_sep.at(i)[0],lightoutCh0_sep.at(i)[1]));

		can_single_ch1_sep->cd();
		loutGraphCh1_sep.at(i)->SetTitle("Ch1, Seperate, lightoutput curve");
		loutGraphCh1_sep.at(i)->SetMarkerStyle(20+i);
		loutGraphCh1_sep.at(i)->SetMarkerColor(i+1);
		loutGraphCh1_sep.at(i)->SetLineColor(i+1);
		if(i==0)loutGraphCh1_sep.at(i)->Draw("AP");
		else loutGraphCh1_sep.at(i)->Draw("P SAME");
		text.SetTextSize(0.04);
		text.SetTextColor(i+1);
		text.DrawLatex(0.2,0.85-i*0.05,Form("LO_S%i_Ch1 = (%.2f +/- %.2f) mm",i,lightoutCh1_sep.at(i)[0],lightoutCh1_sep.at(i)[1]));
	}
	
	TString fname ="../results/lighoutput_all.root";
	TFile *file = new TFile(fname,"RECREATE");
	can_averaged_ch_ave->Write();
	can_single_ch0_ave->Write();
	can_single_ch1_ave->Write();
	can_averaged_ch_sep->Write();
	can_single_ch0_sep->Write();
	can_single_ch1_sep->Write();
	file->Close();

	for(int i=3;i>-1;i--){
		delete data[i];
		delete lout[i];
	}

  }	
  return 1;
} 
