
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
#include "SFEnergyResolution.hh"
#include "TCanvas.h"
#include "TLatex.h"

int main(int argc, char **argv){
	if(argc!=2){
		cout << "to run type: ./energyres seriesNo" << endl;
		return 0;
	}

	int seriesNo = atoi(argv[1]);

	SFData *data; 
	try{
		data = new SFData(seriesNo);
	}
	catch(const char* message){
		cout << message << endl;
		cout << "##### Exception in energyres.cc!" << endl;
		return 0;
	}

	data->Print();

	SFEnergyResolution* enres;
	try{
		enres= new SFEnergyResolution(seriesNo);
	}
	catch(const char* message){
		cout << message << endl;
		cout << "##### Exception in energyres.cc!" << endl;
		return 0;
	}
	int npoints = data->GetNpoints();
	TGraphErrors* enresGraphCh0= enres->GetEnergyResolutionGraph(0); 
	TGraphErrors* enresGraphCh1= enres->GetEnergyResolutionGraph(1); 
	vector<double> EnergyResCh0= enres->GetEnergyResolution(0); 
	vector<double> EnergyResCh1 = enres->GetEnergyResolution(1); 
	vector<TH1D*> Spec_Ch0= enres->GetSpectra(0);
	vector<TH1D*> Spec_Ch1= enres->GetSpectra(1);


	//-----drawing channels
	TLatex text;
	text.SetNDC(true);
	// ----------- ch 0 
	TCanvas *can_ch_0 = new TCanvas("can_ch_0","can_ch_0",700,500);
	can_ch_0->cd();
	gPad->SetGrid(1,1);
	enresGraphCh0->SetTitle("");
	enresGraphCh0->Draw();
	text.SetTextSize(0.04);
	text.DrawLatex(0.2,0.8,Form("ER = (%.2f +/- %.2f) ",EnergyResCh0[0],EnergyResCh0[1]));

	// ----------- ch 1 
	TCanvas *can_ch_1 = new TCanvas("can_ch_1","can_ch_1",700,500);
	can_ch_1->cd();
	gPad->SetGrid(1,1);
	enresGraphCh1->SetTitle("");
	enresGraphCh1->Draw();
	text.SetTextSize(0.04);
	text.DrawLatex(0.2,0.8,Form("ER = (%.2f +/- %.2f) ",EnergyResCh1[0],EnergyResCh1[1]));
	
	// ----------- spectra 
  	TCanvas *can_spec_ch0 = new TCanvas("can_spec_ch0","can_spec_Ch0",1200,1200);
  	can_spec_ch0->Divide(3,3);
	for(int i=0;i< npoints;i++){
		can_spec_ch0->cd(i+1);
		Spec_Ch0[i]->Draw();
	}
	
  	TCanvas *can_spec_ch1 = new TCanvas("can_spec_ch1","can_spec_Ch1",1200,1200);
  	can_spec_ch1->Divide(3,3);
	for(int i=0;i< npoints;i++){
		can_spec_ch1->cd(i+1);
		Spec_Ch1[i]->Draw();
	}
	//----- saving
	TString fname = Form("../results/energyres_series_%i.root",seriesNo);
	TFile *file = new TFile(fname,"RECREATE");
	can_ch_0->Write();
	can_ch_1->Write();
	can_spec_ch0->Write();
	can_spec_ch1->Write();
    	file->Close();


	delete data;
	delete enres;
	return 1;
} 
