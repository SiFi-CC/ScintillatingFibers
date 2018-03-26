// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               TFit.cc                 *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "TFit.hh"
#include "SFPeakFinder.hh"

#include "TObjArray.h"
#include "TMinuit.h"
using namespace std;

ClassImp(TFit);

//------------------------------------------------------------------
///Default constructor.
///If this constructer is used the details need to be called with SetDetails().
TFit::TFit():seriesNo(0)
{
 cout << "##### Warning in TFit constructor!" << endl;
 cout << "You are using the default constructor. No Fit possible without defining the measurement spectra." <<endl;
}

//------------------------------------------------------------------
///Constructor that takes the series number of the measurement series that is analysed.
///\param series_No Series number of the measurement series that is analysed. For numbering se SFData class.

TFit::TFit(int series_No)
{
	SetDetails(series_No);
	FitSpectra();
}

//------------------------------------------------------------------
TFit::~TFit(){
	
}

//------------------------------------------------------------------
void TFit::Reset(){
	if(seriesNo>0 && seriesNo<5){
		seriesNo=0;
		Nspectra=0;
		spectra.clear();
		templates_else.clear();
		templates_c511.clear();
		templates_c1275.clear();
		templates_call.clear();
		templates_pall.clear();
		templates_sum.clear();
		pp511.clear();
		pp1275.clear();
		pall.clear();
		compton511.clear();
		compton1275.clear();
		call.clear();
		fittedtemplates.clear();
		templateweights.clear();
		if(background!=NULL) delete background;
		if(data!=NULL) delete data;
		if(template_data!=NULL) delete template_data;
		if(bgdata!=NULL) delete bgdata;
	}
}

//-----------------------------------------------------------------
///Function that is called by the Constructors or by the user after the default constructor. Sets all needed parameters.
///\param series_No Series number of the measurement series that is analysed. For numbering se SFData class.

bool TFit::SetDetails(int series_No){
	Reset();
	seriesNo= series_No;
	if(seriesNo<1 || seriesNo >4){
		cout << "this Measurement series are not suited to be analyse with the Template fit" << endl;
		return false;
	}
	data = new SFData(seriesNo);
	Nspectra=data->GetNpoints();
	spectra=data->GetSpectra(0,"fPE","ch_0.fT0>0 && ch_0.fT0< 590");
	
	Nbins= spectra[0]->GetNbinsX();
	
	if(seriesNo == 1 || seriesNo==3) bgdata= new SFData(7);
	else if(seriesNo == 2 || seriesNo==4) bgdata= new SFData(8);
	background = bgdata->GetSpectrum(0,"fPE","ch_0.fT0>0 && ch_0.fT0< 590",1);
	
	template_data =new SFMC(1);
	templates_else= template_data->GetSpectra("Else","");
	templates_c511= template_data->GetSpectra("Compton","511");
	templates_c1275= template_data->GetSpectra("Compton","1275");
	templates_p511= template_data->GetSpectra("Photon","511");
	templates_p1275= template_data->GetSpectra("Photon","1275");
	templates_call= template_data->GetSpectra("Compton","All");
	templates_pall= template_data->GetSpectra("Photon","All");
	for(int i =0;i<6;i++){
		templates_sum.push_back((TH1D*)templates_call[i]->Clone());
		templates_sum[i]->Sumw2();
		templates_sum[i]->Add(templates_pall[i]);	cout << "Test2" << endl;
		templates_sum[i]->Add(templates_else[i]);
	}
	Peaks.reserve(Nspectra);
	
	pp511.reserve(Nspectra);
	pp1275.reserve(Nspectra);
	pall.reserve(Nspectra);
	compton511.reserve(Nspectra);
	compton1275.reserve(Nspectra);
	call.reserve(Nspectra);
	ebg.reserve(Nspectra);
	
	residual.reserve(Nspectra);
	fittedtemplates.reserve(Nspectra);
	templateweights.reserve(Nspectra);
	
	Chi2Map.reserve(Nspectra);
	
	anglegenerator = new TRandom3();
	resolutiongenerator = new TRandom3();
	comptongenerator =  new TRandom3();
	
	for(int i=0;i<Nspectra;i++){
		Peaks.push_back(new SFPeakFinder(spectra[i], "511"));
	}
	
	return true;
	
}
//------------------------------------------------------------------
THStack* TFit::FitSingleSpectrumSplit(int position, double* &weights){
	resgenerator =  new TRandom3();
	THStack* thisstack=new THStack(Form("%s_stack",spectra[position]->GetName()),"");
	
	//~ Peaks.push_back(new SFPeakFinder(spectra[position], "511"));
	vector<Double_t> cali = Peaks[position]->GetParameter();
	//~ vector<Double_t> cali = Calibration(spectra[position]);
	cur_spec = (TH1D*)spectra[position]->Clone();
	cur_c511 = (TH1D*)spectra[position]->Clone();
	cur_c1275 = (TH1D*)spectra[position]->Clone();
	cur_pp511 = (TH1D*)spectra[position]->Clone();
	cur_pp1275 = (TH1D*)spectra[position]->Clone();
	cur_ebg = (TH1D*)spectra[position]->Clone();
	cur_c511->Reset();
	cur_pp511->Reset();
	cur_c1275->Reset();
	cur_pp1275->Reset();
	cur_ebg->Reset();
	
	cur_c511->SetTitle("c511");
	cur_pp511->SetTitle("p511");
	cur_c1275->SetTitle("c1275");
	cur_pp1275->SetTitle("pp1275");
	cur_ebg->SetTitle("ebg");
	
	
	int steps_par0=20;
	int steps_par1=10;
	int steps_par2=10;
	
	Double_t StartPar[3]={200-(double(steps_par0)/2),0.5-(double(steps_par1)/2*0.02),(cali[0]/511)-(double(steps_par2)/2*0.02)};
	Double_t EnergyPar[3];
	Double_t EParBest[3];
	
	
	TMinuit *ptMinuit = new TMinuit(10);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel();
	// set the user function that calculates chi_square (the value to minimize)
	ptMinuit->SetFCN(calc_chi_square);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[6] = {0.5, 2 , 2 , 2, 2, 1};
	static Double_t step[6] = {0.001 , 0.001 , 0.001 , 0.001, 0.001, 0.001};
	ptMinuit->mnparm(0, "w_bg", vstart[0], step[0], 0,1,ierflg);
	ptMinuit->mnparm(1, "w_p511", vstart[1], step[1], 0,5,ierflg);
	ptMinuit->mnparm(2, "w_p1275", vstart[2], step[2], 0,5,ierflg);
	ptMinuit->mnparm(3, "w_c511", vstart[3], step[3], 0,5,ierflg);
	ptMinuit->mnparm(4, "w_c1275", vstart[4], step[4], 0,5,ierflg);
	ptMinuit->mnparm(5, "w_elsebg", vstart[5], step[5], 0,5,ierflg);

	// Now ready for minimization step
	arglist[0] = 1000;
	arglist[1] = 1.;
	Double_t min_chi= 1e9;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
		
	double fParamVal[6];
	double fParamErr[6];
	
	int failed=0;
	int suc=0;
	
	//~ for(int i=0;i<fpt;i++){
		//~ for(int j=0;j<fpt;j++){
			//~ cur_bg= (TH1D*)background->Clone();
			//~ cur_ebg= CaliandRes(spectra[position],templates_else[position],cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			//~ cur_pp511= PhotoPeak(spectra[position],cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			//~ cur_pp1275= PhotoPeak(spectra[position],(cali[0]+(-fpt/2+i)*2)/511*1275,cali[1]+(-fpt/2+j)*0.5);
			//~ cur_c511= CaliandRes(spectra[position],templates_c511[position],cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			//~ cur_c1275= CaliandRes(spectra[position],templates_c1275[position],cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			//~ ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
			//~ ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
			//~ cout << "The Chi2 for " << cali[0]+(-fpt/2+i)*2 << " , " << cali[1]+(-fpt/2+j)*0.5<< " is " << amin << endl;
			//~ cur_Chi2Map->Fill(cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5,amin);
	cur_bg= (TH1D*)background->Clone();
	
	for(int a=0;a<steps_par0;a++){
		EnergyPar[0]=StartPar[0]+a;
		for(int b=0;b<steps_par1;b++){
			EnergyPar[1]=StartPar[1]+b*0.02;
			for(int c=0;c<steps_par2;c++){
				cur_ebg->Reset();
				cur_c511->Reset();
				cur_pp511->Reset();
				cur_c1275->Reset();
				cur_pp1275->Reset();
				EnergyPar[2]=StartPar[2]+c*0.02;
				cout << "Par[0]: " << EnergyPar[0] << endl;
				cout << "Par[1]: " << EnergyPar[1] << endl;
				cout << "Par[2]: " << EnergyPar[2] << endl;
				for (int i=0;i<templates_else[position]->GetNbinsX(); i++){
				
					for(int j=0;j<templates_else[position]->GetBinContent(i);j++){
						cur_ebg->Fill(Energy_model(templates_else[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_c511[position]->GetBinContent(i);j++){
						cur_c511->Fill(Energy_model(templates_c511[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_c1275[position]->GetBinContent(i);j++){
						cur_c1275->Fill(Energy_model(templates_c1275[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_p1275[position]->GetBinContent(i);j++){
						cur_pp1275->Fill(Energy_model(templates_p1275[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_p511[position]->GetBinContent(i);j++){
						cur_pp511->Fill(Energy_model(templates_p511[position]->GetBinCenter(i),EnergyPar)); 
					}
				}
				ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
				ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);	
				if(ierflg==0) suc++;
				else failed++;
				if(amin<min_chi && ierflg==0){
					cout << "New minimum " << endl;
					ebg[position]= (TH1D*)cur_ebg->Clone();
					pp511[position]= (TH1D*)cur_pp511->Clone();
					pp1275[position]= (TH1D*)cur_pp1275->Clone();
					compton511[position]= (TH1D*)cur_c511->Clone();
					compton1275[position]= (TH1D*)cur_c1275->Clone();
				
					ptMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
					ptMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
					ptMinuit->GetParameter(2,fParamVal[2],fParamErr[2]);
					ptMinuit->GetParameter(3,fParamVal[3],fParamErr[3]);
					ptMinuit->GetParameter(4,fParamVal[4],fParamErr[4]);
					ptMinuit->GetParameter(5,fParamVal[5],fParamErr[5]);
					min_chi= amin;
					EParBest[0]=EnergyPar[0];
					EParBest[1]=EnergyPar[1];
					EParBest[2]=EnergyPar[2];
				}
			}
		}
	}
	//~ Chi2Map.push_back(cur_Chi2Map);

	
	cur_bg->Scale(fParamVal[0]);
	pp511[position]->Scale(fParamVal[1]);
	pp1275[position]->Scale(fParamVal[2]);
	compton511[position]->Scale(fParamVal[3]);
	compton1275[position]->Scale(fParamVal[4]);
	ebg[position]->Scale(fParamVal[5]);
	//~ 
	cur_bg->SetFillColor(1);
	//~ cur_bg->SetLineColor(1);
	cur_bg->SetMarkerColor(1);
	//~ pp511[position]->SetLineColor(2);
	pp511[position]->SetFillColor(2);
	pp511[position]->SetMarkerColor(2);
	//~ pp1275[position]->SetLineColor(3);
	pp1275[position]->SetFillColor(3);
	pp1275[position]->SetMarkerColor(3);
	//~ compton511[position]->SetLineColor(4);
	compton511[position]->SetFillColor(4);
	compton511[position]->SetMarkerColor(4);
	//~ compton1275[position]->SetLineColor(5);
	compton1275[position]->SetFillColor(5);
	compton1275[position]->SetMarkerColor(5);
	//~ ebg[position]->SetLineColor(8);
	ebg[position]->SetFillColor(8);
	ebg[position]->SetMarkerColor(8);
		

	thisstack->Add(cur_bg);
	thisstack->Add(pp511[position]);
	thisstack->Add(pp1275[position]);
	thisstack->Add(compton511[position]);
	thisstack->Add(compton1275[position]);
	thisstack->Add(ebg[position]);
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << EParBest[0] <<endl << EParBest[1] <<endl << EParBest[2] << endl;
	cout << "The fitting weights of the histograms are: bg" << fParamVal[0] << " , p511 " << fParamVal[1]<< " , pp1275 " << fParamVal[2]<< " , c511 " << fParamVal[3]<< " , c1275 " << fParamVal[4] << " , else " << fParamVal[5] << endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(pp511[position]);
	tempsum->Add(pp1275[position]);
	tempsum->Add(compton511[position]);
	tempsum->Add(compton1275[position]);
	tempsum->Add(ebg[position]);
	
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > 70 && tempsum->GetBinContent(i)!=0 && cur_spec->GetBinContent(i)>10){
			if(start_bin==0) start_bin=i;
			x=cur_spec->GetBinCenter(i);
			y=(cur_spec->GetBinContent(i)-tempsum->GetBinContent(i))/cur_spec->GetBinContent(i);
			thisresidues->SetPoint(i,x,y);
			thisresidues->SetPointError(i,0,abs(y)*sqrt(((cur_spec->GetBinError(i)*cur_spec->GetBinError(i))/(cur_spec->GetBinContent(i)*cur_spec->GetBinContent(i))+((tempsum->GetBinError(i)*tempsum->GetBinError(i))/(tempsum->GetBinContent(i)*tempsum->GetBinContent(i))))));
		}
		else{
			thisresidues->SetPoint(i,cur_spec->GetBinCenter(i),0);
			thisresidues->SetPointError(i,0,0);
		}
	}
	
	residual.push_back(thisresidues);
	
	double total_sum=tempsum->Integral(start_bin,tempsum->GetNbinsX());
	double tempweights[6];

	tempweights[0]=cur_bg->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[1]=pp511[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[2]=pp1275[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[3]=compton511[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[4]=compton1275[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[5]=ebg[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;

	cout << "The fitting weights of the histograms are: bg" << tempweights[0] << " , p511 " << tempweights[0]<< " , pp1275 " << tempweights[2]<< " , c511 " << tempweights[3]<< " , c1275 " << tempweights[4] << " , else " << tempweights[5] << endl;
	
	cout << "Failed :" << failed << " success: " << suc << endl; 

	weights=tempweights;
	
	
	return thisstack;
}

//------------------------------------------------------------------
THStack* TFit::FitSingleSpectrumTwoSplit(int position, double* &weights){
	resgenerator =  new TRandom3();
	THStack* thisstack=new THStack(Form("%s_stack",spectra[position]->GetName()),"");
	
	vector<Double_t> cali = Peaks[position]->GetParameter();
	cur_spec = (TH1D*)spectra[position]->Clone();
	cur_call = (TH1D*)spectra[position]->Clone();
	cur_pall = (TH1D*)spectra[position]->Clone();
	cur_ebg = (TH1D*)spectra[position]->Clone();
	cur_call->Reset();
	cur_pall->Reset();
	cur_ebg->Reset();
	
	cur_call->SetTitle("compton");
	cur_pall->SetTitle("p511");
	cur_ebg->SetTitle("ebg");
	
	
	int steps_par0=50;
	int steps_par1=20;
	int steps_par2=20;
	
	Double_t StartPar[3]={200-double(steps_par0)/2,0.5-(double(steps_par1)/2*0.02),(cali[0]/511)-(double(steps_par2)/2*0.02)};
	Double_t EnergyPar[3];
	Double_t EParBest[3];
	
	
	//~ TH2D* cur_Chi2Map = new TH2D(Form("%s_Chi2Map",spectra[position]->GetName()),"Chi2Map;Calibration;Resolution at 511",steps_par1,,,steps_par2,StartPar[3],StartPar[3]+stepspar_2*0.02);
	//~ 
	
	TMinuit *ptMinuit = new TMinuit(10);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel();
	// set the user function that calculates chi_square (the value to minimize)
	ptMinuit->SetFCN(calc_chi_square_twosplit);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[4] = {0.5, 2 , 2 , 1};
	static Double_t step[4] = {0.001 , 0.001 , 0.001, 0.001};
	ptMinuit->mnparm(0, "w_bg", vstart[0], step[0], 0,1,ierflg);
	ptMinuit->mnparm(1, "w_p", vstart[1], step[1], 0,5,ierflg);
	ptMinuit->mnparm(2, "w_c", vstart[2], step[2], 0,5,ierflg);
	ptMinuit->mnparm(3, "w_elsebg", vstart[3], step[3], 0,5,ierflg);

	// Now ready for minimization step
	arglist[0] = 1000;
	arglist[1] = 1.;
	Double_t min_chi= 1e9;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	
	double fParamVal[4];
	double fParamErr[4];

	cur_bg= (TH1D*)background->Clone();
	
	for(int a=0;a<steps_par0;a++){
		EnergyPar[0]=StartPar[0]+a;
		for(int b=0;b<steps_par1;b++){
			EnergyPar[1]=StartPar[1]+b*0.02;
			for(int c=0;c<steps_par2;c++){
				cur_ebg->Reset();
				cur_call->Reset();
				cur_pall->Reset();
				EnergyPar[2]=StartPar[2]+c*0.02;
				cout << "Par[0]: " << EnergyPar[0] << endl;
				cout << "Par[1]: " << EnergyPar[1] << endl;
				cout << "Par[2]: " << EnergyPar[2] << endl;
				for (int i=0;i<templates_else[position]->GetNbinsX(); i++){
				
					for(int j=0;j<templates_else[position]->GetBinContent(i);j++){
						cur_ebg->Fill(Energy_model(templates_else[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_call[position]->GetBinContent(i);j++){
						cur_call->Fill(Energy_model(templates_call[position]->GetBinCenter(i),EnergyPar)); 
					}
					for(int j=0;j<templates_pall[position]->GetBinContent(i);j++){
						cur_pall->Fill(Energy_model(templates_pall[position]->GetBinCenter(i),EnergyPar)); 
					}
				}
				ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
				ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);	
				if(amin<min_chi && ierflg==0){
					cout << "New minimum " << endl;
					ebg[position]= (TH1D*)cur_ebg->Clone();
					pall[position]= (TH1D*)cur_pall->Clone();
					call[position]= (TH1D*)cur_call->Clone();
				
					ptMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
					ptMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
					ptMinuit->GetParameter(2,fParamVal[2],fParamErr[2]);
					ptMinuit->GetParameter(3,fParamVal[3],fParamErr[3]);
					min_chi= amin;
					EParBest[0]=EnergyPar[0];
					EParBest[1]=EnergyPar[1];
					EParBest[2]=EnergyPar[2];
				}
			}
		}
	}
	//~ Chi2Map.push_back(cur_Chi2Map);

	
	cur_bg->Scale(fParamVal[0]);
	pall[position]->Scale(fParamVal[1]);
	call[position]->Scale(fParamVal[2]);
	ebg[position]->Scale(fParamVal[3]);
	//~ 
	cur_bg->SetFillColor(1);
	cur_bg->SetMarkerColor(1);
	pall[position]->SetFillColor(2);
	pall[position]->SetMarkerColor(2);
	call[position]->SetFillColor(3);
	call[position]->SetMarkerColor(3);
	ebg[position]->SetFillColor(4);
	ebg[position]->SetMarkerColor(4);
		

	thisstack->Add(cur_bg);
	thisstack->Add(pall[position]);
	thisstack->Add(call[position]);
	thisstack->Add(ebg[position]);
	//~ weights=fParamVal;
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << EParBest[0] <<endl << EParBest[1] <<endl << EParBest[2] << endl;
	cout << "The weights of the histograms are: bg" << fParamVal[0] << " , p " << fParamVal[1]<< " , c " << fParamVal[2]<< " , else " << fParamVal[3]<< endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(pall[position]);
	tempsum->Add(call[position]);
	tempsum->Add(ebg[position]);
		
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > 70 && tempsum->GetBinContent(i)!=0 &&cur_spec->GetBinContent(i)!=0){
			if(start_bin==0) start_bin=i;
			x=cur_spec->GetBinCenter(i);
			y=(cur_spec->GetBinContent(i)-tempsum->GetBinContent(i))/cur_spec->GetBinContent(i);
			thisresidues->SetPoint(i,x,y);
			thisresidues->SetPointError(i,0,abs(y)*sqrt(((cur_spec->GetBinError(i)*cur_spec->GetBinError(i))/(cur_spec->GetBinContent(i)*cur_spec->GetBinContent(i))+(tempsum->GetBinError(i)*tempsum->GetBinError(i))/(tempsum->GetBinContent(i)*tempsum->GetBinContent(i)))));
		}
		else{
			thisresidues->SetPoint(i,cur_spec->GetBinCenter(i),0);
			thisresidues->SetPointError(i,0,0);
		}
	}
	
	residual.push_back(thisresidues);
	
	double total_sum=tempsum->Integral(start_bin,tempsum->GetNbinsX());
	double tempweights[4];

	tempweights[0]=cur_bg->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[1]=pall[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[2]=call[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[3]=ebg[position]->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	
	cout << "The fitting weights of the histograms are: bg " << tempweights[0] << " , pall " << tempweights[0]<< " , call " << tempweights[2]<< " , else " << tempweights[3] << endl;
	

	weights=tempweights;
	
	
	return thisstack;
}



//------------------------------------------------------------------
THStack* TFit::FitSSSumEnergy(int position, double* & weight){
	resgenerator =  new TRandom3();
	vector<Double_t> cali = Peaks[position]->GetParameter();
	
	THStack* thisstack=new THStack(Form("%s_stacksum",spectra[position]->GetName()),"");
	
	
	cur_spec = (TH1D*)spectra[position]->Clone();
	cur_bg = (TH1D*)background->Clone();
	cur_sum = (TH1D*)templates_sum[position]->Clone();
	
	TMinuit *ptMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel(3);
	// set the user function that calculates chi_square (the value to minimize)
	ptMinuit->SetFCN(calc_chi_square_sum);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	//~ arglist[0] = 2;
	//~ ptMinuit->mnexcm("SET STR", arglist ,1,ierflg);
	// Set starting values and step sizes for parameters
	static Double_t vstart[5] = {200, 0.5, cali[0]/511 , 1.5, 0.5};
	static Double_t step[5] = {0.01 , 0.001 , 0.001 , 0.001, 0.001};
	ptMinuit->mnparm(0, "energy_const", vstart[0], step[0], 100,300,ierflg);
	ptMinuit->mnparm(1, "energy_stat", vstart[1], step[1], 0.1,0.7,ierflg);
	ptMinuit->mnparm(2, "calibration", vstart[2], step[2], 0.4,0.6,ierflg);
	ptMinuit->mnparm(3, "weight_sum", vstart[3], step[3], 1,2,ierflg);
	ptMinuit->mnparm(4, "weight_bg", vstart[4], step[4], 0.3,0.6,ierflg);

	// Now ready for minimization step
	arglist[0] = 1000;
	arglist[1] = 1;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
		
	double fParamVal[5];
	double fParamErr[5];
	
	ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

				
	ptMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
	ptMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
	ptMinuit->GetParameter(2,fParamVal[2],fParamErr[2]);
	ptMinuit->GetParameter(3,fParamVal[3],fParamErr[3]);
	ptMinuit->GetParameter(4,fParamVal[4],fParamErr[4]);

	TH1D* thissum= (TH1D*)cur_spec->Clone();
	thissum->Reset();
	
	//~ double testpar[3];
	//~ 
	//~ testpar[0]=5;
	//~ testpar[1]=3;
	//~ testpar[2]=cali[0]/511;

	for (int i=0;i<Nbins; i++){
	  for(int j=0;j<cur_sum->GetBinContent(i);j++){
		 thissum->Fill(Energy_model(cur_sum->GetBinCenter(i),fParamVal)); 
	  }
	}

	
	cur_bg->Scale(fParamVal[4]);
	thissum->Scale(fParamVal[3]);
	
	cur_bg->SetFillColor(2);
	thissum->SetFillColor(3);
	
	thisstack->Add(cur_bg);
	thisstack->Add(thissum);
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << fParamVal[0] <<endl << fParamVal[1] <<endl << fParamVal[2] << endl;
	cout << "The fit weights of the histograms are: sum " << fParamVal[3] << " , bg " << fParamVal[4] << endl;
	return thisstack;
}


//------------------------------------------------------------------
THStack* TFit::FitSSSumWeights(int position, double* & weights){
	resgenerator =  new TRandom3();
	vector<Double_t> cali = Peaks[position]->GetParameter();
	
	THStack* thisstack=new THStack(Form("%s_stacksum",spectra[position]->GetName()),"");
	
	
	cur_spec = (TH1D*)spectra[position]->Clone();
	cur_sum = (TH1D*)spectra[position]->Clone();
	cur_sum->Reset();
	cur_bg = (TH1D*)background->Clone();
	TH1D* fit_sum = (TH1D*)templates_sum[position]->Clone();
	
	int steps_par0=50;
	int steps_par1=20;
	int steps_par2=20;
	
	Double_t StartPar[3]={200-double(steps_par0)/2,0.5-(double(steps_par1)/2*0.02),(cali[0]/511)-(double(steps_par2)/2*0.02)};
	Double_t EnergyPar[3];
	Double_t EParBest[3];
	
	
	TMinuit *ptMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel();
	// set the user function that calculates chi_square (the value to minimize)
	ptMinuit->SetFCN(calc_chi_square_weights);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	// Set starting values and step sizes for parameters
	static Double_t vstart[2] = {1.5, 0.5};
	static Double_t step[2] = {0.001, 0.001};
	ptMinuit->mnparm(0, "weight_sum", vstart[0], step[0], 0,2,ierflg);
	ptMinuit->mnparm(1, "weight_bg", vstart[1], step[1], 0.1,1,ierflg);

	// Now ready for minimization step
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	
	arglist[0] = 300;
	arglist[1] = 1;
	


	double fParamVal[2];
	double fParamErr[2];
	double min_chi=1e9;
	for(int a=0;a<steps_par0;a++){
		EnergyPar[0]=StartPar[0]+a;
		for(int b=0;b<steps_par1;b++){
			EnergyPar[1]=StartPar[1]+b*0.02;
			for(int c=0;c<steps_par2;c++){
				cur_sum->Reset();
				EnergyPar[2]=StartPar[2]+c*0.02;
				cout << "Par[0]: " << EnergyPar[0] << endl;
				cout << "Par[1]: " << EnergyPar[1] << endl;
				cout << "Par[2]: " << EnergyPar[2] << endl;
				for (int i=0;i<fit_sum->GetNbinsX(); i++){
					for(int j=0;j<fit_sum->GetBinContent(i);j++){
						cur_sum->Fill(Energy_model(fit_sum->GetBinCenter(i),EnergyPar)); 
					}
				}
				//~ arglist[0]=1;
				//~ arglist[1]=1;
				//~ ptMinuit->mnexcm("CLE", arglist ,1,ierflg);
				
				ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
				ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
		
				if(amin < min_chi  && ierflg==0){
					cout << "New minimum" << endl;
					min_chi=amin;		
					ptMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
					ptMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
					EParBest[0]=EnergyPar[0];
					EParBest[1]=EnergyPar[1];
					EParBest[2]=EnergyPar[2];
				}
			}
		}
	}

	TH1D* thissum= (TH1D*)cur_spec->Clone();
	thissum->Reset();

	for (int i=0;i<Nbins; i++){
	  for(int j=0;j<fit_sum->GetBinContent(i);j++){
		 thissum->Fill(Energy_model(fit_sum->GetBinCenter(i),EParBest)); 
	  }
	}

	
	cur_bg->Scale(fParamVal[1]);
	thissum->Scale(fParamVal[0]);
	
	cur_bg->SetFillColor(2);
	thissum->SetFillColor(3);
	
	thisstack->Add(cur_bg);
	thisstack->Add(thissum);
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << EParBest[0] <<endl << EParBest[1] <<endl << EParBest[2] << endl;
	cout << "The fit weights of the histograms are: sum " << fParamVal[0] << " , bg " << fParamVal[1] << endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(thissum);
	
		
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > 70 && tempsum->GetBinContent(i)!=0 &&cur_spec->GetBinContent(i)!=0){
			if(start_bin==0) start_bin=i;
			x=cur_spec->GetBinCenter(i);
			y=(cur_spec->GetBinContent(i)-tempsum->GetBinContent(i))/cur_spec->GetBinContent(i);
			thisresidues->SetPoint(i,x,y);
			thisresidues->SetPointError(i,0,abs(y)*sqrt(((cur_spec->GetBinError(i)*cur_spec->GetBinError(i))/(cur_spec->GetBinContent(i)*cur_spec->GetBinContent(i))+(tempsum->GetBinError(i)*tempsum->GetBinError(i))/(tempsum->GetBinContent(i)*tempsum->GetBinContent(i)))));
		}
		else{
			thisresidues->SetPoint(i,cur_spec->GetBinCenter(i),0);
			thisresidues->SetPointError(i,0,0);
		}
	}
	
	residual.push_back(thisresidues);
	
	double total_sum=tempsum->Integral(start_bin,tempsum->GetNbinsX());
	double tempweights[2];

	tempweights[0]=cur_bg->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	tempweights[1]=thissum->Integral(start_bin,tempsum->GetNbinsX())/total_sum;
	
	cout << "The fitting weights of the histograms are: bg " << tempweights[0] << " , sum " << tempweights[0] << endl;
	
	weights=tempweights;
	
	
	
	return thisstack;
}


//------------------------------------------------------------------
void TFit::FitSpectra(){
	
	THStack* test =FitSingleSpectrumSplit(0,templateweights[0]);
	test->Print();
	fittedtemplates.push_back(test);
	//~ fittedsums.push_back(FitSSSumEnergy(0,templateweights[0]));
	//~ fittedsums.push_back(FitSSSumWeights(0,templateweights[0]));
	
}

//------------------------------------------------------------------

vector<TH1D*> TFit::GetSpectra(){
	if(seriesNo<1 && seriesNo>5) cout << "The spectrum is not properly defined" << endl;
	return spectra;
}

//------------------------------------------------------------------

vector<THStack*> TFit::GetFittedTemplates(){
	return fittedtemplates;
}
//------------------------------------------------------------------

//~ vector<THStack*> TFit::GetFittedSums(){
	//~ return fittedsums;
//~ }
//------------------------------------------------------------------

vector<TH2D*> TFit::GetChi2Map(){
	return Chi2Map;
}

//------------------------------------------------------------------

vector<TGraphErrors*> TFit::GetResiduals(){
	return residual;
}

//------------------------------------------------------------------

TH1D* TFit::GetBackground(){
	if(seriesNo<1 && seriesNo>5) cout << "The background is not properly defined" << endl;
	return background;
}

//------------------------------------------------------------------

Double_t Templatefit_functionsplit(int bin,Double_t *par)
{
 double value= par[0]*cur_bg->GetBinContent(bin)+par[1]*cur_pp511->GetBinContent(bin)+par[2]*cur_pp1275->GetBinContent(bin)+par[3]*cur_c511->GetBinContent(bin)+par[4]*cur_c1275->GetBinContent(bin)+par[5]*cur_ebg->GetBinContent(bin);
 return value;
}

void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //calculate chisquare
  double chisq = 0;
  for (int i=0;i<Nbins; i++) {
    // chi square is the quadratic sum of the distance from the point to the function weighted by its error
    double delta  = 0;
    if(cur_spec->GetBinCenter(i)>70){
		if(cur_spec->GetBinContent(i)>10) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsplit(i,par))/cur_spec->GetBinError(i);
		chisq += delta*delta;
	}
  }
  f = chisq;
  return;
}

//------------------------------------------------------------------

Double_t Templatefit_functionsumenergy(Double_t sum_cont,Double_t bg_cont,Double_t *par)
{
 Double_t value= par[5]*bg_cont+par[4]*sum_cont;
 return value;
}
//------------------------------------------------------------------

Double_t Energy_model(Double_t energy,Double_t *par){
	Double_t ressqr = par[0]+(par[1]*(par[2]*energy));
	Double_t calcenergy= resgenerator->Gaus(par[2]*energy,sqrt(ressqr));
	return calcenergy;
} 
//------------------------------------------------------------------


void calc_chi_square_sum(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Generate cali templates
  TH1D* fit_sum = (TH1D*)cur_spec->Clone();
  //~ TH1D* fit_bg = (TH1D*)cur_bg->Clone();
  
  fit_sum->Reset();
  //~ fit_bg->Reset();
  
  for (int i=0;i<Nbins; i++){
	  for(int j=0;j<cur_sum->GetBinContent(i);j++){
		 fit_sum->Fill(Energy_model(cur_sum->GetBinCenter(i),par)); 
	  }
	  //~ for(int j=0;j<cur_bg->GetBinContent(i);j++){
		 //~ fit_bg->Fill(Energy_model(cur_bg->GetBinCenter(i),par)); 
	  //~ }
  }
  //calculate chisquare
  double chisq = 0;
  for (int i=0;i<Nbins; i++) {
    // chi square is the quadratic sum of the distance from the point to the function weighted by its error
    double delta  = 0;
    if(cur_spec->GetBinCenter(i)>100){
		if(cur_spec->GetBinContent(i)>10) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsumenergy(fit_sum->GetBinContent(i),cur_bg->GetBinContent(i),par))/cur_spec->GetBinError(i);
		chisq += delta*delta;
	}
  }
  f = chisq;
  delete fit_sum;
  //~ delete fit_bg;
  return;
}

//------------------------------------------------------------------

Double_t Templatefit_functionsumweights(Double_t sum_cont,Double_t bg_cont,Double_t *par)
{
 Double_t value= par[1]*bg_cont+par[0]*sum_cont;
 return value;
}

//------------------------------------------------------------------

void calc_chi_square_weights(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Generate cali templates
  
  //calculate chisquare
  double chisq = 0;
  for (int i=0;i<Nbins; i++) {
    // chi square is the quadratic sum of the distance from the point to the function weighted by its error
    double delta  = 0;
    if(cur_spec->GetBinCenter(i)>100){
		if(cur_spec->GetBinContent(i)>10) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsumweights(cur_sum->GetBinContent(i),cur_bg->GetBinContent(i),par))/cur_spec->GetBinError(i);
		chisq += delta*delta;
	}
  }
  f = chisq;
  return;
}

//------------------------------------------------------------------


Double_t Templatefit_functiontwosplit(int bin,Double_t *par)
{
 double value= par[0]*cur_bg->GetBinContent(bin)+par[1]*cur_pall->GetBinContent(bin)+par[2]*cur_call->GetBinContent(bin)+par[3]*cur_ebg->GetBinContent(bin);
 return value;
}
//------------------------------------------------------------------

void calc_chi_square_twosplit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Generate cali templates
  
  //calculate chisquare
  double chisq = 0;
  for (int i=0;i<Nbins; i++) {
    // chi square is the quadratic sum of the distance from the point to the function weighted by its error
    double delta  = 0;
    if(cur_spec->GetBinCenter(i)>70){
		if(cur_spec->GetBinContent(i)>10) delta= (cur_spec->GetBinContent(i)-Templatefit_functiontwosplit(i,par))/cur_spec->GetBinError(i);
		chisq += delta*delta;
	}
  }
  f = chisq;
  return;
}


