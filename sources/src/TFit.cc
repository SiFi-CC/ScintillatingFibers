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

TFit::TFit(int series_No, TString Fitmode)
{
	SetDetails(series_No,Fitmode);
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

bool TFit::SetDetails(int series_No, TString FitMode){
	
	Reset();
	Fitboarder=70;
	minentry=10;
	
	seriesNo= series_No;
	
	if(FitMode=="TotalSplit") fmode=0;
	else if(FitMode=="ProcessSplit") fmode=1;
	else if(FitMode=="NoSplit") fmode=2;
	
	if(seriesNo<1 || seriesNo >4){
		cout << "this Measurement series are not suited to be analyse with the Template fit" << endl;
		return false;
	}
	data = new SFData(seriesNo);
	Nspectra=data->GetNpoints();
	//Nspectra=1;
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
	for(int i =0;i<9;i++){
		templates_sum.push_back((TH1D*)templates_call[i]->Clone());
		templates_sum[i]->Sumw2();
		templates_sum[i]->Add(templates_pall[i]);	
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
	energypara.reserve(Nspectra);
	chindf.reserve(Nspectra);
	
	Chi2Map.reserve(Nspectra);
	
	resolutiongenerator = new TRandom3();
	
	for(int i=0;i<Nspectra;i++){
		Peaks.push_back(new SFPeakFinder(spectra[i], "511"));
	}
	
	int nx;
	if(fmode==0){ nx=6;}
	else if(fmode==1){ nx=4;}
	else if(fmode==2){ nx=2;}	
	string stringlist[12]={"IBG","PP511","PP1275","C511","C1275","EBG","IBG","PP""C","EBG","IBG","Sim"};
	const char* labels[nx];

	if(fmode==0){
		nx=6;
		GraphicWeights = new TH2D(Form("Weights_snr%i",seriesNo),Form("Weights_snr%i;Position;Contribution",seriesNo),9,5,95,nx,0,nx);
		for(int j=0;j<nx;j++){
			labels[j]=stringlist[j].c_str();
			GraphicWeights->GetYaxis()->SetBinLabel(j+1,labels[j]);
		}
	}
	if(fmode==1){
		nx=4;
		GraphicWeights = new TH2D(Form("Weights_snr%i",seriesNo),Form("Weights_snr%i;Position;Contribution",seriesNo),9,5,95,nx,0,nx);
		for(int j=0;j<nx;j++){
			labels[j]=stringlist[j+6].c_str();
			GraphicWeights->GetYaxis()->SetBinLabel(j+1,labels[j]);
		}
	}
	if(fmode==2){
		nx=2;
		GraphicWeights = new TH2D(Form("Weights_snr%i",seriesNo),Form("Weights_snr%i;Position;Contribution",seriesNo),9,5,95,nx,0,nx);
		for(int j=0;j<nx;j++){
			labels[j]=stringlist[j+10].c_str();
			GraphicWeights->GetYaxis()->SetBinLabel(j+1,labels[j]);
		}
	}
	
	const char* energylabel[3]={"a","b","c"};
	EnergyConstants = new TH2D(Form("EnergyConstant_snr%i",seriesNo),Form("EnergyConstant_snr%i;Position;Parameter",seriesNo),9,5,95,3,0,3);
	for(int j=0;j<3;j++){
		
		EnergyConstants->GetYaxis()->SetBinLabel(j+1,energylabel[j]);
	}
	
	ec_a= new TGraph(Nspectra);
	ec_a->SetName("EC_a");
	ec_b= new TGraph(Nspectra);
	ec_b->SetName("EC_b");
	ec_c= new TGraph(Nspectra);
	ec_c->SetName("EC_c");
	chigraph= new TGraph(Nspectra);
	chigraph->SetName("Chi/ndf");
	
	g_pp_511= new TGraphErrors(Nspectra);
	g_pp_511->SetName("PP511_Weights");
	g_pp_1275= new TGraphErrors(Nspectra);
	g_pp_1275->SetName("PP1275_Weights");
	g_c_511= new TGraphErrors(Nspectra);
	g_c_511->SetName("C511_Weights");
	g_c_1275= new TGraphErrors(Nspectra);
	g_c_1275->SetName("C1275_Weights");
	g_ibg= new TGraphErrors(Nspectra);
	g_ibg->SetName("IBG_Weights");
	g_ebg= new TGraphErrors(Nspectra);
	g_ebg->SetName("EBG_Weights");
	g_sum= new TGraphErrors(Nspectra);
	g_sum->SetName("Sum_Weights");
	
	return true;
	
}
//------------------------------------------------------------------
THStack* TFit::FitSingleSpectrumSplit(int position, double* weights, double* enerconst){
	resgenerator =  new TRandom3();
	THStack* thisstack=new THStack(Form("%s_stack",spectra[position]->GetName()),"");

	vector<Double_t> cali = Peaks[position]->GetParameter();
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
	
	
	int steps_par0=50;
	int steps_par1=20;
	int steps_par2=20;
	
	Double_t StartPar[3]={200-(double(steps_par0)/2)*2,0.5-(double(steps_par1)/2*0.01),(cali[0]/511)-(double(steps_par2)/2*0.01)};
	Double_t EnergyPar[3];
	
	Int_t ierflg = 0;
	
	Double_t min_chi= 1e9;
	Double_t old_chi= 2e9;
	Double_t ndf;
	for(int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > Fitboarder){
			ndf=(cur_spec->GetNbinsX()-i-1)*6;
			break;
		}
	}

	double *fParamVal  = new double [6];
	double *fParamErr= new double [6];
	
	int failed=0;
	int suc=0;

	cur_bg= (TH1D*)background->Clone();
	
	for(int a=0;a<steps_par0;a++){
		EnergyPar[0]=StartPar[0]+a*2;
		for(int b=0;b<steps_par1;b++){
			EnergyPar[1]=StartPar[1]+b*0.01;
			for(int c=0;c<steps_par2;c++){
				cur_ebg->Reset();
				cur_c511->Reset();
				cur_pp511->Reset();
				cur_c1275->Reset();
				cur_pp1275->Reset();
				EnergyPar[2]=StartPar[2]+c*0.01;
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
				ierflg = SetTMiniut(fParamVal,fParamErr,min_chi,0);
				if(ierflg==0) suc++;
				else failed++;
				if(min_chi<old_chi && ierflg==0){
					cout << "New minimum " << endl;
					ebg[position]= (TH1D*)cur_ebg->Clone();
					pp511[position]= (TH1D*)cur_pp511->Clone();
					pp1275[position]= (TH1D*)cur_pp1275->Clone();
					compton511[position]= (TH1D*)cur_c511->Clone();
					compton1275[position]= (TH1D*)cur_c1275->Clone();
				
					old_chi= min_chi;
					
					enerconst[0]=EnergyPar[0];
					enerconst[1]=EnergyPar[1];
					enerconst[2]=EnergyPar[2];
				
				}
				
			}
		}
	}
	
	chindf.push_back(old_chi/ndf);

	cur_bg->Scale(fParamVal[0]);
	pp511[position]->Scale(fParamVal[1]);
	pp1275[position]->Scale(fParamVal[2]);
	compton511[position]->Scale(fParamVal[3]);
	compton1275[position]->Scale(fParamVal[4]);
	ebg[position]->Scale(fParamVal[5]);
	//~ 
	cur_bg->SetFillColor(1);
	cur_bg->SetLineColor(1);
	cur_bg->SetMarkerColor(1);
	pp511[position]->SetLineColor(2);
	pp511[position]->SetFillColor(2);
	pp511[position]->SetMarkerColor(2);
	pp1275[position]->SetLineColor(3);
	pp1275[position]->SetFillColor(3);
	pp1275[position]->SetMarkerColor(3);
	compton511[position]->SetLineColor(4);
	compton511[position]->SetFillColor(4);
	compton511[position]->SetMarkerColor(4);
	compton1275[position]->SetLineColor(5);
	compton1275[position]->SetFillColor(5);
	compton1275[position]->SetMarkerColor(5);
	ebg[position]->SetLineColor(8);
	ebg[position]->SetFillColor(8);
	ebg[position]->SetMarkerColor(8);
		

	thisstack->Add(cur_bg);
	thisstack->Add(pp511[position]);
	thisstack->Add(pp1275[position]);
	thisstack->Add(compton511[position]);
	thisstack->Add(compton1275[position]);
	thisstack->Add(ebg[position]);
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << enerconst[0] <<endl << enerconst[1] <<endl << enerconst[2] << endl;
	cout << "The fitting weights of the histograms are: bg" << fParamVal[0] << " , p511 " << fParamVal[1]<< " , pp1275 " << fParamVal[2]<< " , c511 " << fParamVal[3]<< " , c1275 " << fParamVal[4] << " , else " << fParamVal[5] << endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	thisresidues->SetName(Form("Residues_pos%i",(position+1)*10));
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(pp511[position]);
	tempsum->Add(pp1275[position]);
	tempsum->Add(compton511[position]);
	tempsum->Add(compton1275[position]);
	tempsum->Add(ebg[position]);
	
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > Fitboarder && tempsum->GetBinContent(i)!=0 && cur_spec->GetBinContent(i)>10){
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
	double total_sum_error;
	double total_sum=tempsum->IntegralAndError(start_bin,tempsum->GetNbinsX(),total_sum_error);

	weights[0]=cur_bg->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[6])/total_sum;
	weights[1]=pp511[position]->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[7])/total_sum;
	weights[2]=pp1275[position]->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[8])/total_sum;
	weights[3]=compton511[position]->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[9])/total_sum;
	weights[4]=compton1275[position]->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[10])/total_sum;
	weights[5]=ebg[position]->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[11])/total_sum;
	weights[6]=weights[0]*sqrt((weights[6]*weights[6]/cur_bg->Integral(start_bin,tempsum->GetNbinsX())/cur_bg->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));
	weights[7]=weights[1]*sqrt((weights[7]*weights[7]/pp511[position]->Integral(start_bin,tempsum->GetNbinsX())/pp511[position]->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));		
	weights[8]=weights[2]*sqrt((weights[8]*weights[8]/pp1275[position]->Integral(start_bin,tempsum->GetNbinsX())/pp1275[position]->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));	
	weights[9]=weights[3]*sqrt((weights[9]*weights[9]/compton511[position]->Integral(start_bin,tempsum->GetNbinsX())/compton511[position]->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));
	weights[10]=weights[4]*sqrt((weights[10]*weights[10]/compton1275[position]->Integral(start_bin,tempsum->GetNbinsX())/compton1275[position]->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));
	weights[11]=weights[5]*sqrt((weights[11]*weights[11]/ebg[position]->Integral(start_bin,tempsum->GetNbinsX())/ebg[position]->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));

	cout << "The weights of the histograms are: bg" << weights[0] << " , p511 " << weights[0]<< " , pp1275 " << weights[2]<< " , c511 " << weights[3]<< " , c1275 " << weights[4] << " , else " << weights[5] << endl;
	
	cout << "Failed fits:" << failed << " converged fits: " << suc << endl; 
	cout << "The Chi2 per ndf is " << chindf[position] << endl; 

	
	return thisstack;
}

//------------------------------------------------------------------
THStack* TFit::FitSingleSpectrumTwoSplit(int position, double* weights, double* enerconst){
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
	
	

	Int_t ierflg = 0;

	Double_t min_chi= 1e9;
	Double_t old_chi=2e9;

	double *fParamVal= new double [4];
	double *fParamErr= new double [4];

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
				ierflg= SetTMiniut(fParamVal,fParamErr,min_chi,1);
				if(min_chi<old_chi && ierflg==0){
					cout << "New minimum " << endl;
					ebg[position]= (TH1D*)cur_ebg->Clone();
					pall[position]= (TH1D*)cur_pall->Clone();
					call[position]= (TH1D*)cur_call->Clone();
				
					old_chi= min_chi;
					EParBest[0]=EnergyPar[0];
					EParBest[1]=EnergyPar[1];
					EParBest[2]=EnergyPar[2];
				}
			}
		}
	}

	
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
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << EParBest[0] <<endl << EParBest[1] <<endl << EParBest[2] << endl;
	cout << "The fitting weights of the histograms are: bg" << fParamVal[0] << " , p " << fParamVal[1]<< " , c " << fParamVal[2]<< " , else " << fParamVal[3]<< endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(pall[position]);
	tempsum->Add(call[position]);
	tempsum->Add(ebg[position]);
		
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > Fitboarder && tempsum->GetBinContent(i)!=0 &&cur_spec->GetBinContent(i)!=0){
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
	
	cout << "The weights of the histograms are: bg " << tempweights[0] << " , pall " << tempweights[0]<< " , call " << tempweights[2]<< " , else " << tempweights[3] << endl;
	

	weights=tempweights;
	enerconst=EParBest;
	
	
	return thisstack;
}

//------------------------------------------------------------------
THStack* TFit::FitSSSumEnergy(int position, double*  weight, double* enerconst){
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
THStack* TFit::FitSSSumWeights(int position, double*  weights, double* enerconst){
	resgenerator =  new TRandom3();
	vector<Double_t> cali = Peaks[position]->GetParameter();
	
	THStack* thisstack=new THStack(Form("%s_stacksum",spectra[position]->GetName()),"");
	
	
	cur_spec = (TH1D*)spectra[position]->Clone();
	cur_sum = (TH1D*)spectra[position]->Clone();
	cur_sum->Reset();
	cur_bg = (TH1D*)background->Clone();
	TH1D* fit_sum = (TH1D*)templates_sum[position]->Clone();
	
	cur_bg->SetFillColor(1);
	cur_bg->SetLineColor(1);
	cur_bg->SetMarkerColor(1);
	cur_sum->SetLineColor(3);
	cur_sum->SetFillColor(3);
	cur_sum->SetMarkerColor(3);
	
	int steps_par0=50;
	int steps_par1=20;
	int steps_par2=20;
	
	Double_t StartPar[3]={200-double(steps_par0)/2,0.5-(double(steps_par1)/2*0.02),(cali[0]/511)-(double(steps_par2)/2*0.02)};
	Double_t EnergyPar[3];
		
	Int_t ierflg = 0;

	double *fParamVal=new double[2];
	double *fParamErr= new double[2];
	
	double min_chi=1e9;
	double old_chi=2e9;
	Double_t ndf;
	for(int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > Fitboarder){
			ndf=(cur_spec->GetNbinsX()-i-1)*6;
			break;
		}
	}
	
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

				ierflg=SetTMiniut(fParamVal,fParamErr,min_chi,2);

				if(min_chi < old_chi  && ierflg==0){
					cout << "New minimum" << endl;
					old_chi=min_chi;		
					enerconst[0]=EnergyPar[0];
					enerconst[1]=EnergyPar[1];
					enerconst[2]=EnergyPar[2];
				}
			}
		}
	}

	chindf.push_back(old_chi/ndf);
	
	TH1D* thissum= (TH1D*)cur_spec->Clone();
	thissum->Reset();

	for (int i=0;i<Nbins; i++){
	  for(int j=0;j<fit_sum->GetBinContent(i);j++){
		 thissum->Fill(Energy_model(fit_sum->GetBinCenter(i),enerconst)); 
	  }
	}

	
	cur_bg->Scale(fParamVal[1]);
	thissum->Scale(fParamVal[0]);
	
	cur_bg->SetFillColor(2);
	thissum->SetFillColor(3);
	
	thisstack->Add(cur_bg);
	thisstack->Add(thissum);
	
	cout << "The Parameter of the energy model res^2 =par[0]+par[1]*par[2]*energy are: " << endl << enerconst[0] <<endl << enerconst[1] <<endl << enerconst[2] << endl;
	cout << "The fit weights of the histograms are: sum " << fParamVal[0] << " , bg " << fParamVal[1] << endl;
	
	TGraphErrors* thisresidues = new TGraphErrors(cur_spec->GetNbinsX());
	thisresidues->SetName(Form("Residues_pos%i",(position+1)*10));
	TH1D* tempsum= (TH1D*)cur_bg->Clone();
	tempsum->Add(thissum);
	
		
	double x,y;
	int start_bin=0;
	for( int i=0;i<cur_spec->GetNbinsX();i++){
		if(cur_spec->GetBinCenter(i) > Fitboarder && tempsum->GetBinContent(i)!=0 &&cur_spec->GetBinContent(i)!=0){
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
	
	double total_sum_error;
	double total_sum=tempsum->IntegralAndError(start_bin,tempsum->GetNbinsX(),total_sum_error);

	weights[0]=cur_bg->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[2])/total_sum;
	weights[1]=thissum->IntegralAndError(start_bin,tempsum->GetNbinsX(),weights[3])/total_sum;
	weights[2]=weights[0]*sqrt((weights[2]*weights[2]/cur_bg->Integral(start_bin,tempsum->GetNbinsX())/cur_bg->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));
	weights[3]=weights[0]*sqrt((weights[3]*weights[3]/thissum->Integral(start_bin,tempsum->GetNbinsX())/thissum->Integral(start_bin,tempsum->GetNbinsX()))+(total_sum_error*total_sum_error/total_sum/total_sum));
	
	cout << "The fitting weights of the histograms are: bg " << weights[0] << " , sum " << weights[0] << endl;
	
	return thisstack;
}


//------------------------------------------------------------------
void TFit::FitSpectra(){
	
	//~ THStack* test =FitSingleSpectrumSplit(0,templateweights[0]);
	//~ test->Print();

	for(int i=0;i<Nspectra;i++){
		double* temp;
        double* tempenergy;
		if(fmode==0){
			temp=new double[12];
			tempenergy=new double[6];
			fittedtemplates.push_back(FitSingleSpectrumSplit(i,temp,tempenergy));
		}
		else if(fmode==1){
			temp=new double[4];
			tempenergy=new double[4];
			fittedtemplates.push_back(FitSingleSpectrumTwoSplit(i,temp,tempenergy));
		}
		else if(fmode==2){
			temp=new double[4];
			tempenergy=new double[6];
			
			fittedtemplates.push_back(FitSSSumWeights(i,temp,tempenergy));
		}
		energypara.push_back(tempenergy);
		templateweights.push_back(temp);
		
	}
	int nx;
	if(fmode==0){ nx=6;}
	else if(fmode==1){ nx=4;}
	else if(fmode==2){ nx=2;}
	for(int i=0;i<Nspectra;i++){
		cout << "The procentual weights for position" << (i+1)*10 <<"mm are " << endl;
		
		for(int j=0;j<nx;j++){
			cout << templateweights[i][j]  << endl;
			GraphicWeights->Fill((i+1)*10,j,templateweights[i][j]);
		}
		cout << "The energy parametes for position" << (i+1)*10 <<"mm are " << endl;
		for(int j=0;j<3;j++){
			cout << energypara[i][j]  << endl;
			EnergyConstants->Fill((i+1)*10,j,energypara[i][j]);
		}
		ec_a->SetPoint(i,(i+1)*10,energypara[i][0]);
		ec_b->SetPoint(i,(i+1)*10,energypara[i][1]);
		ec_c->SetPoint(i,(i+1)*10,energypara[i][2]);
		chigraph->SetPoint(i,(i+1)*10,chindf[i]);		
		
		if(fmode==0){	
			g_pp_511->SetPoint(i,(i+1)*10,templateweights[i][1]);
			g_c_511->SetPoint(i,(i+1)*10,templateweights[i][3]);
			g_pp_1275->SetPoint(i,(i+1)*10,templateweights[i][2]);
			g_c_1275->SetPoint(i,(i+1)*10,templateweights[i][4]);
			g_ibg->SetPoint(i,(i+1)*10,templateweights[i][0]);
			g_ebg->SetPoint(i,(i+1)*10,templateweights[i][5]);
			g_pp_511->SetPointError(i,0,templateweights[i][7]);
			g_c_511->SetPointError(i,0,templateweights[i][9]);
			g_pp_1275->SetPointError(i,0,templateweights[i][8]);
			g_c_1275->SetPointError(i,0,templateweights[i][10]);
			g_ibg->SetPointError(i,0,templateweights[i][6]);
			g_ebg->SetPointError(i,0,templateweights[i][11]);
		}
		else if(fmode==2){
			g_ibg->SetPoint(i,(i+1)*10,templateweights[i][0]);
			g_ibg->SetPointError(i,0,templateweights[i][2]);
			g_sum->SetPoint(i,(i+1)*10,templateweights[i][1]);
			g_sum->SetPointError(i,0,templateweights[i][3]);

		}
	}
	
	
}

//------------------------------------------------------------------
int TFit::SetTMiniut(double * fParam,double *fParErr,double &chi ,int fitmode){
	
	TMinuit *ptMinuit = new TMinuit(10);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel();
	// set the user function that calculates chi_square (the value to minimize)
	if(fitmode==0){
		ptMinuit->SetFCN(calc_chi_square);
	}
	else if(fitmode==1){
		ptMinuit->SetFCN(calc_chi_square_twosplit);

	}
	else if(fitmode==2){
		ptMinuit->SetFCN(calc_chi_square_weights);
	}
	
	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	
	int Parleng;

	// Set starting values and step sizes for parameters
	if(fitmode==0){
		Double_t vstart[6] = {0.5, 2 , 2 , 2, 2, 1};
		Double_t step[6] = {0.001 , 0.001 , 0.001 , 0.001, 0.001, 0.001};
		ptMinuit->mnparm(0, "w_bg", vstart[0], step[0], 0,1,ierflg);
		ptMinuit->mnparm(1, "w_p511", vstart[1], step[1], 0,5,ierflg);
		ptMinuit->mnparm(2, "w_p1275", vstart[2], step[2], 0,5,ierflg);
		ptMinuit->mnparm(3, "w_c511", vstart[3], step[3], 0,5,ierflg);
		ptMinuit->mnparm(4, "w_c1275", vstart[4], step[4], 0,5,ierflg);
		ptMinuit->mnparm(5, "w_elsebg", vstart[5], step[5], 0,5,ierflg);		
		Parleng=6;

	}
	else if(fitmode==1){
		Double_t vstart[4] = {0.5, 2 , 2 , 1};
		Double_t step[4] = {0.001 , 0.001 , 0.001, 0.001};
		ptMinuit->mnparm(0, "w_bg", vstart[0], step[0], 0,1,ierflg);
		ptMinuit->mnparm(1, "w_p", vstart[1], step[1], 0,5,ierflg);
		ptMinuit->mnparm(2, "w_c", vstart[2], step[2], 0,5,ierflg);
		ptMinuit->mnparm(3, "w_elsebg", vstart[3], step[3], 0,5,ierflg);		
		Parleng=4;
		
	}
	
	else if (fitmode==2){
		Double_t vstart[2] = {1.5, 0.5};
		Double_t step[2] = {0.001, 0.001};
		ptMinuit->mnparm(0, "weight_sum", vstart[0], step[0], 0,2,ierflg);
		ptMinuit->mnparm(1, "weight_bg", vstart[1], step[1], 0.1,1,ierflg);		
		Parleng=2;
	}
	
			
	arglist[0] = 1000;
	arglist[1] = 1.;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	
	ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

	
	if(amin<chi && ierflg==0){
		for(int i=0;i<Parleng;i++){	
			ptMinuit->GetParameter(i,fParam[i],fParErr[i]);
			//~ cout << "Test " << fParam[i] << endl;
		}
		
	}
	chi= amin;
	delete ptMinuit;
	return ierflg;			
	
}





//------------------------------------------------------------------
//------------------------------------------------------------------
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

TH2D* TFit::GetGraphicWeights(){
	return GraphicWeights;
}

//------------------------------------------------------------------

TGraph* TFit::GetEC_c(){
	return ec_a;
}

//------------------------------------------------------------------

TGraph* TFit::GetEC_a(){
	return ec_b;
}

//------------------------------------------------------------------

TGraph* TFit::GetEC_b(){
	return ec_c;
}

//------------------------------------------------------------------

TGraph* TFit::GetChi(){
	return chigraph;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetIBG_W(){
	return g_ibg;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetEBG_W(){
	return g_ebg;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetPP511_W(){
	return g_pp_511;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetPP1275_W(){
	return g_pp_1275;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetC511_W(){
	return g_c_511;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetC1275_W(){
	return g_c_1275;
}

//------------------------------------------------------------------

TGraphErrors* TFit::GetSum_W(){
	return g_sum;
}

//------------------------------------------------------------------

TH2D* TFit::GetEnergyConstants(){
	return EnergyConstants;
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
    if(cur_spec->GetBinCenter(i)>Fitboarder){
		if(cur_spec->GetBinContent(i)>minentry) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsplit(i,par))/cur_spec->GetBinError(i);
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
    if(cur_spec->GetBinCenter(i)>Fitboarder){
		if(cur_spec->GetBinContent(i)>minentry) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsumenergy(fit_sum->GetBinContent(i),cur_bg->GetBinContent(i),par))/cur_spec->GetBinError(i);
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
    if(cur_spec->GetBinCenter(i)>Fitboarder){
		if(cur_spec->GetBinContent(i)>minentry) delta= (cur_spec->GetBinContent(i)-Templatefit_functionsumweights(cur_sum->GetBinContent(i),cur_bg->GetBinContent(i),par))/cur_spec->GetBinError(i);
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
    if(cur_spec->GetBinCenter(i)>Fitboarder){
		if(cur_spec->GetBinContent(i)>minentry) delta= (cur_spec->GetBinContent(i)-Templatefit_functiontwosplit(i,par))/cur_spec->GetBinError(i);
		chisq += delta*delta;
	}
  }
  f = chisq;
  return;
}


