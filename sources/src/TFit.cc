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
#include "BGFit.hh"
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

TFit::TFit(int series_No, int fitpoints)
{
	SetDetails(series_No,fitpoints);
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
		pp511.clear();
		pp1275.clear();
		compton511.clear();
		compton1275.clear();
		fittedtemplates.clear();
		templateweights.clear();
		if(background!=NULL) delete background;
		if(data!=NULL) delete data;
		if(bgdata!=NULL) delete bgdata;
	}
}

//-----------------------------------------------------------------
///Function that is called by the Constructors or by the user after the default constructor. Sets all needed parameters.
///\param series_No Series number of the measurement series that is analysed. For numbering se SFData class.

bool TFit::SetDetails(int series_No, int fitpoints){
	Reset();
	seriesNo= series_No;
	fpt= fitpoints;
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
	background = bgdata->GetSpectrum(0,"fPE","ch_0.fT0>0",1);
	
	pp511.reserve(Nspectra);
	pp1275.reserve(Nspectra);
	compton511.reserve(Nspectra);
	compton1275.reserve(Nspectra);
	
	fittedtemplates.reserve(Nspectra);
	templateweights.reserve(Nspectra);
	
	Chi2Map.reserve(Nspectra);
	
	anglegenerator = new TRandom3();
	resolutiongenerator = new TRandom3();
	comptongenerator =  new TRandom3();
	cbggenerator =  new TRandom3();
	
	return true;
	
}
//------------------------------------------------------------------

Double_t TFit::KNFormular(Double_t PhotoPeakEnergy, Double_t Angle){
	
	Double_t a= 1/(1+((PhotoPeakEnergy/511)*(1-TMath::Cos(Angle*TMath::DegToRad()))));

	Double_t kn = a*a*(a+(1/a)-TMath::Sin(Angle*TMath::DegToRad())*TMath::Sin(Angle*TMath::DegToRad()))/2/511/511/137.036/137.036;
	
	Double_t kn_0 = 1./511./511./137.036/137.036;
	
	return (kn/kn_0);
	
}
//------------------------------------------------------------------

TH1D* TFit::Compton(TH1D* spec,Double_t PhotoPeakEnergy, Double_t califac, Double_t Resolution){
	//~ cout << "Creating the compton spectrum histogram according to an energy of " << PhotoPeakEnergy << " and a resolution of " << Resolution << endl;
	TString histname= Form("cs_%f_%f",PhotoPeakEnergy/511*califac,Resolution);
	TH1D* cs = new TH1D(histname,histname,Nbins,spec->GetBinCenter(1)-spec->GetBinWidth(1)/2,spec->GetBinCenter(spec->GetNbinsX())+spec->GetBinWidth(spec->GetNbinsX())/2);
	

	double angle;
	double prob;
	double energy_e;
	double energy_g;
	for( int i=0; i<spec->GetEntries();i++){
		angle= anglegenerator->Uniform(0,180);
		prob = KNFormular(PhotoPeakEnergy,angle);
		if(comptongenerator->Uniform(0,1)<prob){
			energy_e= PhotoPeakEnergy*(1-1/(1+((PhotoPeakEnergy/511)*(1-TMath::Cos(angle*TMath::DegToRad())))));
			energy_g= PhotoPeakEnergy-energy_e;
			angle= anglegenerator->Uniform(0,180);
			prob = KNFormular(energy_g,angle);
			if(comptongenerator->Uniform(0,1)<prob) energy_e+= energy_g*(1-1/(1+((energy_g/511)*(1-TMath::Cos(angle*TMath::DegToRad())))));
			if(energy_e > 25){
				energy_e=resolutiongenerator->Gaus(energy_e, Resolution)/511*califac;
				cs->Fill(energy_e);
			}
		}
	}
	
	return cs;
}
//------------------------------------------------------------------

TH1D* TFit::PhotoPeak(TH1D* spec, Double_t PhotoPeakEnergy, Double_t Resolution){
	//~ cout << "Creating the photopeak histogram with an energy of " << PhotoPeakEnergy << " and a resolution of " << Resolution << endl;
	TString histname= Form("%f_%f",PhotoPeakEnergy,Resolution);
	TH1D* pp = new TH1D(histname,histname,Nbins,spec->GetBinCenter(1)-spec->GetBinWidth(1)/2,spec->GetBinCenter(spec->GetNbinsX())+spec->GetBinWidth(spec->GetNbinsX())/2);
	
	for( int i=0; i<spec->GetEntries();i++){
		double temp = resolutiongenerator->Gaus(PhotoPeakEnergy,Resolution);
		pp->Fill(temp);
	}
	
	return pp;
}
//------------------------------------------------------------------

TH1D* TFit::ColiBg(TH1D* spec, Double_t Resolution){
	//~ cout << "Creating the photopeak histogram with an energy of " << PhotoPeakEnergy << " and a resolution of " << Resolution << endl;
	TString histname= Form("ColiBg_%f",Resolution);
	TH1D* cbg = new TH1D(histname,histname,Nbins,spec->GetBinCenter(1)-spec->GetBinWidth(1)/2,spec->GetBinCenter(spec->GetNbinsX())+spec->GetBinWidth(spec->GetNbinsX())/2);
	
	for( int i=0; i<spec->GetEntries();i++){
		double temp= spec->GetMaximum() * cbggenerator->Exp(100);
		temp = resolutiongenerator->Gaus(temp,Resolution);
		if (temp > 10) cbg->Fill(temp);
	}
	
	return cbg;
}

//------------------------------------------------------------------

vector <Double_t> TFit::Calibration(TH1D* spectrum){
	TH1D* spec = (TH1D*)spectrum->Clone("Calspec");
	TH1D* bgs= new TH1D("bgs","bgs",spec->GetNbinsX(),spec->GetBinCenter(1),spec->GetBinCenter(spectrum->GetNbinsX()));
  
	BGFit * mybgfit = new BGFit(215,325);
	TF1 *bgfun = new TF1("bgfun",mybgfit,&BGFit::EvaluatePol3,175,450,4,"BGFit","Evaluate");
	
	spec->Fit("bgfun","R");

	Int_t nbins = spec->GetXaxis()->GetNbins();
	Double_t x = 0.;
	Double_t y = 0.;
  
	for(Int_t i=1; i<nbins+1; i++){
		x = spec->GetBinCenter(i);
		if(x>175 && x<450){
			y = spec->GetBinContent(i)-bgfun->Eval(x);
			if(y<0) y=0;
			bgs->SetBinContent(i,y);
		}
		else{
			bgs->SetBinContent(i,0);
		}
	}
	
	TF1* fun_gaus = new TF1("fungaus","gaus",215,325);

	bgs->Fit("fungaus","R");
	
	cout << "511 equals " << fun_gaus->GetParameter(1) << " P.E." << endl;
	cout << "Resolution: " << fun_gaus->GetParameter(2) << " P.E." << endl;
	vector <Double_t> temp;
	temp.push_back(fun_gaus->GetParameter(1));
	temp.push_back(fun_gaus->GetParameter(2));
	return temp;
 
	
}

//------------------------------------------------------------------
THStack* TFit::FitSingleSpectrum(int position, double* &weights){
	THStack* thisstack=new THStack(Form("%s_stack",spectra[position]->GetName()),"");
	
	
	vector<Double_t> cali = Calibration(spectra[position]);
	cur_spec = (TH1D*)spectra[position]->Clone();

	TH2D* cur_Chi2Map = new TH2D(Form("%s_Chi2Map",spectra[position]->GetName()),"Chi2Map;Position;Resolution",fpt,cali[0]-fpt,cali[0]+fpt,fpt,cali[1]-(fpt/2*0.5),cali[1]+(fpt/2*0.5));
	
	
	TMinuit *ptMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	ptMinuit->SetPrintLevel();
	// set the user function that calculates chi_square (the value to minimize)
	ptMinuit->SetFCN(calc_chi_square);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	static Double_t vstart[5] = {0.4, 0.1 , 0.01 , 1, 3};
	static Double_t step[5] = {0.001 , 0.001 , 0.001 , 0.001, 0.001};
	ptMinuit->mnparm(0, "w_bg", vstart[0], step[0], 0,10,ierflg);
	ptMinuit->mnparm(1, "w_p511", vstart[1], step[1], 0,10,ierflg);
	ptMinuit->mnparm(2, "w_p1275", vstart[2], step[2], 0,10,ierflg);
	ptMinuit->mnparm(3, "w_c511", vstart[3], step[3], 0,10,ierflg);
	ptMinuit->mnparm(4, "w_c1275", vstart[4], step[4], 0,10,ierflg);
	//~ ptMinuit->mnparm(5, "w_cbg", vstart[5], step[5], 0,10,ierflg);

	// Now ready for minimization step
	arglist[0] = 1000;
	arglist[1] = 1.;
	Double_t min_chi= 1e9;
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
		
	double fParamVal[5];
	double fParamErr[5];
	
	for(int i=0;i<fpt;i++){
		for(int j=0;j<fpt;j++){
			cur_bg= (TH1D*)background->Clone();
			//~ cur_cbg= ColiBg(spectra[position],cali[1]+(-fpt/2+j)*0.5);
			cur_pp511= PhotoPeak(spectra[position],cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			cur_pp1275= PhotoPeak(spectra[position],(cali[0]+(-fpt/2+i)*2)/511*1275,cali[1]+(-fpt/2+j)*0.5);
			cur_c511= Compton(spectra[position],511,cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5);
			cur_c1275= Compton(spectra[position],1275,(cali[0]+(-fpt/2+i)*2),cali[1]+(-fpt/2+j)*0.5);
			ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
			ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
			cout << "The Chi2 for " << cali[0]+(-fpt/2+i)*2 << " , " << cali[1]+(-fpt/2+j)*0.5<< " is " << amin << endl;
			cur_Chi2Map->Fill(cali[0]+(-fpt/2+i)*2,cali[1]+(-fpt/2+j)*0.5,amin);
			if(amin<min_chi && ierflg==0){
					
				pp511[position]= (TH1D*)cur_pp511->Clone();
				pp1275[position]= (TH1D*)cur_pp1275->Clone();
				compton511[position]= (TH1D*)cur_c511->Clone();
				compton1275[position]= (TH1D*)cur_c1275->Clone();
				
				ptMinuit->GetParameter(0,fParamVal[0],fParamErr[0]);
				ptMinuit->GetParameter(1,fParamVal[1],fParamErr[1]);
				ptMinuit->GetParameter(2,fParamVal[2],fParamErr[2]);
				ptMinuit->GetParameter(3,fParamVal[3],fParamErr[3]);
				ptMinuit->GetParameter(4,fParamVal[4],fParamErr[4]);
				//~ ptMinuit->GetParameter(5,fParamVal[5],fParamErr[5]);
				
			}
		}
	}
	Chi2Map.push_back(cur_Chi2Map);

	
	cur_bg->Scale(fParamVal[0]);
	cur_bg->SetLineColor(1);
	cur_bg->SetMarkerColor(1);
	//~ cur_cbg->Scale(fParamVal[5]);
	//~ cur_cbg->SetLineColor(6);
	//~ cur_cbg->SetMarkerColor(6);
	pp511[position]->Scale(fParamVal[1]);
	pp511[position]->SetLineColor(2);
	pp511[position]->SetMarkerColor(2);
	pp1275[position]->Scale(fParamVal[2]);
	pp1275[position]->SetLineColor(3);
	pp1275[position]->SetMarkerColor(3);
	compton511[position]->Scale(fParamVal[3]);
	compton511[position]->SetLineColor(4);
	compton511[position]->SetMarkerColor(4);
	compton1275[position]->Scale(fParamVal[4]);
	compton1275[position]->SetLineColor(5);
	compton1275[position]->SetMarkerColor(5);
		

	thisstack->Add(cur_bg);
	thisstack->Add(cur_cbg);
	thisstack->Add(pp511[position]);
	thisstack->Add(pp1275[position]);
	thisstack->Add(compton511[position]);
	thisstack->Add(compton1275[position]);
	weights=fParamVal;
	return thisstack;
}

//------------------------------------------------------------------
void TFit::FitSpectra(){
	
	THStack* test =FitSingleSpectrum(0,templateweights[0]);
	test->Print();
	fittedtemplates.push_back(test);
	
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

vector<TH2D*> TFit::GetChi2Map(){
	return Chi2Map;
}

//------------------------------------------------------------------

TH1D* TFit::GetBackground(){
	if(seriesNo<1 && seriesNo>5) cout << "The background is not properly defined" << endl;
	return background;
}

//------------------------------------------------------------------

Double_t Templatefit_function(int bin,Double_t *par)
{
 double value= par[0]*cur_bg->GetBinContent(bin)+par[1]*cur_pp511->GetBinContent(bin)+par[2]*cur_pp1275->GetBinContent(bin)+par[3]*cur_c511->GetBinContent(bin)+par[4]*cur_c1275->GetBinContent(bin);//+par[5]*cur_cbg->GetBinContent(bin);
 return value;
}

void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  //calculate chisquare
  double chisq = 0;
  for (int i=0;i<Nbins; i++) {
    // chi square is the quadratic sum of the distance from the point to the function weighted by its error
    double delta  = 0;
    if(cur_spec->GetBinContent(i)>0) delta= (cur_spec->GetBinContent(i)-Templatefit_function(i,par))/cur_spec->GetBinError(i);
    chisq += delta*delta;
  }
  f = chisq;
  return;
}
