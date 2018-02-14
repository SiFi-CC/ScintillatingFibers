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
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
using namespace std;

ClassImp(TFit);

//------------------------------------------------------------------
///Default constructor.
///If this constructer is used the details need to be called with SetDetails().
TFit::TFit():name(""),position(0),spectrum(NULL),background(NULL)
{
 cout << "##### Warning in TFit constructor!" << endl;
 cout << "You are using the default constructor. No Fit possible without the respective spectra." <<endl;
}

//------------------------------------------------------------------
TFit::TFit(TH1D* Spectrum, TH1D* Background)
{
	SetDetails(Spectrum,Background,"DefaultName",0);

	pp511= PhotoPeak(Calipar[0],Calipar[1]);
	pp1275= PhotoPeak(Calipar[0]/511*1275,Calipar[1]);
	compton511 = Compton(Calipar[0],Calipar[1]);
	compton1275= Compton(Calipar[0]/511*1275,Calipar[1]);
	Fitresult= Fitting();
}

//------------------------------------------------------------------
TFit::TFit(TH1D* Spectrum, TH1D* Background,TString Name,Double_t Position)
{
 SetDetails(Spectrum,Background,Name,Position);
 
 pp511 = PhotoPeak(511,50);
 pp1275= PhotoPeak(1275,50);
 compton511 = Compton(511,50);
 compton1275= Compton(1275,50);
}

//------------------------------------------------------------------
TFit::~TFit(){
	
}

//------------------------------------------------------------------
void TFit::SetDetails(TH1D* Spectrum, TH1D* Background,TString Name,Double_t Position){
	
	spectrum = Spectrum;
	background= Background;
	name = Name;
	position = Position;
	
	anglegenerator = new TRandom3();
	Nbins= spectrum->GetNbinsX();
	resolutiongenerator = new TRandom3();
	comptongenerator =  new TRandom3();
	Calipar= Calibration();
	
}
//------------------------------------------------------------------

Double_t TFit::KNFormular(Double_t PhotoPeakEnergy, Double_t Angle){
	
	Double_t a= 1/(1+((PhotoPeakEnergy/511)*(1-TMath::Cos(Angle*TMath::DegToRad()))));

	Double_t kn = a*a*(a+(1/a)-TMath::Sin(Angle*TMath::DegToRad())*TMath::Sin(Angle*TMath::DegToRad()))/2/511/511/137.036/137.036;
	
	Double_t kn_0 = 1./511./511./137.036/137.036;
	
	return (kn/kn_0);
	
}
//------------------------------------------------------------------

TH1D* TFit::Compton(Double_t PhotoPeakEnergy, Double_t Resolution){
	cout << "Creating the compton spectrum histogram according to an energy of " << PhotoPeakEnergy << " and a resolution of " << Resolution << endl;
	TString histname= Form("cs_%f_%f",PhotoPeakEnergy,Resolution);
	TH1D* cs = new TH1D(histname,histname,Nbins,start,end);
	

	double angle;
	double prob;
	double energy_e;
	for( int i=0; i<50000;i++){
		angle= anglegenerator->Uniform(0,180);
		prob = KNFormular(PhotoPeakEnergy,angle);
		if(comptongenerator->Uniform(0,1)<prob){
			energy_e= PhotoPeakEnergy*(1-1/(1+((PhotoPeakEnergy/511)*(1-TMath::Cos(angle*TMath::DegToRad())))));
			cs->Fill(resolutiongenerator->Gaus(energy_e, Resolution));
		}
	}
	
	return cs;
}
//------------------------------------------------------------------

TH1D* TFit::PhotoPeak(Double_t PhotoPeakEnergy, Double_t Resolution){
	cout << "Creating the photopeak histogram with an energy of " << PhotoPeakEnergy << " and a resolution of " << Resolution << endl;
	TString histname= Form("%f_%f",PhotoPeakEnergy,Resolution);
	TH1D* pp = new TH1D(histname,histname,Nbins,start,end);
	
	for( int i=0; i<50000;i++){
		double temp = resolutiongenerator->Gaus(PhotoPeakEnergy,Resolution);
		pp->Fill(temp);
	}
	
	return pp;
}

//------------------------------------------------------------------

vector <Double_t> TFit::Calibration(){
	TH1D* spec = (TH1D*)spectrum->Clone("Calspec");
	TH1D* bgs= new TH1D("bgs","bgs",spec->GetNbinsX(),spec->GetBinCenter(1),spec->GetBinCenter(spectrum->GetNbinsX()));
  
	BGFit * mybgfit = new BGFit(215,325);
	TF1 *bgfun = new TF1("bgfun",mybgfit,&BGFit::Evaluate,175,450,4,"BGFit","Evaluate");
	
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
vector<Double_t> TFit::Fitting(){
	
	TObjArray *mc = new TObjArray(5);        // MC histograms are put in this array
	mc->Add(background);
	mc->Add(pp511);
	mc->Add(compton511);
	mc->Add(pp1275);
	mc->Add(compton1275);
	TFractionFitter* finalfit = new TFractionFitter(spectrum, mc); // initialise
	finalfit->Constrain(0,0.1,5.0);               // constrain fraction 1 to be between 0 and 1
	finalfit->Constrain(1,0.1,5.0);               // constrain fraction 1 to be between 0 and 1
	finalfit->Constrain(2,0.1,10.0);               // constrain fraction 1 to be between 0 and 1
	finalfit->Constrain(3,0.01,10.0);               // constrain fraction 1 to be between 0 and 1
	finalfit->Constrain(4,0.01,10.0);               // constrain fraction 1 to be between 0 and 1
	finalfit->SetRangeX(250,800);  
	
	Int_t status = finalfit->Fit();               // perform the fit
	std::cout << "fit status: " << status << std::endl;	
	vector <Double_t> temp(5);
	vector <Double_t> temperr(5);
	if(status==0){
		for(int i=0;i<5;i++){
			finalfit->GetResult(i,temp[i],temperr[i]);
		}
	}
	else{
		for(int i=0;i<5;i++){
			temp[i]=1;
		}
	}
	return temp;
}

//------------------------------------------------------------------

TH1D* TFit::GetSpectrum(){
	if(spectrum!=NULL) return spectrum;
	else{
		cout << "The spectrum is not properly defined" << endl;
		return NULL;
	}
}

//------------------------------------------------------------------

TH1D* TFit::GetFittedBackground(){
	if(background!=NULL){
		background->Scale(Fitresult[0]);
		return background;
	}
	else{
		cout << "The background is not properly defined" << endl;
		return NULL;
	}
}

//------------------------------------------------------------------

TH1D* TFit::GetFittedPhotoPeak(Double_t PhotoPeakEnergy){
	if(PhotoPeakEnergy == 511.){
		if(pp511!=NULL){
			pp511->Scale(Fitresult[1]);
			return pp511;
		}
		else{
			cout << "The photopeak at 511 keV is not properly defined" << endl;
			return NULL;
		}
	}
	if(PhotoPeakEnergy == 1275.){
		if(pp1275!=NULL){
			pp1275->Scale(Fitresult[3]);
			return pp1275;
		}
		else{
			cout << "The photopeak at 1275 keV is not properly defined" << endl;
			return NULL;
		}
	}
}

//------------------------------------------------------------------

TH1D* TFit::GetFittedCompton(Double_t PhotoPeakEnergy){
	if(PhotoPeakEnergy == 511.){
		if(compton511!=NULL){
			compton511->Scale(Fitresult[2]);
			return compton511;
		}
		else{
			cout << "The compton spectrum originating from 511 keV is not properly defined" << endl;
			return NULL;
		}
	}
	if(PhotoPeakEnergy == 1275.){
		if(compton1275!=NULL){
			compton1275->Scale(Fitresult[4]);
			return compton1275;
		}
		else{
			cout << "The compton spectrum originating from 1275 keV is not properly defined" << endl;
			return NULL;
		}
	}
}


//------------------------------------------------------------------

vector <TH1D*> TFit::GetFittedSpectra(){
	bool error = false;
	vector <TH1D*> temp(6);
	temp[0]=GetSpectrum();
	temp[1]=GetFittedBackground();
	temp[2]=GetFittedPhotoPeak(511);
	temp[3]=GetFittedCompton(511);
	temp[4]=GetFittedPhotoPeak(1275);
	temp[5]=GetFittedCompton(1275);
	
	for(int i=0; i< temp.size();i++){
		if (temp[i]==NULL) error=true;
	}
	
	if(!error) 	return temp;
	else{
		cout << "ATTENTION FROM THE FIT" << endl;
		return temp;
	} 
}
 
