// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             TimeConst.cc              *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "TimeConst.hh"

using namespace std;

ClassImp(TimeConst);

//------------------------------------------------------------------
///Default constructor. If this constructor is used the signal and the fitting option needs to be set in set details 
///and the fitting option needs to be set in SetDetails(TProfile signal*, string Option) 
TimeConst::TimeConst():tmax(0),tsplit(0),signal(NULL),option(0)
{
	cout << "##### Warning in TimeConst constructor!" << endl;
	cout << "You are using the default constructor. Determination of the timeconstants not possible without the respective spectra." <<endl;
}

//------------------------------------------------------------------
///Standard constructor (recommended). 
///\param Signal (TProfile*) is the average signal that shall be analyzed 
///\param Option (string) defines mode of analysis ( single or one assumes single decay mode, double or two assumes two decay modes
TimeConst::TimeConst(TProfile* Signal, std::string Option):tmax(0),tsplit(0)
{
	SetDetails(Signal,Option);
}

//------------------------------------------------------------------
///Default deconstructor. 
TimeConst::~TimeConst(){
	if(fitresults!=NULL) delete fitresults;
	if(signal!= NULL) delete signal;
}

//------------------------------------------------------------------
///Sets all details of the signal that is analyzed, all parameters are set and the signal is fitted.
/// The following attributes are st within this function 
///\param Signal (TProfile*) is the average signal that shall be analyzed 
///\param Option (string) defines mode of analysis ( single or one assumes single decay mode, double or two assumes two decay modes)


void TimeConst::SetDetails(TProfile* Signal, string Option){
	///- the signal
	signal=Signal;
	
	///- fitting option
	if(Option.find("one")!=string::npos || Option.find("single")!=string::npos) option=1;
	if(Option.find("two")!=string::npos || Option.find("double")!=string::npos) option=2;
	
	///- tmax	
	tmax = signal->GetMaximumBin();
	
	///- tsplit
	if(option==2){
		for( int i = tmax; i<signal->GetNbinsX(); i++){
			if (signal->GetBinContent(i) < (signal->GetMaximum()/10)) {
				tsplit = i;
				break;
			}
		}
	}
	
	cout << "The Maximum of the signal is placed at " << tmax << "ns." << endl;
	if(option==2)cout << "The signal is split at the position " << tsplit << "ns." << endl;
	///-  The fit functions
	if(option==1)singleexp= new TF1("singleexp","expo",tmax+10,1024);
	else if(option==2){
		doubleexp= new TF1("doubleexp","expo(0)+expo(2)",tmax+10,1024);
		fastexp= new TF1("fastexp","expo",tmax+10,tsplit);
		slowexp= new TF1("slowexp","expo",tsplit,1024);
	}
	Fitting();
}

void TimeConst::Fitting(){
	if(option==1){
		signal->Fit("singleexp","R");
		fitresults = new Double_t[2];
		fitresults[0]= singleexp->Eval(tmax);
		fitresults[1]= 1/singleexp->GetParameter(1);
	}
	else if(option==2){
		signal->Fit("fastexp","R");
		signal->Fit("slowexp","R");
		doubleexp->SetParameters(fastexp->GetParameter(0),fastexp->GetParameter(1),slowexp->GetParameter(0),slowexp->GetParameter(1));
		signal->Fit("doubleexp","R");
		fastexp->SetParameters(doubleexp->GetParameter(0),doubleexp->GetParameter(1));
		slowexp->SetParameters(doubleexp->GetParameter(2),doubleexp->GetParameter(3));
		fitresults = new Double_t[4];
		fitresults[0]= fastexp->Eval(tmax);
		fitresults[1]= 1/doubleexp->GetParameter(1);
		fitresults[2]= slowexp->Eval(tmax);
		fitresults[3]= 1/doubleexp->GetParameter(3);
	}  
}
///Print method. 
///Prints the result of the fitting process according to the used methode 

void TimeConst::Print(){
	if(option==1){
		cout << "Assuming a single decay constant it is determined to " << fitresults[1] << "ns." << endl;
		cout << "The amplitude is determined to " << fitresults[0] << "mV." << endl; 
	}
	else if (option==2){
		cout << "Assuming a decay  build by two decay constants" << endl; 
		cout << "The first one is determined to " << fitresults[1] << "ns." << endl;
		cout << "The  corresponding amplitude is determined to " << fitresults[0] << "mV." << endl; 
		cout << "The second one is determined to " << fitresults[3] << "ns." << endl;
		cout << "The  corresponding amplitude is determined to " << fitresults[2] << "mV." << endl; 
	}
	else {
		cout << "The option was not right defined so no constants are determined." << endl; 	
	}
}
///Methode that returns the results of the parameters 
///\return Double_t*[2] or Double_t*[4] 

Double_t* TimeConst::GetFitData(){
	return fitresults;	
}
///Methode that draws the fitted signal including the fits 
///\param name (string) of the canvas
///\return TCanvas* 

TCanvas* TimeConst::DrawFittedSignal(string name){
	TCanvas* can = new TCanvas(name.c_str(),name.c_str(),1000,800);
	signal->Draw("");
	if (option==1) singleexp->Draw("same");	
	else if (option==2)  doubleexp->Draw("same");
	
	return can;	
}
