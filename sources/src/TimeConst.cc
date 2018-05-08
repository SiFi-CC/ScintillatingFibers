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
TimeConst::TimeConst():tmax(0),tsplit(0),option(0)
{
	cout << "##### Warning in TimeConst constructor!" << endl;
	cout << "You are using the default constructor. Determination of the timeconstants not possible without the respective spectra." <<endl;
}

//------------------------------------------------------------------
///Standard constructor (recommended). 
///\param series_No Series number of the measurement series that is analyzed. For the numbering see SFData. 
///\param Option (string) defines mode of analysis ( single or one assumes single decay mode, double or two assumes two decay modes

TimeConst::TimeConst(int series_No, string Option):tmax(0),tsplit(0)
{
	SetDetails(series_No,Option);
}

//------------------------------------------------------------------
///Default deconstructor. 
TimeConst::~TimeConst(){
}

//------------------------------------------------------------------
///Sets all details of the signal that is analyzed, all parameters are set and the signal is fitted.
/// The following attributes are st within this function 
///\param series_No Series number of the measurement series that is analyzed. For the numbering see SFData. 
///\param Option (string) defines mode of analysis ( single or one assumes single decay mode, double or two assumes two decay modes


bool TimeConst::SetDetails(int series_No, string Option){
	tmax=0;
	tsplit=0;
	NperPoint=3;
	
	///- fitting option
	if(Option.find("one")!=string::npos || Option.find("single")!=string::npos) option=1;
	if(Option.find("two")!=string::npos || Option.find("double")!=string::npos) option=2;
	
	
	///- the signal
	seriesNo=series_No;
	if (seriesNo>0 && seriesNo<10) data = new SFData(seriesNo);
	else return false;
	
	Npoints=data->GetNpoints();
	signals.reserve(Npoints*NperPoint);
	
	for(int i=0;i<Npoints*NperPoint;i++){
		signals[i]=data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>59.5 && ch_0.fPE<60.5",100,true);
		signals[Npoints+i]= data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>119.5 && ch_0.fPE<120.5",100,true);
		signals[Npoints*2+i]= data->GetSignalAverage(0,(i+1)*10,"ch_0.fPE>199.5 && ch_0.fPE<200.5",100,true);
	}
	
	GraphicSlowComp= new TGraphErrors(Npoints*NperPoint);
	GraphicFastComp= new TGraphErrors(Npoints*NperPoint);
	

	
}

Double_t* TimeConst::FitSingleSignal(TProfile* Signal){
	Double_t* fitresults;
	double lastchi=0;
	///- tmax	
	tmax = Signal->GetMaximumBin();
	
	if(option==1){
		for(int i=0;i<5;i++){
			singleexp= new TF1("singleexp","expo",tmax+(i+1)*10,1024);
			Signal->Fit("singleexp","R");
			if(lastchi==0 || lastchi> singleexp->GetChisquare()){ 
				lastchi=singleexp->GetChisquare();
				fitresults = new Double_t[4];
				fitresults[0]= singleexp->Eval(tmax);
				fitresults[1]= 1/singleexp->GetParameter(1);
				fitresults[2]= singleexp->Eval(tmax);
				fitresults[3]= singleexp->GetParameter(1)/singleexp->GetParError(1)/singleexp->GetParError(1);
			}
		}
	}
	
	else if(option==2){
		for(int i=0;i<5;i++){
			for(int j=0;j<5;j++){
				for( int k = tmax; k<Signal->GetNbinsX(); k++){
					if (Signal->GetBinContent(k) < (Signal->GetMaximum()/(j+2))) {
						tsplit = k;
						break;
					}
				}
				doubleexp= new TF1("doubleexp","expo(0)+expo(2)",tmax+(i+1)*10,1024);
				fastexp= new TF1("fastexp","expo",tmax+(i+1)*10,tsplit);
				slowexp= new TF1("slowexp","expo",tsplit,1024);
				Signal->Fit("fastexp","R");
				Signal->Fit("slowexp","R");
				doubleexp->SetParameters(fastexp->GetParameter(0),fastexp->GetParameter(1),slowexp->GetParameter(0),slowexp->GetParameter(1));
				Signal->Fit("doubleexp","R");
				if(lastchi==0 || lastchi> doubleexp->GetChisquare()){
					fastexp->SetParameters(doubleexp->GetParameter(0),doubleexp->GetParameter(1));
					slowexp->SetParameters(doubleexp->GetParameter(2),doubleexp->GetParameter(3));
					lastchi=doubleexp->GetChisquare();
					fitresults = new Double_t[8];
					fitresults[0]= fastexp->Eval(tmax);
					fitresults[2]= 1/doubleexp->GetParameter(1);
					fitresults[4]= slowexp->Eval(tmax);
					fitresults[6]= 1/doubleexp->GetParameter(3);
					fitresults[1]= fastexp->Eval(tmax);
					fitresults[3]= TMath::Abs(doubleexp->GetParError(1)/doubleexp->GetParameter(1)/doubleexp->GetParameter(1));
					fitresults[5]= slowexp->Eval(tmax);
					fitresults[7]= TMath::Abs(doubleexp->GetParError(3)/doubleexp->GetParameter(3)/doubleexp->GetParameter(3));
				}
			}
		}
	}
	return fitresults;
}
void TimeConst::FitSignals(){
	for(int i=0;i<signals.size();i++){
		fitresults.push_back(FitSingleSignal(signals[i]));
		GraphicFastComp->SetPoint(i,i,fitresults[i][2]);
		GraphicFastComp->SetPointError(i,0,fitresults[i][3]);
		GraphicSlowComp->SetPoint(i,i,fitresults[i][6]);
		GraphicSlowComp->SetPointError(i,0,fitresults[i][7]);
	}
	
}


TGraphErrors* TimeConst::GetSlowComponent(){
	return GraphicSlowComp;		
}

TGraphErrors* TimeConst::GetFastComponent(){
	return GraphicFastComp;		
}
///Print method. 
///Prints the result of the fitting process according to the used methode 

//~ void TimeConst::Print(){
	//~ if(option==1){
		//~ cout << "Assuming a single decay constant it is determined to " << fitresults[1] << "ns." << endl;
		//~ cout << "The amplitude is determined to " << fitresults[0] << "mV." << endl; 
	//~ }
	//~ else if (option==2){
		//~ cout << "Assuming a decay  build by two decay constants" << endl; 
		//~ cout << "The first one is determined to " << fitresults[1] << "ns." << endl;
		//~ cout << "The  corresponding amplitude is determined to " << fitresults[0] << "mV." << endl; 
		//~ cout << "The second one is determined to " << fitresults[3] << "ns." << endl;
		//~ cout << "The  corresponding amplitude is determined to " << fitresults[2] << "mV." << endl; 
	//~ }
	//~ else {
		//~ cout << "The option was not right defined so no constants are determined." << endl; 	
	//~ }
//~ }
