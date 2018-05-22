// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              TFit.hh                  *
// *            Jonas Kasper               *
// *    kasper@physik.rwth-aachen.de       *        
// *          Created in 2018              *
// *                                       *
// *****************************************


#ifndef __TFit_H_
#define __TFit_H_ 1
#include "SFData.hh"
#include "SFMC.hh"
#include "SFPeakFinder.hh"
#include "TObject.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TGraphErrors.h" 
#include "TF1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <stdlib.h>
#include <vector>


static TH1D* cur_spec;

static TH1D* cur_bg;
static TH1D* cur_ebg;
static TH1D* cur_pp511;
static TH1D* cur_pp1275;
static TH1D* cur_c511;
static TH1D* cur_c1275;
static TH1D* cur_pall;
static TH1D* cur_call;
static TH1D* cur_sum;
static TRandom3* resgenerator;
static int Nbins;

static int Fitboarder_low;
static int Fitboarder_up;
static int minentry;
static int steps_par0;
static int steps_par1;
static int steps_par2;


static void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Double_t Templatefit_functionsplit(int bin,Double_t* par);

static void calc_chi_square_twosplit(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Double_t Templatefit_functiontwosplit(int bin,Double_t* par);

static void calc_chi_square_sum(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Double_t Templatefit_functionsumenergy(Double_t sum_cont,Double_t bg_cont,Double_t* par);

static void calc_chi_square_weights(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Double_t Templatefit_functionsumweights(Double_t sum_cont,Double_t bg_cont,Double_t* par);

static Double_t Energy_model(Double_t energy,Double_t* par);

class TFit : public TObject{
  
private:
  
	Int_t seriesNo; ///< Series number that is analysed by the Fit
	Int_t fmode;///< Inticates the fit that is performed: 0 Total Split, 1 Split in Compton and Photopeak, 2 No split
	Int_t rebin;

	
	SFData* data; ///< Data of the measurement series
	int Nspectra; ///< Number of spectra in the measurement series
	
	SFData* bgdata; ///< Data of the background measurement series
	
	SFMC* template_data; ///< Data of the simulated templates
	
	TRandom3* resolutiongenerator; ///<Used to smear the templates 
	
	std::vector<TH1D*> spectra; ///< Spectra that will be analysed 
	
	// Needed Templates
	TH1D* background; ///< Internal activity background template
	
	std::vector<TH1D*> templates_else; ///< raw templates of other depositions within the fibres
	std::vector<TH1D*> templates_c511; ///< raw templates of the 511 keV  compton spectra
	std::vector<TH1D*> templates_c1275; ///< raw templates of the 1275 keV compton spectra
	std::vector<TH1D*> templates_p511; ///< raw templates of the 511 keV photon peak
	std::vector<TH1D*> templates_p1275; ///< raw templates of the 1275 keV compton spectra
	std::vector<TH1D*> templates_pall; ///< raw templates of the sumed photon peak
	std::vector<TH1D*> templates_call; ///< raw templates of the sumed compton spectra
	std::vector<TH1D*> templates_sum; ///< raw templates of the sumed depositions within the fibre
	
	std::vector<TH1D*> pp511; ///< scaled 511 keV photonpeak templates
	std::vector<TH1D*> compton511; ///< scaled 1275 keV photonpeak templates
	std::vector<TH1D*> pp1275; ///< scaled 511 keV compton spectra templates
	std::vector<TH1D*> compton1275; ///< scaled 1275 keV compton spectra templates
	std::vector<TH1D*> pall; ///< scaled sumed photonpeak templates
	std::vector<TH1D*> call; ///< scaled sumed compton spectra templates
	std::vector<TH1D*> ebg; ///< background of other depositions within the fibre
	
	
	
	std::vector<THStack*> fittedtemplates; ///<THStacks that contain the fitted, weighted templates 
	std::vector<TGraphErrors*> residual; ///< residuals of the fitted spectra
	TH2D* GraphicWeights;///<Contains all Weights of each position
	
	TGraph* ec_a; ///<Graph containing the Constant [0] of the energy model for the best fits of the positions
	TGraph* ec_b;///<Graph containing the Constant [1] of the energy model for the best fits of the positions
	TGraph* ec_c;///<Graph containing the Constant [2] of the energy model for the best fits of the positions
	TGraph* chigraph;///<Graph containing the chi/ndf of the best fits 
	TGraphErrors* g_pp_511;///<Graph containing the weights of the pp_511 templates if this fit is executed 
	TGraphErrors* g_c_511;///<Graph containing the weights of the c_511 templates if this fit is executed  
	TGraphErrors* g_pp_1275;///<Graph containing the  weights of the pp_1275 templates if this fit is executed 
	TGraphErrors* g_c_1275;///<Graph containing the weights of the cc_1275 templates if this fit is executed  
	TGraphErrors* g_ibg;///<Graph containing the weights of the ibg templates if this fit is executed  
	TGraphErrors* g_ebg;///<Graph containing the weights of the ebg templates if this fit is executed 
	TGraphErrors* g_sum;///<Graph containing the weights of the geant4_sum templates if this fit is executed 
	
	std::vector<double*> templateweights; ///< Prozentual weights of the different templates of the spectra
	std::vector<double*> energypara; ///< energy calibration parameter
	std::vector<double> chindf; ///< Chi2 der fits
	
	std::vector<SFPeakFinder*> Peaks; ///< Peaks of the 511 keV Peak, needed for the start parameter of the energy model 

	// Fitting function called in Constructer
	//------------------------------------------------------------------
	///Function that performes a fit of all templates to the spectrum that is defined by the position
	///\param position - defining spectrum that the templates are fitted to
	///\param *weights - the fit weights are returned here
	///\param *enerconst - the parameter of the energy model are returned here
	THStack* FitSingleSpectrumSplit(int position, double *weights, double* enerconst); 
	///Function that performes a fit of a compton and photon template to the spectrum that is defined by the position
	///\param position - defining spectrum that the templates are fitted to
	///\param *weights - the fit weights are returned here
	///\param *enerconst - the parameter of the energy model are returned here
	THStack* FitSingleSpectrumTwoSplit(int position, double *weights, double* enerconst);
	///Function that performes a fit of sum template and the energy constants to the spectrum that is defined by the position
	///\param position - defining spectrum that the templates are fitted to
	///\param *weight - the fit weights are returned here
	///\param *enerconst - the parameter of the energy model are returned here
	THStack* FitSSSumEnergy(int position, double* weight, double* enerconst);
	///Function that performes a fit of a sum of all templates to the spectrum that is defined by the position
	///\param position - defining spectrum that the templates are fitted to
	///\param *weight - the fit weights are returned here
	///\param *enerconst - the parameter of the energy model are returned here
	THStack* FitSSSumWeights(int position, double* weight, double* enerconst);
	
	///Function to set up the TMinuit for the different fitting modes 
	///\param *fParam - returns fitting parameter
	///\param *fParErr - returns errors on fitting parameter
	///\param &chi - retruns the chi/ndf 
	///\param fitmode - sets the fitmode for the TMinuit
	int SetTMiniut(double *fParam,double *fParErr,double &chi,int fitmode);
	
	///Function to call the Fit funcitons for each Position of a series
	void FitSpectra();
	
	/// Function to reset everything
	void Reset();

	
public:
	///Default constructor.
	///If this constructer is used the details need to be called with SetDetails().
	TFit();
	///Constructor that takes the series number of the measurement series that is analysed.
	///\param series_No Series number of the measurement series that is analysed. For numbering se SFData class.
	///\param FitMode - specifies the Fit type that is done 
	///\param re_bin - rebining factor for the spectra
	TFit(int series_No,TString FitMode, int re_bin);
	///Deconstructer
	~TFit();
	///Function that sets up the series and the templates specifications that are requested
	///\param series_No - Series to be analyzed
	///\param FitMode - Defining the fit that is performed on the spectra of the serie
	bool SetDetails(int serie_No,TString FitMode, int re_bin);

	std::vector<TH1D*> GetSpectra();///Returns the spectra of the series
	TH1D* GetBackground();///Returns the internal background for the series
	std::vector <THStack*> GetFittedTemplates();///returns the THStacks containing the fitted templates
	TH2D* GetGraphicWeights(); ///Returns visualisation of all weights of each position 
	TGraph* GetEC_a();///Returns Graph with the [0] parameter of the energy model for each position
	TGraph* GetEC_b();///Returns Graph with the [1] parameter of the energy model for each position
	TGraph* GetEC_c();///Returns Graph with the [2] parameter of the energy model for each position
	TGraph* GetChi();///Returns Graph with the Chi/ndf of the fits at each position
	TGraphErrors* GetIBG_W();///Returns TGraphErrors with the weights of the IBG for each position
	TGraphErrors* GetPP511_W();///Returns TGraphErrors with the weights of the PP511 for each position
	TGraphErrors* GetPP1275_W();///Returns TGraphErrors with the weights of the PP1275 for each position
	TGraphErrors* GetEBG_W();///Returns TGraphErrors with the weights of the EBG for each position
	TGraphErrors* GetC511_W();///Returns TGraphErrors with the weights of the C511 for each position
	TGraphErrors* GetC1275_W();///Returns TGraphErrors with the weights of the C1275 for each position
	TGraphErrors* GetSum_W();///Returns TGraphErrors with the weights of the sum_geant4 for each position
	std::vector <TGraphErrors*> GetResiduals();///Returns TGraphErrors with the resiudlas to the spectrum of each position
  
	ClassDef(TFit,1)
  
};

#endif
