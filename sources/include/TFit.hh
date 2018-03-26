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
	
	SFData* data; ///< Data of the measurement series
	int Nspectra; ///< Number of spectra in the measurement series
	
	SFData* bgdata; ///< Data of the background measurement series
	
	SFMC* template_data; ///< Data of the simulated templates
	
	TRandom3* resolutiongenerator; 
	TRandom3* anglegenerator;
	TRandom3* comptongenerator;
	
	// Spectra that will analysed
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
	
	std::vector<TGraphErrors*> residual; ///< residuals of the fitted spectra
	
	
	std::vector<THStack*> fittedtemplates; ///<THStacks that contain the fitted, weighted templates 
	std::vector<double*> templateweights; ///< Prozentual weights of the different templates of the spectra
	
	std::vector<SFPeakFinder*> Peaks; ///< Peaks of the 511 keV Peak, needed for the start parameter of the energy model 

	std::vector<TH2D*> Chi2Map;
	
	// Fitting function called in Constructer
	THStack* FitSingleSpectrumSplit(int position, double *&weights); 
	THStack* FitSingleSpectrumTwoSplit(int position, double *&weights);
	THStack* FitSSSumEnergy(int position, double* &weight);
	THStack* FitSSSumWeights(int position, double* &weight);
	void FitSpectra();
	
	void Reset();

	
public:
	TFit();
	TFit(int series_No);
	~TFit();

	bool SetDetails(int serie_No);

	std::vector<TH1D*> GetSpectra();
	TH1D* GetBackground();
	std::vector <THStack*> GetFittedTemplates();
	std::vector <TH2D*> GetChi2Map();
	std::vector <TGraphErrors*> GetResiduals();
  
	ClassDef(TFit,1)
  
};

#endif
