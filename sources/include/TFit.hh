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
#include "TObject.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TF1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <stdlib.h>
#include <vector>


static TH1D* cur_spec;

static TH1D* cur_bg;
static TH1D* cur_cbg;
static TH1D* cur_pp511;
static TH1D* cur_pp1275;
static TH1D* cur_c511;
static TH1D* cur_c1275;
static int Nbins;

static void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Double_t Templatefit_function(int bin,Double_t* par);

class TFit : public TObject{
  
private:
  
	Int_t seriesNo; ///< Series number that is analysed by the Fit
	
	SFData* data; ///< Data of the measurement series
	int Nspectra; ///< NUmber of spectra in the measurement series
	int fpt; ///< NUmber of variasions that are fitted (for resolution and position each)
	
	SFData* bgdata; ///< Data of the measurement series
	
	
	
	TRandom3* resolutiongenerator; 
	TRandom3* anglegenerator;
	TRandom3* comptongenerator;
	TRandom3* cbggenerator;
	
	// Spectrum that is analysed
	std::vector<TH1D*> spectra; ///< Spectra that will be analysed 
	
	// Needed Templates
	TH1D* background; ///< Internal activity background template
	
	std::vector<TH1D*> pp511; ///< Photopeaks of 511 keV template
	std::vector<TH1D*> compton511; ///< Compton spectra of 511 keV template
	std::vector<TH1D*> pp1275; ///< Photopeaks of 1275 keV template
	std::vector<TH1D*> compton1275; ///< Compton spectra of 1275 keV template
	
	std::vector<THStack*> fittedtemplates; ///<THStacks that contain the fitted templates 
	std::vector<double*> templateweights;
	
	// Fitting function called in Constructer
	THStack* FitSingleSpectrum(int position, double *&weights);
	void FitSpectra();
	
	void Reset();

	
	TH1D* Compton(TH1D* spec,Double_t PhotoPeakEnergy,Double_t califac, Double_t Resolution);
	TH1D* PhotoPeak(TH1D* spec,Double_t PhotoPeakEnergy, Double_t Resolution);
	TH1D* ColiBg(TH1D* spec,Double_t Resolution);
	Double_t KNFormular(Double_t PhotoPeakEnergy, Double_t Angle);
	std::vector<Double_t> Calibration(TH1D* spectrum);

	std::vector<TH2D*> Chi2Map;
	
public:
	TFit();
	TFit(int series_No, int fitpoints);
	~TFit();

	bool SetDetails(int serie_No, int fitpoints);

	std::vector<TH1D*> GetSpectra();
	TH1D* GetBackground();
	std::vector <THStack*> GetFittedTemplates();
	std::vector <TH2D*> GetChi2Map();

  
	ClassDef(TFit,1)
  
};

#endif
