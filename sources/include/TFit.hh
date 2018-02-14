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
#include "TObject.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <stdlib.h>
#include <vector>

class TFit : public TObject{
  
private:
  
	TString name; ///< Name for the fit, for identification purpose
	Double_t position; ///< position of the spectrum, for identification purpose
	
	Int_t Nbins; ///< Number of bins of the spectra
	Double_t start; ///< Position of the first bin of the spectra
	Double_t end; ///< Position of the last bin of the spectra
	std::vector<Double_t> Calipar;
	std::vector<Double_t> Fitresult;
	
	TRandom3* resolutiongenerator;
	TRandom3* anglegenerator;
	TRandom3* comptongenerator;
	
	// Spectrum that is analysed
	TH1D* spectrum; ///< Spectrum that will be analysed 
	
	// Needed Templates
	TH1D* background; ///< Internal activity background template
	TH1D* pp511; ///< Photopeak of 511 keV template
	TH1D* compton511; ///< Compton spectrum of 511 keV template
	TH1D* pp1275; ///< Photopeak of 1275 keV template
	TH1D* compton1275; ///< Compton spectrom of 1275 keV template
	
	// Fitting function called in Constructer
	std::vector<Double_t> Fitting();
	TH1D* Compton(Double_t PhotoPeakEnergy, Double_t Resolution);
	TH1D* PhotoPeak(Double_t PhotoPeakEnergy, Double_t Resolution);
	Double_t KNFormular(Double_t PhotoPeakEnergy, Double_t Angle);
	std::vector<Double_t> Calibration();
	
public:
	TFit();
	TFit(TH1D* Spectrum,TH1D* Background);
	TFit(TH1D* Spectrum,TH1D* Background, TString Name, Double_t position);
	~TFit();

	void SetDetails(TH1D* SPectrum, TH1D* Background, TString Name, Double_t Position);
	
	TH1D* GetFittedPhotoPeak(Double_t PhotoPeakEnergy);
	TH1D* GetFittedCompton(Double_t PhotoPeakEnergy);
	TH1D* GetSpectrum();
	TH1D* GetFittedBackground();
	std::vector <TH1D*> GetFittedSpectra();
	Double_t* GetFitData();
  
	ClassDef(TFit,1)
  
};

#endif
