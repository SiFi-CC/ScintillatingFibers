// *****************************************
// *                                       *
// *        ScintillatingFibers            *
// *           TimeConst.hh               *
// *            Jonas Kasper               *
// *    kasper@physik.rwth-aachen.de       *        
// *          Created in 2018              *
// *                                       *
// *****************************************


#ifndef __TimeConst_H_
#define __TimeConst_H_ 1
#include "SFData.hh"
#include "TObject.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"

#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

class TimeConst : public TObject{
  
private:

  Int_t seriesNo; 	///< Series number of the measearument series that is analysed 
  SFData* data; 	///< Data of the measurement series 
  Int_t Npoints; 	///< Number of Points in measurement series
  Int_t NperPoint; 	///< Number of Tprofils that are analysed per Point in the measurement series

  Int_t tmax; 			///< Time of maximal signal amplitude
  Int_t tsplit; 		///< Time, where the slow decay process becomes dominant
  vector <double> fitresults; 	///< Array containing the determined fit data

  vector <TProfile*> signals; 	///< Averraged signals that are analysed
  
  int option; 			///< option of different fitting modes
  
  TF1* singleexp;		///< function for the assumption of a single decay mode
  TF1* doubleexp;		///< function for the assumption of a double decay mode
  TF1* fastexp;			///< function for the assumption of a double decay mode
  TF1* slowexp;			///< function for the assumption of a double decay mode
  
  Double_t* FitSingleSignal(TProfile* Signal);
  void FitSignals();
  
  void Reset();
  
public:
  TimeConst();
  //~ TimeConst(TProfile* Signal, std::string Option, int start, double fraction);
  TimeConst(int series_No, std::string Option);
  ~TimeConst();
  
  bool SetDetails(int series_No, std::string Option );
  std::vector <TProfile*> GetSignals();
  std::vector <Double_t*> GetFitResults();
  std::vector <Double_t*> GetAveragedSignals();
  void Print(); 
    
  ClassDef(TimeConst,1)
  
};

#endif
