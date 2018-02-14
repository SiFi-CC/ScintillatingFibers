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
#include "TObject.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"

#include <string>
#include <iostream>
#include <stdlib.h>

class TimeConst : public TObject{
  
private:

  Int_t tmax; ///< Time of maximal signal amplitude
  Int_t tsplit; ///< Time, where the slow decay process becomes dominant
  Double_t* fitresults; ///< Array containing the determined fit data
  Double_t* fitresultserror; ///< Array containing the determined fit data uncertainties 
  
  

  TProfile* signal; ///< Averraged signal that is analysed
  
  int option; ///< option of different fitting modes
  
  TF1* singleexp; ///< function for the assumption of a single decay mode
  TF1* doubleexp;///< function for the assumption of a double decay mode
  TF1* fastexp;///< function for the assumption of a double decay mode
  TF1* slowexp;///< function for the assumption of a double decay mode
  
  void Fitting();
  
public:
  TimeConst();
  TimeConst(TProfile* Signal, std::string Option, int start, double fraction);
  ~TimeConst();
  
  void SetDetails(TProfile* Signal,std::string Option, int start, double fraction);
  TCanvas* DrawFittedSignal(std::string name);
  Double_t* GetFitData();
  Double_t* GetFitDataError();
  void Print(); 
    
  ClassDef(TimeConst,1)
  
};

#endif
