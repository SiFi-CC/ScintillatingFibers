// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFAttenuation.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFAttenuation_H_
#define __SFAttenuation_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include <iostream>

using namespace std;

class SFAttenuation : public TObject{
 
private:
  int    fSeriesNo;		///< Number of experimental series to be analyzed
  SFData *fData;		///< SFData object of the analyzed series
  
public:
  SFAttenuation();
  SFAttenuation(int seriesNo);
  ~SFAttenuation();
  
  bool AttAveragedCh(void);
  bool AttSeparateCh(int ch);
  void Print(void);
  
  ClassDef(SFAttenuation,1)
};

#endif