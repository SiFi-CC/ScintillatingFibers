// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFData.hh                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFData_H_
#define __SFData_H_ 1
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

class SFData : public TObject{
  
private:
  Int_t fSeriesNo;
  
public:
  SFData();
  SFData(Int_t seriesNo);
  ~SFData();
  
  void Print(void);
  
  ClassDef(SFData,1)
  
};

#endif