// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFData.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
using namespace std;

ClassImp(SFData);

//------------------------------------------------------------------
SFData::SFData(){
 fSeriesNo = 0; 
}
//------------------------------------------------------------------
SFData::SFData(Int_t seriesNo){
  fSeriesNo = seriesNo;
}
//------------------------------------------------------------------
SFData::~SFData(){
}
//------------------------------------------------------------------
void SFData::Print(void){
 cout << "\n\n------------------------------------------------" << endl;
 cout << "This is Print() for SFData class object" << endl;
 cout << "Number of the experimental series: " << fSeriesNo << endl;
}
//------------------------------------------------------------------