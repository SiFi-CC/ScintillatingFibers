// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFAttenuation.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// ***************************************** 

#include "SFAttenuation.hh"

ClassImp(SFAttenuation);

//------------------------------------------------------------------
///Default constructor.
SFAttenuation::SFAttenuation(){
  cout << "#### Warning in SFAttenuation constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  fSeriesNo = -1;
  fData = NULL;
}
//------------------------------------------------------------------
///Standard constructor (reccommended)
///\param seriesNo is number of experimental series to be analyzed. 
SFAttenuation::SFAttenuation(int seriesNo){
  fSeriesNo = seriesNo;
  fData = new SFData(fSeriesNo);
}
//------------------------------------------------------------------
///Default destructor.
SFAttenuation::~SFAttenuation(){
  if(fData!=NULL) delete fData;
}
//------------------------------------------------------------------

//------------------------------------------------------------------
///Prints details of SFAttenuation class object.
void SFAttenuation::Print(void){
 cout << "This is print out of SFAttenuation class object" << endl;
 cout << "Experimental series number " << fSeriesNo << endl;
}
//------------------------------------------------------------------