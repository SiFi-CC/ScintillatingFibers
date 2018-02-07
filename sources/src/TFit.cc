// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               TFit.cc                 *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "TFit.hh"
using namespace std;

ClassImp(TFit);

//------------------------------------------------------------------
TFit::TFit():name(""),position(0),spectrum(NULL),background(NULL)
{
 cout << "##### Warning in TFit constructor!" << endl;
 cout << "You are using the default constructor. No Fit possible without the respective spectra." <<endl;
}
TFit::~TFit(){
	
}
