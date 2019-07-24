// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFTools.hh               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFTools_H_
#define __SFTools_H_ 1
#include "TObject.h"
#include "TString.h"
#include <iostream>

class SFTools : public TObject{
    
public:
    /// Standard constructor.
    SFTools()  {};
    /// Standard destructor.
    ~SFTools() {};
    
    static int    GetIndex(std::vector <double> positions, double position);
    static int    GetSeriesNo(TString hname_tstr);
    static double GetPosError(TString collimator, TString testBench);
    
    ClassDef(SFTools,1)
};

#endif
