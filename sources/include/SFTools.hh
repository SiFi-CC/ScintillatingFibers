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
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include <sqlite3.h>

class SFTools : public TObject{
    
public:
    /// Standard constructor.
    SFTools()  {};
    /// Standard destructor.
    ~SFTools() {};
    
    static int     GetIndex(std::vector <double> positions, double position);
    static int     GetSeriesNo(TString hname_tstr);
    static int     GetChannel(TString hname_tstr);
    static double  GetPosition(TString hname_tstr);
    static double  GetPosError(TString collimator, TString testBench);
    static void    CheckDBStatus(int status, sqlite3 *database); 
    static bool    SaveResultsDB(TString database, TString table, 
                                 TString query, int seriesNo);
    static bool    CreateTable(TString database, TString table);  
    static double  GetMean(std::vector <double> vec);
    static double  GetStandardErr(std::vector <double> vec);
    static TString FindData(TString directory);
    
    ClassDef(SFTools,1)
};

#endif
