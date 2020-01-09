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
#include "TFile.h"
#include "TF1.h"
#include "SFData.hh"
#include <iostream>
#include <sqlite3.h>

namespace SFTools{
    
    int     GetIndex(std::vector <int> measurementsIDs, int id);
    int     GetSeriesNo(TString hname_tstr);
    int     GetChannel(TString hname_tstr);
    double  GetPosition(TString hname_tstr);
    int     GetMeasurementID(TString hname_tstr);
    int     GetMeasurementID(int seriesNo, double position);
    double  GetPosError(TString collimator, TString testBench);
    void    CheckDBStatus(int status, sqlite3 *database); 
    bool    SaveResultsDB(TString database, TString table, 
                                 TString query, int seriesNo);
    bool    CreateTable(TString database, TString table);  
    double  GetMean(std::vector <double> vec);
    double  GetStandardDev(std::vector <double> vec);
    double  GetStandardErr(std::vector <double> vec);
    std::vector <double> GetFWHM(TH1D* h);
    TString FindData(TString directory);
    
};

#endif
