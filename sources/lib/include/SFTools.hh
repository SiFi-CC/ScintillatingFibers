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

#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>

#include "SFData.hh"

#include <iostream>
#include <sqlite3.h>
#include <stdlib.h>
#include <time.h>

/// Namespace containing functions useful for data analysis, e.g.
/// getting details of the measurement series, accessing and modyfing
/// data bases and data, calculating statistical parameters, fitting, etc.

namespace SFTools
{
    
int                 GetIndex(std::vector<int> measurementsIDs, int id);
int                 GetSeriesNo(TString hname_tstr);
int                 GetChannel(TString hname_tstr);
double              GetPosition(TString hname_tstr);
int                 GetMeasurementID(TString hname_tstr);
int                 GetMeasurementID(int seriesNo, double position);
double              GetPosError(TString collimator, TString testBench);
double              GetSigmaBL(TString SiPM);
bool                CheckDBStatus(int status, sqlite3* database);
bool                SaveResultsDB(TString database, TString table, TString query, int seriesNo);
bool                CreateTable(TString database, TString table);
double              GetMean(std::vector<double> vec);
double              GetStandardDev(std::vector<double> vec);
double              GetStandardErr(std::vector<double> vec);
double              FindMaxXaxis(TH1D* h);
double              FindMaxYaxis(TH1D* h);
bool                RatiosFitGauss(std::vector<TH1D*>& vec, float range_in_RMS = 1);
bool                RatiosFitDoubleGauss(std::vector<TH1D*>& vec, float range_in_RMS = 1);
bool                FitGaussSingle(TH1D* h, float range_in_RMS);
TString             FindData(TString directory);
std::vector<double> GetFWHM(TH1D* h);

};

#endif /* __SFTools_H_ */
