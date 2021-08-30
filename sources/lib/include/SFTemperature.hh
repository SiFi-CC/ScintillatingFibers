// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFTemperature.hh            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#ifndef __SFTemperature_H_
#define __SFTemperature_H_ 1
#include "SFData.hh"
#include "SFTools.hh"
#include "SFResults.hh"
#include <TGraphErrors.h>
#include <TMath.h>
#include <TObject.h>
#include <iostream>
#include <map>
#include <sqlite3.h>
#include <stdlib.h>
#include <string>
#include <vector>

/// This class prepares temperature plots and calculates average
/// temperature for the experimental series. Temperature data is 
/// read from the data base. Plots of temperature per minute and 
/// average temperature during the measurement are prepared. Plots
/// and calculations of average temeratures are done for requested 
/// temperature sensors, identified by their ID.

class SFTemperature : public TObject
{

  private:
    int     fSeriesNo; ///< Series number
    SFData* fData;     ///< Analyzed experimental series

    std::map<TString, TGraphErrors*> fTempPlot;   ///< Map containing temperature plots
                                                  ///< (key - sensor ID)
    std::map<TString, TGraphErrors*> fTempPlotAv; ///< Map containing plots of average
                                                  ///< temperature during measurement
                                                  ///< (key - sensor ID)
    std::map<TString, SFResults*>    fAvTemp;     ///< Map containing average temperature 
                                                  ///< for experimental series (key - 
                                                  ///< sensor ID)

    std::vector<int>    fTime;         ///< Vector containing timestamp of
                                       ///< the temperature measurement
    std::vector<double> fTemperatures; ///< Vector containing temperature values

    sqlite3* fDB;                      ///< SQLite3 data base

    bool       OpenDataBase(void);
    bool       LoadFromDB(TString sensor, TString name = "all");
    SFResults* CalcAverageTempMeasure(TString sensor, TString name);

  public:
    SFTemperature(int seriesNo);
    ~SFTemperature();

    bool CalcAverageTempSeries(TString sensor);
    bool BuildTempPlotAverage(TString sensor);
    bool BuildTempPlot(TString sensor);

    SFResults*    GetAverageTempSeries(TString sensor);
    TGraphErrors* GetTempPlot(TString sensor);
    TGraphErrors* GetTempPlotAverage(TString sensor);

    void Print(void);

    ClassDef(SFTemperature, 1)
};

#endif
