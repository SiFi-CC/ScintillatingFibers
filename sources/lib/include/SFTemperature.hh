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

class SFTemperature : public TObject
{

  private:
    int     fSeriesNo;
    SFData* fData;

    std::map<TString, TGraphErrors*> fTempPlot;
    std::map<TString, TGraphErrors*> fTempPlotAv;
    std::map<TString, SFResults*>    fAvTemp;

    std::vector<int>    fTime;
    std::vector<double> fTemperatures;

    sqlite3* fDB;

    bool       OpenDataBase(void);
    bool       LoadFromDB(TString sensor, TString name = "none");
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

    ClassDef(SFTemperature, 1);
};

#endif
