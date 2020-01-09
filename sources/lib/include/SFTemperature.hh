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
#include "TObject.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "SFData.hh"
#include "SFTools.hh"
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <map>
#include <sqlite3.h>


struct TemperatureResults{
  
    double fTemp    = -1;
    double fTempErr = -1;
};

class SFTemperature : public TObject{
    
private:
    int     fSeriesNo;
    SFData *fData;
    
    std::map <TString, TGraphErrors*> fTempPlot;
    std::map <TString, TGraphErrors*> fTempPlotAv;
    std::map <TString, TemperatureResults> fAvTemp;
    
    std::vector <int> fTime;
    std::vector <double> fTemperatures;
    
    sqlite3 *fDB;
    
    bool               OpenDataBase(void);
    bool               LoadFromDB(TString sensor, TString name = "none");
    TemperatureResults CalcAverageTempMeasure(TString sensor, TString name);
    
public:
    SFTemperature(int seriesNo);
    ~SFTemperature();
    
    bool CalcAverageTempSeries(TString sensor);
    bool BuildTempPlotAverage(TString sensor);
    bool BuildTempPlot(TString sensor);
    
    TemperatureResults GetAverageTempSeries(TString sensor);
    TGraphErrors      *GetTempPlot(TString sensor);
    TGraphErrors      *GetTempPlotAverage(TString sensor);
    
    void Print(void);
    
    ClassDef(SFTemperature,1);
};

#endif
