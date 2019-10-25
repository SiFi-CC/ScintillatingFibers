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
#include <sqlite3.h>

enum class SFTempSensor{
     Sensor_Ch01,
     Sensor_Ch10,
     Sensor_Ref,
     Sensor_Out
};

class SFTemperature : public TObject{
    
private:
    int     fSeriesNo;
    double  fAvTemp_446D45;
    double  fAvTemp_044F45;
    double  fAvTemp_2BAD44;
    double  fAvTemp_8F1F46;
    double  fAvTempErr_446D45;
    double  fAvTempErr_044F45;
    double  fAvTempErr_2BAD44;
    double  fAvTempErr_8F1F46;
    SFData *fData;
    TGraphErrors *fTempPlot_446D45;
    TGraphErrors *fTempPlot_044F45;
    TGraphErrors *fTempPlot_2BAD44;
    TGraphErrors *fTempPlot_8F1F46;
    TGraphErrors *fTempPlotAv_446D45;
    TGraphErrors *fTempPlotAv_044F45;
    TGraphErrors *fTempPlotAv_2BAD44;
    TGraphErrors *fTempPlotAv_8F1F46;
    
    std::vector <int> fTime;
    std::vector <double> fTemperatures;
    
    sqlite3 *fDB;
    
    bool OpenDataBase(void);
    TString GetSensorID(SFTempSensor sensor);
    bool LoadFromDB(SFTempSensor sensor, TString name = "none");
    std::vector <double> CalcAverageTempMeasure(SFTempSensor sensor, TString name);
    
public:
    SFTemperature(int seriesNo);
    ~SFTemperature();
    
    bool CalcAverageTempSeries(SFTempSensor sensor);
    bool BuildTempPlotAverage(SFTempSensor sensor);
    bool BuildTempPlot(SFTempSensor sensor);
    
    std::vector <double> GetAverageTempSeries(SFTempSensor sensor);
    TGraphErrors *GetTempPlot(SFTempSensor sensor);
    TGraphErrors *GetTempPlotAverage(SFTempSensor sensor);
    
    void Print(void);
    
    ClassDef(SFTemperature,1);
};

#endif
