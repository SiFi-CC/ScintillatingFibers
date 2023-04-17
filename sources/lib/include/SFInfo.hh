// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFInfo.hh                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#ifndef __SFInfo_H_
#define __SFInfo_H_ 1

#include <TObject.h>
#include <TString.h>

#include <iostream>
#include <iomanip>
#include <sqlite3.h>
#include <vector>
#include <stdlib.h>

class SFInfo : public TObject
{

private:
    int      fSeriesNo;      ///< Experimental series number
    int      fNpoints;       ///< Number of measurements in the series
    double   fFiberLength;   ///< Length of the scintillating fiber [mm]
    double   fOvervoltage;   ///< Overvoltage [V]
    TString  fFiber;         ///< Scintillating fiber type e.g. LuAG:Ce Crytur (1)
    TString  fSource;        ///< Type of the radioactive source
    TString  fCollimator;    ///< Type of the used collimator: Lead or Electronic
    TString  fDesc;          ///< Description of the measurement series
    TString  fTestBench;     ///< Type of test bench: PL/DE/Simulation
    TString  fSiPM;          ///< SiPM type: Hamamatsu or SensL
    TString  fCoupling;      ///< Coupling type: silicone gel/silicone pads
    TString  fLogFile;       ///< Name of measurment log file
    TString  fTempFile;      ///< Name of temperature log file
    TString  fDAQ;           ///< Data acquisition system
    
    sqlite3* fDB;            ///< SQLite3 data base

    std::vector<TString> fNames;     ///< Vector containing names of measurements
    std::vector<double>  fPositions; ///< Vector containing positions of radioactive source [mm]
    std::vector<int>     fMeasureID; ///< Vector containing IDs of measurements
    std::vector<int>     fTimes;     ///< Vector containing times of measurements [s]
    std::vector<int>     fStart;     ///< Vector containing starting times of 
                                     ///< measurements (in UNIX time)
    std::vector<int>     fStop;      ///< Vector containing stopping times of 
                                     ///< measurements (in UNIX time)
                                     
public:
    SFInfo(int seriesNo);
    SFInfo();
    ~SFInfo();
    
    /// Returns number of measurements in the series.
    int GetNpoints(void) { return fNpoints; };
    /// Returns length of the fiber [mm]
    double GetFiberLength(void) { return fFiberLength; };
    /// Returns overvoltage [V].
    double GetOvervoltage(void) { return fOvervoltage; };
    /// Returns fiber type.
    TString GetFiber(void) { return fFiber; };
    /// Returns source type.
    TString GetSource(void) { return fSource; };
    /// Returns collimator type.
    TString GetCollimator(void) { return fCollimator; };
    /// Returns description of the series.
    TString GetDescription(void) { return fDesc; };
    /// Returns test bench type.
    TString GetTestBench(void) { return fTestBench; };
    /// Returns SiPM type.
    TString GetSiPM(void) { return fSiPM; };
    /// Returns coupling type.
    TString GetCoupling(void) { return fCoupling; };
    /// Returns name of the measurment log file.
    TString GetLogFile(void) { return fLogFile; };
    /// Returns name of the temperature log file.
    TString GetTempFile(void) { return fTempFile; };
    /// Returns DAQ type.
    TString GetDAQ(void) { return fDAQ; };
    /// Returns a vector containing names of all measurements in the series.
    std::vector<TString> GetNames(void) { return fNames; };
    /// Returns a vector containing source positions for all measurements in the series.
    std::vector<double> GetPositions(void) { return fPositions; };
    /// Returns a vector containing measurement times for all measurements in the series.
    std::vector<int> GetTimes(void) { return fTimes; };
    /// Returns a vector containing starting times of all measurements in the series.
    std::vector<int> GetStartTimes(void) { return fStart; };
    /// Returns a vector containing stopping times of all measurements in the series.
    std::vector<int> GetStopTimes(void) { return fStop; };
    /// Returns a vector containing IDs of all measurements in the series.
    std::vector<int> GetMeasurementsIDs(void) { return fMeasureID; };
    
    void SetNpoints (int npoints) { fNpoints = npoints; };
    void SetFiberLength (double fiberLen) { fFiberLength = fiberLen; };
    void SetOvervoltage (double overvol) { fOvervoltage = overvol; };
    void SetFiber (TString fiber) { fFiber = fiber; };
    void SetSource (TString source) { fSource = source; };
    void SetCollimator (TString collimator) { fCollimator = collimator; };
    void SetDescription (TString desc) { fDesc = desc; };
    void SetTestBench (TString testBench) { fTestBench = testBench; };
    void SetSiPM (TString sipm) { fSiPM = sipm; };
    void SetCoupling (TString coupling) { fCoupling = coupling; };
    void SetLogFile (TString logfile) { fLogFile = logfile; };
    void SetTempFile (TString tempfile) { fTempFile = tempfile; };
    void SetDAQ (TString daq) { fDAQ = daq; };
    void SetNames (std::vector<TString> names) { fNames = names; };
    void SetPositions (std::vector<double> positions) { fPositions = positions; };
    void SetTimes (std::vector<int> times) { fTimes = times; };
    void SetStartTimes (std::vector<int> start) { fStart = start; };
    void SetStopTimes (std::vector<int> stop) { fStop = stop; };
    void SetMeasurementIDs (std::vector<int> ids) { fMeasureID = ids; };
     
    void Print(void);

//     ClassDef(SFInfo, 1)
};

#endif /* __SFInfo_H_ */
