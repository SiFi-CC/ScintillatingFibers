// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFData.hh                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#ifndef __SFData_H_
#define __SFData_H_ 1

#include "DDSignal.hh"
#include "SCategoryManager.h"
#include "SDDSamples.h"
#include "SFDrawCommands.hh"
#include "SFTools.hh"
#include "SFibersCal.h"
#include "SFibersRaw.h"
#include "SLoop.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObject.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>
#include <TVectorT.h>

#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sqlite3.h>
#include <stdlib.h>
#include <string>
#include <vector>

/// \file
/// Structure representing recorded signal. It allows to convert
/// DDSignal and SDDSignal type obcjects into universal type for
/// the purpuse of the analysis. 

/// Structure representing recorded signal. It allows to convert
/// DDSignal and SDDSignal type obcjects into universal type for
/// the purpuse of the analysis. 

struct SFSignal
{

    float fAmp       = -100; ///< Amplitude
    float fCharge    = -100; ///< Uncalibrated charge (signal integral)
    float fPE        = -100; ///< Calibrated charge 
    float fT0        = -100; ///< Time T0
    float fTOT       = -100; ///< Time over threshold (TOT)
    float fBL        = -100; ///< Base line 
    float fBLsig     = -100; ///< Sigma of the base line
    int   fPileUp    = -100; ///< Pile-up flag
    int   fVeto      = -100; ///< Veto flag
    bool  fSDDSignal = 0;    ///< Flag indicatig whether contained information
                             ///< comes from SDDSignal (true) or DDSignal (false)

    /// Prints details of the SFSignal class object.
    void Print() const
    {
        printf("Amp=%f  Charge=%f  PE=%f  T0=%f  TOT=%f "
               " BL=%f BLsig=%f FU=%d VETO=%d SDDSig=%d\n",
               fAmp, fCharge, fPE, fT0, fTOT, fBL, fBLsig,
               fPileUp, fVeto, fSDDSignal);
    }
};

/// Class to access experiemntal data. Information about an experimental
/// series and all measurements is loaded from the SQLite3 data base.
/// Subsequently requested data is accessed from ROOT files and binary
/// files stored on the computer or remote scratch drive. It is possible
/// to get single histograms as well as vectors of histograms of requested
/// type. All possible selections for histograms are defined in SFDrawCommands
/// class. Single signals and averaged signals are also possible to access.
/// Cutting functionality based on ROOT's TTree for all has been implemented.

class SFData : public TObject
{

  private:
    int      fSeriesNo;      ///< Experimental series number
    int      fNpoints;       ///< Number of measurements in the series
    int      fAnalysisGroup; ///< Analysis group
    TString  fFiber;         ///< Scintillating fiber type e.g. LuAG:Ce Crytur (1)
    double   fFiberLength;   ///< Length of the scintillating fiber [mm]
    TString  fSource;        ///< Type of the radioactive source
    TString  fCollimator;    ///< Type of the used collimator: Lead or Electronic
    TString  fDesc;          ///< Description of the measurement series
    TString  fTestBench;     ///< Type of test bench: PL/DE/Simulation
    TString  fSiPM;          ///< SiPM type: Hamamatsu or SensL
    double   fOvervoltage;   ///< Overvoltage [V]
    TString  fCoupling;      ///< Coupling type: silicone gel/silicone pads
    TString  fLogFile;       ///< Name of measurment log file
    TString  fTempFile;      ///< Name of temperature log file
    sqlite3* fDB;            ///< SQLite3 data base

    std::vector<TFile*>  fFiles;     ///< Vector containing ROOT files with experimental data
    std::vector<TString> fNames;     ///< Vector containing names of measurements
    std::vector<double>  fPositions; ///< Vector containing positions of radioactive source [mm]
    std::vector<int>     fMeasureID; ///< Vector containing IDs of measurements
    std::vector<int>     fTimes;     ///< Vector containing times of measurements [s]
    std::vector<int>     fStart;     ///< Vector containing starting times of measurements
                                     ///< (in UNIX time)
    std::vector<int> fStop;          ///< Vector containing stopping times of measurements 
                                     ///< (in UNIX time)

    int gUnique = 0.;                ///< Unique flag to identify temporary histograms

    SFSignal* ConvertSignal(DDSignal* sig);
    SFSignal* ConvertSignal(SDDSignal* sig);
    bool      InterpretCut(SFSignal* sig, TString cut);
    TProfile* GetSignalAverageKrakow(int ch, int ID, TString cut, int number, bool bl);
    TProfile* GetSignalAverageAachen(int ch, int ID, TString cut, int number);
    TH1D*     GetSignalKrakow(int ch, int ID, TString cut, int number, bool bl);
    TH1D*     GetSignalAachen(int ch, int ID, TString cut, int number);

  public:
    SFData();
    SFData(int seriesNo);
    ~SFData();

    bool               OpenFiles(void);
    bool               OpenDataBase(TString name);
    bool               SetDetails(int seriesNo);
    SLoop*             GetTree(int ID);
    TH1D*              GetSpectrum(int ch, SFSelectionType sel_type, TString cut, int ID);
    TH1D*              GetCustomHistogram(SFSelectionType sel_type, TString cut, int ID,
                                          std::vector<double> customNum = {});
    TH1D*              GetCustomHistogram(int ch, SFSelectionType sel_type, TString cut, int ID,
                                          std::vector<double> customNum);
    TH2D*              GetCorrHistogram(SFSelectionType sel_type, TString cut, int ID, int ch = -1);
    TH2D*              GetRefCorrHistogram(int ID, int ch);
    std::vector<TH1D*> GetSpectra(int ch, SFSelectionType sel_type, TString cut);
    std::vector<TH1D*> GetCustomHistograms(SFSelectionType sel_type, TString cut);
    std::vector<TH2D*> GetCorrHistograms(SFSelectionType sel_type, TString cut, int ch = -1);
    std::vector<TH2D*> GetRefCorrHistograms(int ch);
    TProfile*          GetSignalAverage(int ch, int ID, TString cut, int number, bool bl);
    TH1D*              GetSignal(int ch, int ID, TString cut, int number, bool bl);
    void               Print(void);

    /// Returns number of measurements in the series.
    int GetNpoints(void) { return fNpoints; };
    /// Returns analysis group number.
    int GetAnalysisGroup(void) { return fAnalysisGroup; };
    /// Returns fiber type.
    TString GetFiber(void) { return fFiber; };
    /// Returns length of the fiber [mm]
    double GetFiberLength(void) { return fFiberLength; };
    /// Returns description of the series.
    TString GetDescription(void) { return fDesc; };
    /// Returns collimator type.
    TString GetCollimator(void) { return fCollimator; };
    /// Returns source type.
    TString GetSource(void) { return fSource; };
    /// Returns test bench type.
    TString GetTestBench(void) { return fTestBench; };
    /// Returns SiPM type.
    TString GetSiPM(void) { return fSiPM; };
    /// Returns overvoltage [V].
    double GetOvervoltage(void) { return fOvervoltage; };
    /// Returns coupling type.
    TString GetCoupling(void) { return fCoupling; };
    /// Returns name of the measurment log file.
    TString GetLogFile(void) { return fLogFile; };
    /// Returns name of the temperature log file.
    TString GetTempFile(void) { return fTempFile; };
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

    ClassDef(SFData, 1)
};

#endif /* __SFData_H_ */
