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
#include "SFibersCal.h"
#include "SFibersRaw.h"
#include "SLoop.h"
#include "SFDrawCommands.hh"
//#include "SFTools.hh"
#include "SFInfo.hh"

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
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>

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

class SFData : public TObject
{
    
private:
    int                 fSeriesNo;
    SFInfo*             fInfo;
    std::vector<TFile*> fFiles; ///< Vector containing ROOT files with experimental data
 
    int gUnique = 0.;           ///< Unique flag to identify temporary histograms

    SFSignal* ConvertSignal(DDSignal* sig);
    SFSignal* ConvertSignal(SDDSignal* sig);
    bool      InterpretCut(SFSignal* sig, TString cut);
//     TProfile* GetSignalAverageKrakow(int ch, int ID, TString cut, int number, bool bl);
//     TProfile* GetSignalAverageAachen(int ch, int ID, TString cut, int number);
//     TH1D*     GetSignalKrakow(int ch, int ID, TString cut, int number, bool bl);
//     TH1D*     GetSignalAachen(int ch, int ID, TString cut, int number);

  public:
    SFData(int seriesNo);
    ~SFData();

    bool               OpenFiles(void);
//     TH1D*              GetSpectrum(SFChAddr ch,
//                                    SFSelectionType sel_type,
//                                    TString cut,
//                                    int ID);
// //     TH1D*              GetCustomHistogram(SFSelectionType sel_type,
// //                                           TString cut,
// //                                           int ID,
// //                                           std::vector<double> customNum = {});
//     TH1D*              GetCustomHistogram(SFChAddr ch,
//                                           SFSelectionType sel_type,
//                                           TString cut,
//                                           int ID,
//                                           std::vector<double> customNum);
//     TH2D*              GetCorrHistogram(SFSelectionType sel_type,
//                                         TString cut,
//                                         int ID, int ch = -1);
// //     TH2D*              GetRefCorrHistogram(int ID, int ch);
//     
//     std::vector<TH1D*> GetSpectra(int ch, SFSelectionType sel_type, TString cut);
//     std::vector<TH1D*> GetCustomHistograms(SFSelectionType sel_type, TString cut);
//     std::vector<TH2D*> GetCorrHistograms(SFSelectionType sel_type, TString cut, int ch = -1);
//     std::vector<TH2D*> GetRefCorrHistograms(int ch);
//     
//     TProfile*          GetSignalAverage(int ch, int ID, TString cut, int number, bool bl);
//     TH1D*              GetSignal(int ch, int ID, TString cut, int number, bool bl);

    ClassDef(SFData, 1)
};

#endif /* __SFData_H_ */
