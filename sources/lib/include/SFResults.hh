// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFResults.hh             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#ifndef __SFResults_H_
#define __SFResults_H_ 1

#include <TObject.h>
#include <TString.h>

#include <iostream>
#include <map>
#include <set>

enum SFResultTypeNum
{
    //----- SFPeakFinder
    kPeakConst,
    kPeakPosition,
    kPeakSigma,

    //----- SFStabilityMon
    kAveragePeakPos,
    
    //----- SFTemperature
    kTemp,

    //----- SFAttenuation
    kLambda,
    kChi2NDF,

    //----- SFEnergyRes
    kEnergyRes,

    //----- SFLightOutput
    kLight,

    //----- SFTimingRes
    kTimeRes,

    //----- SFPositionRes
    kPositionRes,
    
    //----- SFTimeConst
    kFastDecay,
    kSlowDecay,
    kIFast,
    kISlow,

    //----- SFReconstruction
    //----- light attenuation model
    kS0,      ///< Model parameter: amplitude
    kMLambda, ///< Model parameter: attenuation length
    kEtaR,    ///< Model parameter: reflection coefficient right (ch1)
    kEtaL,    ///< Model parameter: reflection coefficient left (ch0)
    kKsi,     ///< Model parameter: asymmetry coefficient
    kLength,  ///< Model parameter: fiber length

    //----- energy reconstruction
    kAlpha,

    //----- position reconstruction
    kMLRSlope,
    kMLROffset,
};

enum SFResultTypeObj
{
    //----- SFPeakFinder
    kSpectrum,
    kPeak,

    //----- SFStabilityMon
    kPeakPosGraph,
    kSMResidualGraph,

    //----- SFAttenuation
    kAttGraph,
    kMLRSigmaGraph,

    //----- SFEnergyRes
    kEnergyResGraph,

    //----- SFLightOutput
    kLightGraph,

    //----- SFTimingRes
    kTimeResGraph,
    kTimeSigGraph,

    //----- SFPositionRes
    kPosRecoVsPosGraph,
    kPosResVsPosGraph,
    kMLRvsPosGraph,
    kPRResidualGraph,

    //----- SFReconstruction
    //----- light attenuation model
    kPlFun,     ///< Function representing primary light emission (left/ch0)
    kPrFun,     ///< Function representing primary light emission (right/ch1)
    kRlFun,     ///< Function representing reflected light component (left/ch0)
    kRrFun,     ///< Function representing reflected light component (right/ch1)
    kSlFun,     ///< Function representing total measured signal (left/ch0)
    kSrFun,     ///< Function representing total measured signal (right/ch1)

    //----- signal & reconstructed direct component
    kPlRecoFun,
    kPrRecoFun,
    kSlVsPosGraph,
    kSrVsPosGraph,
    kPlVsPosGraph,
    kPrVsPosGraph,
    kSMLRVsPosGraph,
    kPMLRVsPosGraph,

    //----- energy reconstruction
    kAlphaGraph,
    kEnergyRecoGraph,
    kEnergyRecoFun,
    kEnergyDiffGraph,

    //----- position reconstruction
    kMPosRecoVsPosGraph
};

class SFResults : public TObject
{

  private:
    TString                             fName;
    std::map<SFResultTypeNum, double>   fValues;
    std::map<SFResultTypeNum, double>   fUncert;
    std::map<SFResultTypeObj, TObject*> fObjects;

  public:
    SFResults(TString name);
    SFResults();
    ~SFResults();

    void AddResult(SFResultTypeNum id, double value, double error);
    void AddValue(SFResultTypeNum id, double vlaue);
    void AddUncertainty(SFResultTypeNum id, double error);
    void AddObject(SFResultTypeObj id, TObject* obj);

    std::vector<double> GetResult(SFResultTypeNum id);
    double              GetValue(SFResultTypeNum id);
    double              GetUncertainty(SFResultTypeNum id);
    TObject*            GetObject(SFResultTypeObj id);

    TString EnumToString(SFResultTypeNum id);
    TString EnumToString(SFResultTypeObj id);

    void    SetName(TString name) { fName = name; };
    TString GetName(void)         { return fName; };

    void Print(void);

    ClassDef(SFResults, 1)
};

#endif
