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

/// \file 
/// Enumeration representing different types of 
/// numerical results.

/// Enumeration representing different types of 
/// numerical results.

enum SFResultTypeNum
{
    //----- SFPeakFinder
    kPeakConst,      ///< Height oh the 511 keV peak
    kPeakPosition,   ///< Position of the 511 keV peak
    kPeakSigma,      ///< Sigma of the 511 keV peak

    //----- SFStabilityMon
    kAveragePeakPos, ///< Average peak position in stability monitoring series
    
    //----- SFTemperature
    kTemp,           ///< Temperature

    //----- SFAttenuation
    kLambda,         ///< Attenuation length
    kChi2NDF,        ///< Chi2/NDF of the fit

    //----- SFEnergyRes
    kEnergyRes,      ///< Energy resolution

    //----- SFLightOutput
    kLight,          ///< Light collection/Light output

    //----- SFTimingRes
    kTimeRes,        ///< Timing resolution

    //----- SFPositionRes
    kPositionRes,    ///< Position resolution
    kPositionReco,   ///< Reconstructed position
    
    //----- SFTimeConst
    kFastDecay,      ///< Fast decay constant
    kSlowDecay,      ///< Slow decay constant
    kIFast,          ///< Intensity of the fast decay
    kISlow,          ///< Intensity of the slow decay

    //----- SFAttenuationModel
    kS0,             ///< Model parameter: amplitude
    kEtaR,           ///< Model parameter: reflection coefficient right (ch1)
    kEtaL,           ///< Model parameter: reflection coefficient left (ch0)
    kKsi,            ///< Model parameter: asymmetry coefficient
    kLength,         ///< Model parameter: fiber length

    //----- SFEnergyReco
    kAlpha,          ///< Alpga coefficient for energy reconstruction

    //----- SFPositionReco
    kMLRSlope,       ///< Slope of the corrected MLR curve
    kMLROffset,      ///< Offset of the corrected MLR curve
    kACoeff,         ///< A coefficient for position reconstruction
    kBCoeff,         ///< B coefficient for position reconstruction
    
    //----- counts map
    kCounts,        ///< Average number of counts 
    kCountsStdDev   ///< Standard deviation of number of counts
};

/// \file 
/// Enumeration representing different types of 
/// object-based results.

/// Enumeration representing different types of 
/// object-based results.

enum SFResultTypeObj
{
    //----- SFPeakFinder
    kSpectrum,          ///< Fitted charge spectrum
    kPeak,              ///< Background-subtracted charge spectrum

    //----- SFStabilityMon
    kPeakPosGraph,      ///< Stability monitoring graph i.e. 511 keV 
                        ///< peak position vs. source position
    kSMResidualGraph,   ///< Stability monitoring residuals graph

    //----- SFAttenuation
    kAttGraph,          ///< Attenuation graph
    kMLRSigmaGraph,     ///< Sigma of the MLR distribution vs. source position

    //----- SFEnergyRes
    kEnergyResGraph,    ///< Energy resolution graph

    //----- SFLightOutput
    kLightGraph,        ///< Light collection/light output graph 

    //----- SFTimingRes
    kTimeResGraph,      ///< Mean of the TD distribution vs. source position 
    kTimeSigGraph,      ///< Sigma of the TD distribution vs. source position

    //----- SFPositionRes
    kPosRecoVsPosGraph, ///< Reconstructed position vs. true source position
    kPosResVsPosGraph,  ///< Position resolution vs. source position
    kPosVsMLRGraph,     ///< Reversed MLR curve
    kResidualGraph,     ///< Reconstructed position residuals graph
    kPositionDist,      ///< Summed reconstructed position distribution

    //----- SFAttenuationModel
    kPlFun,             ///< Function representing primary light emission (left/ch0)
    kPrFun,             ///< Function representing primary light emission (right/ch1)
    kRlFun,             ///< Function representing reflected light component (left/ch0)
    kRrFun,             ///< Function representing reflected light component (right/ch1)
    kSlFun,             ///< Function representing total measured signal (left/ch0)
    kSrFun,             ///< Function representing total measured signal (right/ch1)
    kPlRecoFun,         ///< Function of the reconstructed primary light component (left)
    kPrRecoFun,         ///< Function of the reconstructed primary light component (right)
    kSlVsPosGraph,      ///< 511 keV peak position (experimental) vs. source position (left)
    kSrVsPosGraph,      ///< 511 keV peak position (experimental) vs. source position (right)
    kPlVsPosGraph,      ///< 511 keV peak position (corrected) vs. source position (left)  
    kPrVsPosGraph,      ///< 511 keV peak position (corrected) vs. source position (right)

    //----- SFEnergyReco
    kAlphaGraph,          ///< Alpha coefficient (for energy reconstruction) vs. source position
    kEnergyRecoGraph,     ///< 511 keV peak position (from attenuation curves) vs. source
                          ///< position
    kEnergyRecoFun,       ///< Reconstructed energy function
    kEnergyRecoSpecGraph, ///< 511 keV peak position (from energy-reconstructed spectra) vs.
                          ///< source position
    kEnergyAllHist,       ///< Summed energy-reconstructed spectrum

    //----- SFPositionReco
    kAGraph,              ///< A coeffcient (for position reconstruction) vs. source position
    
    //----- SFCountMap
    kCountsGraph          ///< Counts vs. fiber number
};

/// Container class which allows to store results of the analysis. Each 
/// SFResults type object contains three maps storing numerical values, 
/// uncertainties of numerical values and objects. Therefore it is possible 
/// to store numerical results along with their uncertainties and object-
/// based results e.g. histograms, graphs, profiles etc. Storing containers 
/// (e.g. vectors or tables) is not supported. In order to save results use
/// keys predefined as enumerators SFResultsTypeNum and SFResultsTypeObj.

class SFResults
{

  private:
    TString                             fName;    ///< Name of the object
    std::map<SFResultTypeNum, double>   fValues;  ///< Map containing numerical results
    std::map<SFResultTypeNum, double>   fUncert;  ///< Map containing uncertainties of 
                                                  ///< numerical results 
    std::map<SFResultTypeObj, TObject*> fObjects; ///< Map containing object-based results

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

    /// Sets the name of the object
    void    SetName(TString name) { fName = name; };
    /// Returns the name of the object
    TString GetName(void)         { return fName; };

    void Print(void);
};

#endif /* __SFResults_H_ */
