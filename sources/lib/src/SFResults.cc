// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              SFResults.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// ***************************************** 

#include "SFResults.hh"

ClassImp(SFResults);

//------------------------------------------------------------------
/// Table storing names of SFResultTypeNum enumerators.
TString gEnumNamesNum[] = {"kPeakConst", "kPeakPosition", "kPeakSigma", //SFPeakFinder
                           "kAveragePeakPos", //SFStability
                           "kTemp", //SFTemperature
                           "kLambda", "kChi2NDF", //SFAttenuation
                           "kEnergyRes", //SFEnergyRes
                           "kLight", //SFLightOutput
                           "kTimeRes", //SFTimingRes
                           "kPositionRes", //SFPositionRes
                           "kFastDecay", "kSlowDecay", "kIFast", "kISlow", //SFTimeConst
                           "kS0", "kEtaR", "kEtaL", "kKsi", "kLength", //SFAttenuationModel
                           "kAlpha", //SFEnergyreco
                           "kMLRSlope", "kMLROffset", "kACoeff", "kBCoeff", //SFPositionReco
                           "kCounts", "kCountsStdDev" //SFCountsMap
                          };
                           
/// Table storing names of SFResultTypeObj enumerators.                       
TString gEnumNamesObj[] = {"kSpectrum", "kPeak", //SFPeakFinder
                           "kPeakPosGraph", "kSMResidualGraph", //SFStability
                           "kAttGraph", "kMLRSigmaGraph", //SFAttenuation
                           "kEnergyResGraph", //SFEnergyRes
                           "kLightGraph", //SFLightOutput
                           "kTimeResGraph", "kTimeSigGraph", //SFTimingRes
                           "kPosRecoVsPosGraph", "kPosResVsPosGraph", "kPosVsMLRGraph", //SFPositionRes
                           "kResidualGraph", "kPositionAllHist",
                           "kPlFun", "kPrFun", "kRlFun", "kRrFun", "kSlFun", "kSrFun", //SFAttenuationModel
                           "kPlRecoFun", "kPrRecoFun", "kSlVsPosGraph", "kSrVsPosGraph",
                           "kPlVsPosGraph", "kPrVsPosGraph", 
                           "kAlphaGraph", "kEnergyRecoGraph", "kEnergyRecoFun", //SFEnergyReco
                           "kEnergyRecoSpecGraph", "kEnergyAllHist",
                           "kAGraph", //SFPositionreco
                           "kCountsGraph" //SFCountsMap
                          };
//------------------------------------------------------------------
/// Standard constructor.
/// \param name - name of the object
SFResults::SFResults(TString name) : fName(name)
{
}
//------------------------------------------------------------------
/// Default constructor. 
SFResults::SFResults() : fName("results")
{
}
//------------------------------------------------------------------
/// Default destructor.
SFResults::~SFResults()
{
    fValues.clear();
    fUncert.clear();
    fObjects.clear();
}
//------------------------------------------------------------------
/// Returns name of the enumerator SFResultTypeNum as a string.
/// \param id - enumerator
TString SFResults::EnumToString(SFResultTypeNum id)
{
    TString result_name = gEnumNamesNum[id];
    return result_name;
}
//------------------------------------------------------------------
/// Returns name of the enumerator SFResultTypeObj as a string.
/// \param id - enumerator
TString SFResults::EnumToString(SFResultTypeObj id)
{
    TString result_name = gEnumNamesObj[id];
    return result_name;
}
//------------------------------------------------------------------
/// Adds numerical result along with its uncertainty. 
/// \param id - enumerator
/// \param value - numerical result
/// \param error - uncertainty
void SFResults::AddResult(SFResultTypeNum id, double value, double error)
{
    if(fValues.count(id) != 0 ||
       fUncert.count(id) != 0)
    {
        std::cerr << "##### Warning in SFResults::AddResults()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " already exists!" << std::endl;
        return;
    }
    
    fValues.insert(std::pair <SFResultTypeNum, double> (id, value));
    fUncert.insert(std::pair <SFResultTypeNum, double> (id, error));
    
    return;
}
//------------------------------------------------------------------
/// Adds numerical result only (no uncertainty).
/// \param id - enumerator
/// \param value - numerical result
void SFResults::AddValue(SFResultTypeNum id, double value)
{
    if(fValues.count(id) != 0)
    {
        std::cerr << "##### Warning in SFResults::AddValue()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " already exists!" << std::endl;
        return;
    }
    
    fValues.insert(std::pair <SFResultTypeNum, double> (id, value));
    
    return;
}
//------------------------------------------------------------------
/// Adds uncertainty of the numerical result. 
/// \param id - enumerator 
/// \param error - uncertainty
void SFResults::AddUncertainty(SFResultTypeNum id, double error)
{
    if(fUncert.count(id) != 0)
    {
        std::cerr << "##### Warning in SFResults::AddUncertainty()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " already exists!" << std::endl;
        return;
    }
    
    fUncert.insert(std::pair <SFResultTypeNum, double> (id, error));
    
    return;
}
//------------------------------------------------------------------
/// Adds object-based result.
/// \param id - enumerator 
/// \param obj - pointer to the object (any object inheriting from TObject)
void SFResults::AddObject(SFResultTypeObj id, TObject *obj)
{
    if(fObjects.count(id) != 0)
    {
        std::cerr << "##### Warning in SFResults::AddObject()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " already exists!" << std::endl;
        return;
    }
    
    fObjects.insert(std::pair <SFResultTypeObj, TObject*> (id, obj));
    
    return;
}
//------------------------------------------------------------------
/// Returns numerical result along with its uncertainty in a std::vector:
/// [0] - result, [1] - uncertainty.
/// \param id - enumerator
std::vector <double> SFResults::GetResult(SFResultTypeNum id)
{
    std::vector <double> vec(2);
    
    if(fValues.count(id) == 0 ||
       fUncert.count(id) == 0)
    {
        std::cerr << "##### Warning in SFResults::GetResult()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " does not exist!" << std::endl;
        vec[0] = -1;
        vec[1] = -1;
    }
    
    vec[0] = fValues.find(id)->second;
    vec[1] = fUncert.find(id)->second;
    
    return vec;
}
//------------------------------------------------------------------
/// Returns numerical result only (no uncertainty).
/// \param id - enumerator
double SFResults::GetValue(SFResultTypeNum id)
{
    if(fValues.count(id) == 0)
    {
        std::cerr << "##### Warning in SFResults::GetValue()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " does not exist!" << std::endl;
        return -1;
    } 
    
    return fValues.find(id)->second;
}
//------------------------------------------------------------------
/// Returns uncertainty of the numerical result.
/// \param id - enumerator
double SFResults::GetUncertainty(SFResultTypeNum id)
{
    if(fUncert.count(id) == 0)
    {
        std::cerr << "##### Warning in SFResults::GetUncertainty()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " does not exist!" << std::endl;
        return -1;
    } 
    
    return fUncert.find(id)->second;
}
//------------------------------------------------------------------
/// Returns object-based result.
/// \param id - enumerator
TObject* SFResults::GetObject(SFResultTypeObj id)
{
    if(fObjects.count(id) == 0)
    {
        std::cerr << "##### Warning in SFResults::GetObject()!" << std::endl;
        std::cerr << "Element " << EnumToString(id) << " does not exist!" << std::endl;
        return nullptr;
    }
    
    return fObjects.find(id)->second;
}
//------------------------------------------------------------------
/// Prints details of the SFresults type object, including stored
/// results. 
void SFResults::Print(void)
{
    
    if(fValues.empty() &&
       fUncert.empty() &&
       fObjects.empty())
    {
        std::cerr << "##### Error in SFResults::Print()!" << std::endl;
        std::cerr << "Resuts not available!" << std::endl;
        return;
    }
    
    std::cout << "\n\n-------------------------------------------\n" << std::endl;
    std::cout << "SFResults object " << fName << "\n" << std::endl;
    std::cout << "Numerical results:" << std::endl;
    
    for(std::map <SFResultTypeNum, double>::iterator i = fValues.begin(); i != fValues.end() ; ++i)
    {
        std::cout << EnumToString(i->first) << " = " << i->second << " +/- " << fUncert.find(i->first)->second << std::endl;
    }
    
    std::cout << "\nObjects: " << std::endl;
    
    for(std::map <SFResultTypeObj, TObject*>::iterator i = fObjects.begin(); i != fObjects.end() ; ++i)
    {
        std::cout << EnumToString(i->first) << ": ";
        i->second->Print();
        std::cout << std::endl;
    }
    
    std::cout << "\n-------------------------------------------\n\n" << std::endl;
    
    return;
}
//------------------------------------------------------------------
