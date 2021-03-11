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
TString gEnumNamesNum[] = {"kPeakConst", "kPeakPosition", "kPeakSigma", //SFPeakFinder
                           "kAveragePeakPos", //SFStability
                           "kTemp", //SFTemperature
                           "kLambda", "kChi2NDF", //SFAttenuation
                           "kEnergyRes", //SFEnergyRes
                           "kLight", //SFLightOutput
                           "kTimeRes", //SFTimingRes
                           "kPositionRes", //SFPositionRes
                           "kFastDecay", "kSlowDecay", "kIFast", "kISlow", //SFTimeConst
                           "kS0", "kMLambda", "kEtaR", "kEtaL", "kKsi", "kLength", //SFReconstrunction - model
                           "kAlpha", //SFReconstruction - energy reco
                           "kMLRSlope", "kMLROffset", "kACoeff", "kBCoeff" //SFReconstruction - position reco
                           };
                       
TString gEnumNamesObj[] = {"kSpectrum", "kPeak", //SFPeakFinder
                           "kPeakPosGraph", "kSMResidualGraph", //SFStability
                           "kAttGraph", "kMLRSigmaGraph", //SFAttenuation
                           "kEnergyResGraph", //SFEnergyRes
                           "kLightGraph", //SFLightOutput
                           "kTimeResGraph", "kTimeSigGraph", //SFTimingRes
                           "kPosRecoVsPosGraph", "kPosResVsPosGraph", "kMLRvsPosGraph", "kResidualGraph", //SFPositionRes
                           "kPlFun", "kPrFun", "kRlFun", "kRrFun", "kSlFun", "kSrFun" //SFReconstrunction - model
                           "kPlRecoFun", "kPrRecoFun", "kSlVsPosGraph", "kSrVsPosGraph",
                           "kPlVsPosGraph", "kPrVsPosGraph", // SFReconstruction - components
                           "kAlphaGraph", "kEnergyRecoGraph", "kEnergyRecoFun", "kEnergyDiffGraph", //SFReconstruction - energy reco
                           "kAGraph" //SFReconstruction - position reco
                       };
//------------------------------------------------------------------
SFResults::SFResults(TString name) : fName(name)
{
}
//------------------------------------------------------------------
SFResults::SFResults() : fName("results")
{
}
//------------------------------------------------------------------
SFResults::~SFResults()
{
    fValues.clear();
    fUncert.clear();
    fObjects.clear();
}
//------------------------------------------------------------------
TString SFResults::EnumToString(SFResultTypeNum id)
{
    TString result_name = gEnumNamesNum[id];
    return result_name;
}
//------------------------------------------------------------------
TString SFResults::EnumToString(SFResultTypeObj id)
{
    TString result_name = gEnumNamesObj[id];
    return result_name;
}
//------------------------------------------------------------------
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
