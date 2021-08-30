// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             SFCountMap.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFCountMap.hh"

ClassImp(SFCountMap);

//------------------------------------------------------------------
/// Standard constructor.
/// \param seriesNo is a number of experimental series to be analyzed.
SFCountMap::SFCountMap(int seriesNo) : fSeriesNo(seriesNo),
                                       fData(nullptr),
                                       fResultsCh0(nullptr),
                                       fResultsCh1(nullptr)
{
    
    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cout << message << std::endl;
        throw "##### Exception in SFCountMap constructor!";
    }
    
    TString testBench   = fData->GetTestBench();
    TString description = fData->GetDescription();
    
    if (testBench != "PMI")
    {
        std::cout << "##### Warning in SFCountMap constructor!" << std::endl;
        std::cout << "This analysis was meant for PMI data!" << std::endl;
        std::cout << "test bench: " << testBench << std::endl;
    }
    
    if(! description.Contains("Module series"))
    {
        std::cout << "##### Warning in SFCountMap constructor!" << std::endl;
        std::cout << "This analysis was meant for module series!" << std::endl;
        std::cout << "description: " << description << std::endl;
    }
    
    TString cutCh0 = SFDrawCommands::GetCut(SFCutType::kPMISpecCh0);
    TString cutCh1 = SFDrawCommands::GetCut(SFCutType::kPMISpecCh1);
    fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::kPMICharge, cutCh0);
    fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::kPMICharge, cutCh1);
    
    fResultsCh0 = new SFResults(Form("CountMapResults_S%i_ch0", fSeriesNo));
    fResultsCh1 = new SFResults(Form("CountMapResults_S%i_ch1", fSeriesNo));
}
//------------------------------------------------------------------
/// Default destructor.
SFCountMap::~SFCountMap()
{
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
/// Creates graph of counts registered in the fibers as a function of 
/// the fiber number for a chosen fiber side (channel 0/1). Number of 
/// counts is obtained as an integral charge spectra. IntegralAndError()
/// function of ROOT framework is used for that purpose. Average number
/// of counts is calulated as weighted mean with standard error of the 
/// weighted mean for each channel. Additionally, standard deviation of 
/// number of counts is calculated.
/// \param ch is a channel number (0 or 1)
bool SFCountMap::DrawCounts(int ch)
{
    std::cout << "\n----- Inside SFCountsMap::DrawCounts()" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << " ch: " << ch << std::endl;
    
    int npoints                 = fData->GetNpoints();
    TString collimator          = fData->GetCollimator();
    TString testBench           = fData->GetTestBench();
    std::vector<double> fiberNo = fData->GetPositions();
    
    std::vector<TH1D*> spectra;
    
    if (ch == 0)
        spectra = fSpectraCh0;
    else if (ch == 1)
        spectra = fSpectraCh1;
    else
    {
        std::cerr << "##### Error in SFCountMap::DrawCounts()!" << std::endl;
        std::cerr << "Incorrect channel number: " << ch << std::endl;
        return false;
    }
    
    TGraphErrors *graph = new TGraphErrors(npoints);
    graph->SetName(Form("CountsMap_S%i_ch%i", fSeriesNo, ch));
    graph->SetTitle(Form("Counts per fiber S%i", fSeriesNo));
    graph->GetXaxis()->SetTitle("fiber number");
    graph->GetYaxis()->SetTitle("counts");
    graph->SetMarkerStyle(4);
    
    double integral    = 0;
    double integralErr = 0;
    int    nbins       = 0;
    
    std::vector<double> integralsAll;
    
    double countsAv    = 0;
    double countsAvErr = 0;
    
    for (int i=0; i<npoints; i++)
    {
        nbins = spectra[i]->GetNbinsX();
        integral = spectra[i]->IntegralAndError(1, nbins, integralErr);
        graph->SetPoint(i, fiberNo[i], integral);
        graph->SetPointError(i, 0, integralErr);
        //std::cout << "nbins: " << nbins << " integral: " << integral 
        //          << " error: " << integralErr << std::endl;
        countsAv += integral * (1./pow(integralErr, 2));
        countsAvErr += 1./pow(integralErr, 2);
        integralsAll.push_back(integral);
    }
    
    countsAv    = countsAv / countsAvErr;
    countsAvErr = sqrt(1. / countsAvErr);
    
    double stdDev = SFTools::GetStandardDev(integralsAll);
    
    if(ch == 0)
    {
        fResultsCh0->AddObject(SFResultTypeObj::kCountsGraph, graph);
        fResultsCh0->AddResult(SFResultTypeNum::kCounts, countsAv, countsAvErr);
        fResultsCh0->AddResult(SFResultTypeNum::kCountsStdDev, stdDev, 0);
    }
    else if (ch ==1)
    {
        fResultsCh1->AddObject(SFResultTypeObj::kCountsGraph, graph);
        fResultsCh1->AddResult(SFResultTypeNum::kCounts, countsAv, countsAvErr);
        fResultsCh1->AddResult(SFResultTypeNum::kCountsStdDev, stdDev, 0);
    }
    
    std::cout << "\n\tAveraged counts: " << countsAv << " +/- " << countsAvErr << std::endl;
    std::cout << "\tStandard deviation: " << stdDev << std::endl << std::endl;
    
    return true;
}
//------------------------------------------------------------------
/// Returns analyzed charge spectra for requested channel.
/// \param ch is a channel number (0 or 1)
std::vector<TH1D*> SFCountMap::GetSpectra(int ch)
{
    if ((ch == 0 && fSpectraCh0.empty()) ||
        (ch == 1 && fSpectraCh1.empty()))
    {
        std::cerr << "##### Error in SFCountMap::GetSpectra() for ch" << ch << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }
        
    if (ch == 0)
        return fSpectraCh0;
    else if(ch == 1)
        return fSpectraCh1;
    else
    {
        std::cerr << "##### Error in SFCountMap::GetSpectra for ch" << ch << std::endl;
        std::cerr << "Incorrect channel number!" << std::endl;
        std::abort();
    }
}
//------------------------------------------------------------------
/// Returns vector containing SFResults type objects with the results
/// of the analysis for both channels (channel 0 at index 0 and channel
/// 1 at index 1) 
std::vector<SFResults*> SFCountMap::GetResults(void)
{
    if(fResultsCh0 == nullptr ||
       fResultsCh1 == nullptr)
    {
        std::cerr << "##### Error in SFCountMap::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(2);
    results[0] = fResultsCh0;
    results[1] = fResultsCh1;
    
    return results;
}
//------------------------------------------------------------------
/// Prints details of the SFCountMap class object.
void SFCountMap::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFCountMap class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
