// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFEnergyRes.cc              *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyRes.hh"
#include <cmath>

ClassImp(SFEnergyRes);

//------------------------------------------------------------------
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed.
SFEnergyRes::SFEnergyRes(int seriesNo) : fSeriesNo(seriesNo),
                                         fData(nullptr),
                                         fEnergyResCh0Graph(nullptr),
                                         fEnergyResCh1Graph(nullptr),
                                         fEnergyResAveGraph(nullptr),
                                         fResultsCh0(nullptr),
                                         fResultsCh1(nullptr),
                                         fResultsAve(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFEnergyRes constructor!";
    }

    TString desc = fData->GetDescription();
    if (!desc.Contains("Regular series"))
    {
        std::cout << "##### Warning in SFEnergyRes constructor!" << std::endl;
        std::cout << "Calculating energy resolution with non-regular series!" << std::endl;
    }

    //--- getting necessary PE spectra
    double              s          = SFTools::GetSigmaBL(fData->GetSiPM());
    std::vector<double> sigmas     = {s, s};
    TString             cut_ch0    = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigmas);
    TString             cut_ch1    = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigmas);
    TString             cut_ch0ch1 = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);
    fSpectraCh0                    = fData->GetSpectra(0, SFSelectionType::kPE, cut_ch0);
    fSpectraCh1                    = fData->GetSpectra(1, SFSelectionType::kPE, cut_ch1);
    fSpectraAve = fData->GetCustomHistograms(SFSelectionType::kPEAverage, cut_ch0ch1);

    fResultsCh0 = new SFResults(Form("EnergyResResults_S%i_ch0", fSeriesNo));
    fResultsCh1 = new SFResults(Form("EnergyResResults_S%i_ch1", fSeriesNo));
    fResultsAve = new SFResults(Form("EnergyResResults_S%i_ave", fSeriesNo));
}
//------------------------------------------------------------------
/// Default destructor.
SFEnergyRes::~SFEnergyRes()
{
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
/// Calculates energy resolution based on charge spectra of requested channel
/// for all measurements in analyzed measurement series. In this function
/// fEnergyResGrapCh0 and fEnergyResGraphCh1 graphs are filled and values of
/// fResults.fEnergyResCh0, fEnergyResCh0Err, fResults.fEnergyResCh1 and
/// fResults.fEnergyResAveCh1Err are assigned.
/// \param ch - channel number
bool SFEnergyRes::CalculateEnergyRes(int ch)
{

    std::cout << "\n----- Inside SFEnergyRes::CalculateEnergyRes()" << std::endl;
    std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;
    std::cout << "----- Analyzing channel: " << ch << std::endl;

    int                        npoints    = fData->GetNpoints();
    TString                    collimator = fData->GetCollimator();
    TString                    testBench  = fData->GetTestBench();
    std::vector<double>        positions  = fData->GetPositions();
    std::vector<SFPeakFinder*> peakFin;

    for (int i = 0; i < npoints; i++)
    {
        if (ch == 0)
            peakFin.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
        else if (ch == 1)
            peakFin.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
        else
        {
            std::cerr << "##### Error in SFEnergyRes::CalculateEnergyRes() for ch" << ch
                      << std::endl;
            std::cerr << "Incorrect channel number!" << std::endl;
            return false;
        }
    }

    TString       gname = Form("ER_s%i_ch%i", fSeriesNo, ch);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("energy resolution [%]");
    graph->SetTitle(gname);
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    SFResults* parameters;
    double     enRes       = 0;
    double     enResErr    = 0;
    double     enResAve    = 0;
    double     enResAveErr = 0;

    for (int i = 0; i < npoints; i++)
    {
        peakFin[i]->FindPeakFit();
        parameters = peakFin[i]->GetResults();
        enRes      = parameters->GetValue(SFResultTypeNum::kPeakSigma) /
                     parameters->GetValue(SFResultTypeNum::kPeakPosition);
        enResErr   = enRes * sqrt(pow(parameters->GetUncertainty(SFResultTypeNum::kPeakPosition), 2) /
                                  pow(parameters->GetValue(SFResultTypeNum::kPeakPosition), 2) +
                                  pow(parameters->GetUncertainty(SFResultTypeNum::kPeakSigma), 2) /
                                  pow(parameters->GetValue(SFResultTypeNum::kPeakSigma), 2));
        enRes    = enRes * 100;
        enResErr = enResErr * 100;
        graph->SetPoint(i, positions[i], enRes);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), enResErr);
        enResAve += enRes * (1. / pow(enResErr, 2));
        enResAveErr += (1. / pow(enResErr, 2));
    }

    enResAve    = enResAve / enResAveErr;
    enResAveErr = sqrt(1. / enResAveErr);

    std::cout << "Average energy resolution for channel " << ch << ": " << enResAve << " +/- "
              << enResAveErr << " % \n"
              << std::endl;

    if (ch == 0)
    {
        fEnergyResCh0Graph = graph;
        fResultsCh0->AddResult(SFResultTypeNum::kEnergyRes, enResAve, enResAveErr);
        fResultsCh0->AddObject(SFResultTypeObj::kEnergyResGraph, fEnergyResCh0Graph);
    }
    else if (ch == 1)
    {
        fEnergyResCh1Graph = graph;
        fResultsCh1->AddResult(SFResultTypeNum::kEnergyRes, enResAve, enResAveErr);
        fResultsCh1->AddObject(SFResultTypeObj::kEnergyResGraph, fEnergyResCh1Graph);
    }

    return true;
}
//------------------------------------------------------------------
/// Calculates energy resolution based on averaged charge spectra for
/// all measurements in analyzed measurement series. In this function
/// fEnergyResGraphAve graph is filled and values of fResults.fEnergyResAve
/// and fResults.fEnergyResAveErr are assigned.
bool SFEnergyRes::CalculateEnergyRes(void)
{

    std::cout << "\n----- Inside SFEnergyRes::CalculateEnergyRes()" << std::endl;
    std::cout << "----- Analyzing series: " << fSeriesNo << std::endl;

    int                        npoints    = fData->GetNpoints();
    TString                    collimator = fData->GetCollimator();
    TString                    testBench  = fData->GetTestBench();
    std::vector<double>        positions  = fData->GetPositions();
    std::vector<int>           measIDs    = fData->GetMeasurementsIDs();
    std::vector<SFPeakFinder*> peakFin;

    TString       gname = Form("ER_s%i_ave", fSeriesNo);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("energy resolution [%]");
    graph->SetTitle(gname);
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    SFResults* parameters;
    double     enRes, enResErr;
    double     enResAve, enResAveErr;

    for (int i = 0; i < npoints; i++)
    {
        peakFin.push_back(new SFPeakFinder(fSpectraAve[i], 0));
        peakFin[i]->FindPeakFit();
        parameters = peakFin[i]->GetResults();

        enRes = parameters->GetValue(SFResultTypeNum::kPeakSigma) /
                parameters->GetValue(SFResultTypeNum::kPeakPosition);
        enResErr = enRes * sqrt(pow(parameters->GetUncertainty(SFResultTypeNum::kPeakPosition), 2) /
                                pow(parameters->GetValue(SFResultTypeNum::kPeakPosition), 2) +
                                pow(parameters->GetUncertainty(SFResultTypeNum::kPeakSigma), 2) /
                                pow(parameters->GetValue(SFResultTypeNum::kPeakSigma), 2));
        enRes    = enRes * 100;
        enResErr = enResErr * 100;
        graph->SetPoint(i, positions[i], enRes);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), enResErr);
        enResAve += enRes * (1. / pow(enResErr, 2));
        enResAveErr += (1. / pow(enResErr, 2));
    }

    enResAve    = enResAve / enResAveErr;
    enResAveErr = sqrt(1. / enResAveErr);

    if (std::isnan(enResAve) || enResAve < 0 || enResAve > 100)
    {
        std::cout << "##### Warning in SFEnergyRes class!" << std::endl;
        std::cout << "Incorrect energy resolution value: " << enResAve << " +/- " << enResAveErr
                  << std::endl;
        enResAve    = 0;
        enResAveErr = 0;
    }

    fEnergyResAveGraph = graph;
    fResultsAve->AddValue(SFResultTypeNum::kEnergyRes, enResAve);
    fResultsAve->AddUncertainty(SFResultTypeNum::kEnergyRes, enResAveErr);
    fResultsAve->AddObject(SFResultTypeObj::kEnergyResGraph, fEnergyResAveGraph);

    std::cout << "Average energy resolution calculated from summed and "
                 "attenuation length corrected histograms is: "
              << fResultsAve->GetValue(SFResultTypeNum::kEnergyRes) << " +/- "
              << fResultsAve->GetUncertainty(SFResultTypeNum::kEnergyRes) << " % \n"
              << std::endl;

    return true;
}
//------------------------------------------------------------------
/// Returns vector containing charge spectra of requested channel.
/// \par ch - channel number
std::vector<TH1D*> SFEnergyRes::GetSpectra(int ch)
{

    if ((ch == 0 && fSpectraCh0.empty()) || (ch == 1 && fSpectraCh1.empty()))
    {
        std::cerr << "##### Error in SFEnergyRes::GetSpectra() fo ch" << ch << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    if (ch == 0)
        return fSpectraCh0;
    else if (ch == 1)
        return fSpectraCh1;
    else
    {
        std::cerr << "##### Error in SFEnergyRes::GetSpectra for ch" << ch << std::endl;
        std::cerr << "Incorrect channel number!" << std::endl;
        std::abort();
    }
}
//------------------------------------------------------------------
/// Returns vector containing averaged charge spectra.
std::vector<TH1D*> SFEnergyRes::GetSpectra(void)
{

    if (fSpectraAve.empty())
    {
        std::cerr << "##### Error in SFEnergyRes::GetSpectraSum()" << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    return fSpectraAve;
}
//------------------------------------------------------------------
/// Returns vector containing background subtracted averaged charge spectra.
std::vector<TH1D*> SFEnergyRes::GetPeaks(void)
{

    if (fPeaksAve.empty())
    {
        std::cerr << "##### Error in SFEnergyRes::GetPeaks()" << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    return fPeaksAve;
}
//------------------------------------------------------------------
/// Returns vector containing background subtracted charge spectra of requested
/// channel.
/// \param ch - channel number
std::vector<TH1D*> SFEnergyRes::GetPeaks(int ch)
{

    if ((ch == 0 && fPeaksCh0.empty()) || (ch == 1 && fPeaksCh1.empty()))
    {
        std::cerr << "##### Error in SFEnergyRes::GetPeaks() for ch" << ch << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    if (ch == 0)
        return fPeaksCh0;
    else if (ch == 1)
        return fPeaksCh1;
    else
    {
        std::cerr << "##### Error in SFEnergyRes::GetPeaks() for ch" << ch << std::endl;
        std::cerr << "Incorrect channel number!" << std::endl;
        std::abort();
    }
}
//------------------------------------------------------------------
std::vector<SFResults*> SFEnergyRes::GetResults(void)
{
    
    if (fResultsCh0 == nullptr ||
        fResultsCh1 == nullptr ||
        fResultsAve == nullptr)
    {
        std::cerr << "##### Error in SFEnergyRes::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(3);
    results[0] = fResultsCh0;
    results[1] = fResultsCh1;
    results[2] = fResultsAve;
    
    return results;
}
//------------------------------------------------------------------
/// Prints details of the SFEnergyResolution class object.
void SFEnergyRes::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFEnergyRes class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
