// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFStabilityMon.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFStabilityMon.hh"

ClassImp(SFStabilityMon);

//------------------------------------------------------------------
SFStabilityMon::SFStabilityMon(int seriesNo) : fSeriesNo(seriesNo),
                                               fData(nullptr),
                                               fCh0PeakPosGraph(nullptr),
                                               fCh1PeakPosGraph(nullptr),
                                               fCh0ResidualGraph(nullptr),
                                               fCh1ResidualGraph(nullptr),
                                               fResultsCh0(nullptr),
                                               fResultsCh1(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFStabilityMon constructor!";
    }

    TString desc = fData->GetDescription();

    if (!desc.Contains("Stability monitoring"))
    {
        std::cerr << "##### Error in SFStabilityMon constructor!" << std::endl;
        std::cerr << "Series " << fSeriesNo << " is NOT stability monitoring series!" << std::endl;
        throw "##### Exception in SFStabilityMon constructor!";
    }

    double              s      = SFTools::GetSigmaBL(fData->GetSiPM());
    std::vector<double> sigma  = {s};
    TString             cutCh0 = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
    TString             cutCh1 = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);

    fSpecCh0 = fData->GetSpectra(0, SFSelectionType::kPE, cutCh0);
    fSpecCh1 = fData->GetSpectra(1, SFSelectionType::kPE, cutCh1);

    fResultsCh0 = new SFResults(Form("StabilityResults_S%i_ch0", fSeriesNo));
    fResultsCh1 = new SFResults(Form("StabilityResults_S%i_ch1", fSeriesNo));
}
//------------------------------------------------------------------
SFStabilityMon::~SFStabilityMon()
{
    delete fData;
}
//------------------------------------------------------------------
bool SFStabilityMon::AnalyzeStability(int ch)
{
    std::cout << "\n\n----- Stability analysis " << std::endl;
    std::cout << "----- Series: " << fSeriesNo << std::endl;
    std::cout << "----- Channel: " << ch << std::endl;

    int npoints = fData->GetNpoints();

    std::vector<SFPeakFinder*> peakFin;
    std::vector<TH1D*>         spec;
    SFResults*                 peakParams;
    std::vector<double>        peakPositions;

    if (ch == 0)
        spec = fSpecCh0;
    else if (ch == 1)
        spec = fSpecCh1;
    else
    {
        std::cerr << "##### Error in SFStabilityMon::AnalyzeStability()!" << std::endl;
        std::cerr << "Incorrect channel number! Please check!" << std::endl;
        return false;
    }

    TGraphErrors* gPeakPos = new TGraphErrors(npoints);
    gPeakPos->SetName(Form("511PeakPosition_S%i_ch%i", fSeriesNo, ch));
    gPeakPos->SetTitle(Form("511PeakPosition_S%i_ch%i", fSeriesNo, ch));
    gPeakPos->GetXaxis()->SetTitle("source position [mm]");
    gPeakPos->GetYaxis()->SetTitle("511 keV peak position [PE]");
    gPeakPos->SetMarkerStyle(4);

    TGraphErrors* gResiduals = new TGraphErrors(npoints);
    gResiduals->SetName(Form("Residuals_S%i_ch%i", fSeriesNo, ch));
    gResiduals->SetTitle(Form("Residuals_S%i_ch%i", fSeriesNo, ch));
    gResiduals->GetXaxis()->SetTitle("source position [mm]");
    gResiduals->GetYaxis()->SetTitle("residual [PE]");
    gResiduals->SetMarkerStyle(4);

    for (int i = 0; i < npoints; i++)
    {
        peakFin.push_back(new SFPeakFinder(spec[i], false));
        peakFin[i]->FindPeakFit();
        peakParams = peakFin[i]->GetResults();
        gPeakPos->SetPoint(i, i, peakParams->GetValue(SFResultTypeNum::kPeakPosition));
        gPeakPos->SetPointError(i, 0, peakParams->GetUncertainty(SFResultTypeNum::kPeakPosition));
        peakPositions.push_back(peakParams->GetValue(SFResultTypeNum::kPeakPosition));
    }

    double mean   = SFTools::GetMean(peakPositions);
    double stdDev = SFTools::GetStandardDev(peakPositions);

    TF1* funPol0 = new TF1("funPol0", "pol0", 0, npoints);
    funPol0->FixParameter(0, mean);
    gPeakPos->Fit(funPol0, "Q");

    std::cout << "\tCalculation results: " << mean << "+/-" << stdDev << std::endl;

    double x, y;
    double res;

    for (int i = 0; i < npoints; i++)
    {
        gPeakPos->GetPoint(i, x, y);
        res = y - mean;
        gResiduals->SetPoint(i, i, res);
    }

    if (ch == 0)
    {
        fCh0PeakPosGraph  = gPeakPos;
        fCh0ResidualGraph = gResiduals;
        fResultsCh0->AddResult(SFResultTypeNum::kAveragePeakPos, mean, stdDev);
        fResultsCh0->AddObject(SFResultTypeObj::kPeakPosGraph, fCh0PeakPosGraph);
        fResultsCh0->AddObject(SFResultTypeObj::kSMResidualGraph, fCh0ResidualGraph);
    }
    else if (ch == 1)
    {
        fCh1PeakPosGraph  = gPeakPos;
        fCh1ResidualGraph = gResiduals;
        fResultsCh1->AddResult(SFResultTypeNum::kAveragePeakPos, mean, stdDev);
        fResultsCh1->AddObject(SFResultTypeObj::kPeakPosGraph, fCh1PeakPosGraph);
        fResultsCh1->AddObject(SFResultTypeObj::kSMResidualGraph, fCh1ResidualGraph);
    }

    return true;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFStabilityMon::GetSpectra(int ch)
{
    std::vector<TH1D*> tmp;

    if (ch == 0)
        tmp = fSpecCh0;
    else if (ch == 1)
        tmp = fSpecCh1;
    else
    {
        std::cerr << "##### Error in SFStabilityMon::GetSpectra()! " << std::endl;
        std::cerr << "Incorrect channel number! Please check!" << std::endl;
        std::abort();
    }

    if (tmp.empty())
    {
        std::cerr << "##### Error in SFStabilityMon::GetSpectra()! " << std::endl;
        std::cerr << "Requested spectra don't exist!" << std::endl;
        std::abort();
    }

    return tmp;
}
//-----------------------------------------------------------------
std::vector<SFResults*> SFStabilityMon::GetResults(void)
{
    if(fResultsCh0 == nullptr ||
       fResultsCh1 == nullptr)
    {
        std::cerr << "##### Error in SFStabilityMon::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(2);
    results[0] = fResultsCh0;
    results[1] = fResultsCh1;
    
    return results;
}
//------------------------------------------------------------------
void SFStabilityMon::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFStabilityMon class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
