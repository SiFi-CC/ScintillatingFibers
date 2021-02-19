// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFAttenuation.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFAttenuation.hh"

ClassImp(SFAttenuation);

//------------------------------------------------------------------
/// Standard constructor (recommended)
/// \param seriesNo is number of experimental series to be analyzed.
SFAttenuation::SFAttenuation(int seriesNo) : fSeriesNo(seriesNo),
                                             fData(nullptr),
                                             fAttGraph(nullptr),
                                             fSigmaGraph(nullptr),
                                             fAttCh0Graph(nullptr),
                                             fAttCh1Graph(nullptr),
                                             fResultsCh0(nullptr),
                                             fResultsCh1(nullptr),
                                             fResultsCombPol1(nullptr),
                                             fResultsCombPol3(nullptr)
{
    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFAttenuation constructor!";
    }

    TString desc = fData->GetDescription();

    if (!desc.Contains("Regular series"))
    {
        std::cout << "##### Error in SFAttenuation constructor! Non-regular series!" << std::endl;
        throw "##### Exception in SFAttenuation constructor!";
    }

    fResultsCh0 = new SFResults(Form("AttenuationResults_S%i_ch0", fSeriesNo));
    fResultsCh1 = new SFResults(Form("AttenuationResults_S%i_ch1", fSeriesNo));
    fResultsCombPol1 = new SFResults(Form("AttenuationResults_S%i_CombPol1", fSeriesNo));
    fResultsCombPol3 = new SFResults(Form("AttenuationResults_S%i_CombPol3", fSeriesNo));
}
//------------------------------------------------------------------
/// Default destructor.
SFAttenuation::~SFAttenuation()
{
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
/// Method to determine attenuation length used in Pauwels et al., JINST 8 (2013) P09019.
/// For both ends of the fiber one value is calculated, since combined signal from both channels
/// is taken into account.
bool SFAttenuation::AttCombinedCh(void)
{
    std::cout << "\n----- Inside SFAttenuation::AttCombinedCh() for series " << fSeriesNo
              << std::endl;

    int                 npoints    = fData->GetNpoints();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             sipm       = fData->GetSiPM();
    std::vector<double> positions  = fData->GetPositions();

    double              s      = SFTools::GetSigmaBL(sipm);
    std::vector<double> sigmas = {s, s};
    TString             cut    = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);
    fRatios                    = fData->GetCustomHistograms(SFSelectionType::kLogSqrtPERatio, cut);

    std::vector<TF1*> fun;

    TString gname = Form("att_s%i", fSeriesNo);
    fAttGraph     = new TGraphErrors(npoints);
    fAttGraph->GetXaxis()->SetTitle("source position [mm]");
    fAttGraph->GetYaxis()->SetTitle("ln(M_{LR})");
    fAttGraph->SetTitle(gname);
    fAttGraph->SetName(gname);
    fAttGraph->SetMarkerStyle(4);

    gname       = Form("sigmas_s%i", fSeriesNo);
    fSigmaGraph = new TGraphErrors(npoints);
    fSigmaGraph->GetXaxis()->SetTitle("source position [mm]");
    fSigmaGraph->GetYaxis()->SetTitle("#sigma M_{LR}");
    fSigmaGraph->SetTitle(gname);
    fSigmaGraph->SetName(gname);
    fSigmaGraph->SetMarkerStyle(4);

    int     parNo    = 1;
    double  const_1  = 0;
    double  const_2  = 0;
    TString fun_name = "";

    for (int i = 0; i < npoints; i++)
    {
        if (collimator.Contains("Lead"))
        {
            SFTools::RatiosFitDoubleGauss(fRatios, 5);
            fun_name = "fDGauss";
            const_1  = fRatios[i]->GetFunction(fun_name)->GetParameter(0);
            const_2  = fRatios[i]->GetFunction(fun_name)->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("SensL"))
        {
            SFTools::RatiosFitDoubleGauss(fRatios, 2);
            fun_name = "fDGauss";
            const_1  = fRatios[i]->GetFunction(fun_name)->GetParameter(0);
            const_2  = fRatios[i]->GetFunction(fun_name)->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("Hamamatsu"))
        {
            SFTools::RatiosFitGauss(fRatios, 1);
            fun_name = "fGauss";
            parNo    = 1;
        }

        fAttGraph->SetPoint(i, positions[i],
                            fRatios[i]->GetFunction(fun_name)->GetParameter(parNo));
        fAttGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                                 fRatios[i]->GetFunction(fun_name)->GetParError(parNo));
        fSigmaGraph->SetPoint(i, positions[i],
                              fRatios[i]->GetFunction(fun_name)->GetParameter(parNo + 1));
        fSigmaGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                                   fRatios[i]->GetFunction(fun_name)->GetParError(parNo + 1));
    }

    Fit1stOrder();
    Fit3rdOrder();

    fResultsCombPol1->AddObject(SFResultTypeObj::kAttGraph, fAttGraph);
    fResultsCombPol1->AddObject(SFResultTypeObj::kMLRSigmaGraph, fSigmaGraph);

    return true;
}
//------------------------------------------------------------------
double myPol3(double* x, double* par)
{
    // double xx = (x[0]-par[3]);
    // double f = par[0] + par[1]*(xx + par[2]*pow(xx,3));
    double f = par[0] + par[1] * (x[0] + par[2] * pow(x[0], 3));
    return f;
}
//------------------------------------------------------------------
bool SFAttenuation::Fit1stOrder(void)
{
    if (fAttGraph == nullptr) AttCombinedCh();

    TF1* fpol1 = new TF1("fpol1", "pol1", -50, 150);
    fpol1->SetParameters(-0.15, 0.005);
    TFitResultPtr ptr = fAttGraph->Fit(fpol1, "SQR+");

    double att     = fabs(1. / fpol1->GetParameter(1));
    double err     = fpol1->GetParError(1) / pow(fpol1->GetParameter(1), 2);
    double Chi2NDF = ptr->Chi2() / ptr->Ndf();

    fResultsCombPol1->AddResult(SFResultTypeNum::kLambda, att, err);
    fResultsCombPol1->AddResult(SFResultTypeNum::kChi2NDF, Chi2NDF, -1);

    std::cout << "Attenuation lenght from pol1 fit is: " << att << " +/- " << err << " mm"
              << std::endl;
    std::cout << "Chi2/NDF = " << Chi2NDF << std::endl << std::endl;
}
//------------------------------------------------------------------
bool SFAttenuation::Fit3rdOrder(void)
{
    if (fAttGraph == nullptr) AttCombinedCh();

    // double fiberLengthHalf = fData->GetFiberLength()/2;

    TF1* fpol3 = new TF1("fpol3", myPol3, -50, 150, 3);
    fpol3->SetParameter(0, fAttGraph->GetFunction("fpol1")->GetParameter(0));
    fpol3->SetParameter(1, fAttGraph->GetFunction("fpol1")->GetParameter(1));
    // fpol3->FixParameter(3, fiberLengthHalf);
    // fpol3->SetParameter(3, fiberLengthHalf);
    // fpol3->SetParLimits(0, 0, 1);

    fpol3->SetLineColor(kBlue - 7);

    TFitResultPtr ptr = fAttGraph->Fit(fpol3, "SQR+");

    double att     = fabs(1. / fpol3->GetParameter(1));
    double err     = fpol3->GetParError(1) / pow(fpol3->GetParameter(1), 2);
    double Chi2NDF = ptr->Chi2() / ptr->Ndf();

    fResultsCombPol3->AddResult(SFResultTypeNum::kLambda, att, err);
    fResultsCombPol3->AddResult(SFResultTypeNum::kChi2NDF, Chi2NDF, -1);

    std::cout << "Attenuation lenght from pol3 fit is: " << att << " +/- " << err << " mm"
              << std::endl;
    std::cout << "Chi2/NDF = " << Chi2NDF << std::endl << std::endl;

    return true;
}
//------------------------------------------------------------------
/// Method to determine attenuation length for both channels independently.
/// If series was measured with lead collimator peak position is determied
/// with the FindPeakNoBackground() method of the SFPeakFinder class. If
/// series was measured with electronic collimator - FindPeakFit() method
/// of the SFPeakFinder class is used.
/// \param ch - channel number
bool SFAttenuation::AttSeparateCh(int ch)
{
    std::cout << "\n----- Inside SFAttenuation::AttSeparateCh() for series " << fSeriesNo
              << std::endl;
    std::cout << "----- Analyzing channel " << ch << std::endl;

    int                 npoints    = fData->GetNpoints();
    double              fiberLen   = fData->GetFiberLength();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             sipm       = fData->GetSiPM();
    std::vector<double> positions  = fData->GetPositions();

    TString             cut   = "";
    double              s     = SFTools::GetSigmaBL(sipm);
    std::vector<double> sigma = {s};

    if (ch == 0)
        cut = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
    else if (ch == 1)
        cut = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);

    std::vector<TH1D*> spectra = fData->GetSpectra(ch, SFSelectionType::kPE, cut);

    TString       gname = Form("att_s%i_ch%i", fSeriesNo, ch);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    graph->SetTitle(gname);
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    std::vector<SFPeakFinder*> peakfin;
    SFResults*                 peakParams;

    // TString fname = Form("/home/kasia/S%ich%i.txt", fSeriesNo, ch);
    // std::ofstream output(fname);

    for (int i = 0; i < npoints; i++)
    {
        peakfin.push_back(new SFPeakFinder(spectra[i], false));
        peakfin[i]->FindPeakFit();
        peakParams = peakfin[i]->GetResults();
        graph->SetPoint(i, positions[i], peakParams->GetValue(SFResultTypeNum::kPeakPosition));
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                             peakParams->GetUncertainty(SFResultTypeNum::kPeakPosition));
        // output << positions[i] << "\t" << fResults->GetValue(SFResultType::fPeakPosition)
        //       << "\t" << fResults->GetUncertainty(SFResultType::fPeakPosition) << "\t"
        //       << fResults->GetValue(SFResultType::fPeakSigma) << std::endl;
    }

    //----- fitting
    TF1* fexp;

    if (ch == 0)
    {
        fexp = new TF1("funCh0", "[0]*exp(-x/[1])", positions[0], positions[npoints - 1]);
        fexp->SetParameters(500, 200);
    }
    else if (ch == 1)
    {
        fexp = new TF1("funCh1", "[0]*exp(-([2]-x)/[1])", positions[0], positions[npoints - 1]);
        fexp->SetParameter(0, 500);
        fexp->SetParameter(1, 200);
        fexp->FixParameter(2, fiberLen);
    }

    TFitResultPtr ptr     = graph->Fit(fexp, "QRS");
    double        Chi2NDF = ptr->Chi2() / ptr->Ndf();

    //----- calculating attenuation length
    std::cout << "\n\tAttenuation for channel " << ch << ": " << fexp->GetParameter(1) << " +/- "
              << fexp->GetParError(1) << " mm" << std::endl;
    std::cout << "Chi2/NDF = " << Chi2NDF << std::endl << std::endl;

    if (ch == 0)
    {
        fResultsCh0->AddResult(SFResultTypeNum::kLambda, fexp->GetParameter(1), fexp->GetParError(1));
        fResultsCh0->AddResult(SFResultTypeNum::kChi2NDF, Chi2NDF, -1);
        fSpectraCh0  = spectra;
        fAttCh0Graph = graph;
        fResultsCh1->AddObject(SFResultTypeObj::kAttGraph, fAttCh0Graph);
    }
    else if (ch == 1)
    {
        fResultsCh1->AddResult(SFResultTypeNum::kLambda, fexp->GetParameter(1), fexp->GetParError(1));
        fResultsCh1->AddResult(SFResultTypeNum::kChi2NDF, Chi2NDF, -1);
        fSpectraCh1  = spectra;
        fAttCh1Graph = graph;
        fResultsCh1->AddObject(SFResultTypeObj::kAttGraph, fAttCh1Graph);
    }

    return true;
}
//------------------------------------------------------------------
/// Returns vector containing PE spectra used in determination of attenuation
/// length with separate channels method i.e. AttSeparateCh().
///\param ch - channel number
std::vector<TH1D*> SFAttenuation::GetSpectra(int ch)
{
    if ((ch == 0 && fSpectraCh0.empty()) || (ch == 1 && fSpectraCh1.empty()))
    {
        std::cerr << "##### Error in SFAttenuation::GetSpectra(). Empty vector!" << std::endl;
        std::abort();
    }
    if (ch == 0)
        return fSpectraCh0;
    else if (ch == 1)
        return fSpectraCh1;
}
//------------------------------------------------------------------
/// Returns vector containing peaks (spectra after background subtraction with
/// SFPeakFinder class) used in determination of attenuation length with separate
/// channels method i.e. AttSeparateCh().
///\param ch - channel number.
std::vector<TH1D*> SFAttenuation::GetPeaks(int ch)
{
    if ((ch == 0 && fPeaksCh0.empty()) || (ch == 1 && fPeaksCh1.empty()))
    {
        std::cerr << "##### Error in SFAttenuation::GetPeaks(). Empty vector!" << std::endl;
        std::abort();
    }
    if (ch == 0)
        return fPeaksCh0;
    else if (ch == 1)
        return fPeaksCh1;
}
//------------------------------------------------------------------
/// Returns vector containing histograms with signal ratios from both channels.
/// Histograms are used during averaged channels analysis i.e. in AttAveragedCh().
std::vector<TH1D*> SFAttenuation::GetRatios(void)
{
    if (fRatios.empty())
    {
        std::cerr << "#### Error in SFAttenuation::GetRatios(). Empty vector!" << std::endl;
        std::abort();
    }
    return fRatios;
}
//------------------------------------------------------------------
std::vector<SFResults*> SFAttenuation::GetResults(void)
{
    if (fResultsCh0 == nullptr ||
        fResultsCh1 == nullptr ||
        fResultsCombPol1 == nullptr ||
        fResultsCombPol3 == nullptr)
    {
        std::cerr << "##### Error in SFAttenuation::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(4);
    results[0] = fResultsCh0;
    results[1] = fResultsCh1;
    results[2] = fResultsCombPol1;
    results[3] = fResultsCombPol3;
    
    return results;    
}
//------------------------------------------------------------------
/// Prints details of SFAttenuation class object.
void SFAttenuation::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFAttenuation class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
