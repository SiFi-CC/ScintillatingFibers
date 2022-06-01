// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimingRes.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFTimingRes.hh"

ClassImp(SFTimingRes);

//------------------------------------------------------------------
/// Standard constructor
/// \param seriesNo is number of experimental series to be analyzed.
SFTimingRes::SFTimingRes(int seriesNo) : fSeriesNo(seriesNo),
                                         fData(nullptr),
                                         fTResGraph(nullptr),
                                         fTResECutGraph(nullptr),
                                         fTSigmaGraph(nullptr),
                                         fTSigmaECutGraph(nullptr),
                                         fResults(nullptr),
                                         fResultsECut(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFTimingRes constructor!";
    }
    
    fResults     = new SFResults(Form("TimingResResults_S%i", fSeriesNo));
    fResultsECut = new SFResults(Form("TimingResResults_S%i_ECut", fSeriesNo));
}
//------------------------------------------------------------------
/// Default destructor.
SFTimingRes::~SFTimingRes()
{
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
/// Private method to get ratio histograms necesarry to impose cuts.
/// Depending on measurement type single or double gaussian function
/// is fitted to the histograms.
bool SFTimingRes::LoadRatios(void)
{
    int     npoints    = fData->GetNpoints();
    TString collimator = fData->GetCollimator();
    TString sipm       = fData->GetSiPM();
    TString testBench  = fData->GetTestBench();

    TString cut = " ";
    
    if (testBench == "PMI")
    {
        cut = SFDrawCommands::GetCut(SFCutType::kPMICombCh0Ch1);
        fRatios = fData->GetCustomHistograms(SFSelectionType::kPMILogSqrtChargeRatio, cut);
    }
    else
    {
        double s = SFTools::GetSigmaBL(sipm);
        std::vector<double> sigmas = {s, s};
        cut = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);
        fRatios = fData->GetCustomHistograms(SFSelectionType::kLogSqrtPERatio, cut);
    }

    std::vector<TF1*> fun;

    for (int i = 0; i < npoints; i++)
    {
        if (collimator.Contains("Lead")) { SFTools::RatiosFitDoubleGauss(fRatios, 5); }
        else if (collimator.Contains("Electronic") && sipm.Contains("SensL"))
        {
            SFTools::RatiosFitDoubleGauss(fRatios, 2);
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("Hamamatsu"))
        {
            SFTools::RatiosFitGauss(fRatios, 1);
        }
    }

    return true;
}
//------------------------------------------------------------------
/// Formula for Lorentzian function (if needed).
double LorentzianFun(double* x, double* par)
{
    // par0 - integral
    // par1 - width (FWHM)
    // par2 - mean
    return (0.5 * par[0] * par[1] / TMath::Pi()) /
           ((x[0] - par[2]) * (x[0] - par[2]) + 0.25 * par[1] * par[1]);
}
//------------------------------------------------------------------
/// Method for simple timing resolution determination. Timing resolution spectrum
/// is build as difference between channel 0 and 1 signals with additional cut
/// on scattered events. Timing resolution and its uncertainty is determined
/// based on Gaussian function as its sigma.
/// If measurement was taken with lead collimator - double Gauss fitted.
/// If measurement was taken with electronic collimator - single Gaussian.
bool SFTimingRes::AnalyzeNoECut(void)
{

    std::cout << "----- Inside SFTimingRes::AnalyzeNoECut()" << std::endl;
    std::cout << "----- Series number: " << fSeriesNo << std::endl;

    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    std::vector<int>    measIDs    = fData->GetMeasurementsIDs();
    TString             collimator = fData->GetCollimator();
    TString             sipm       = fData->GetSiPM();
    TString             testBench  = fData->GetTestBench();
    TString             fiber      = fData->GetFiber();

    TString cut;
    double  mean, sigma;

    if (fRatios.empty()) LoadRatios();
    std::vector<TF1*> fun;

    TString gname = Form("timeDiff_S%i", fSeriesNo);
    fTResGraph    = new TGraphErrors(npoints);
    fTResGraph->GetXaxis()->SetTitle("source position [mm]");
    fTResGraph->GetYaxis()->SetTitle("#mu_{T_{D}} [ns]");
    fTResGraph->SetTitle(Form("#mu of T_{D} Distribution S%i", fSeriesNo));
    fTResGraph->SetName(gname);
    fTResGraph->SetMarkerStyle(4);
    
    gname        = Form("timeSigma_S%i", fSeriesNo);
    fTSigmaGraph = new TGraphErrors(npoints);
    fTSigmaGraph->GetYaxis()->SetTitle("#sigma_{T_D} [ns]");
    fTSigmaGraph->GetXaxis()->SetTitle("source position [mm]");
    fTSigmaGraph->SetName(gname);
    fTSigmaGraph->SetTitle(Form("#sigma of T_{D} Distribution S%i", fSeriesNo));
    fTSigmaGraph->SetMarkerStyle(4);

    double s = SFTools::GetSigmaBL(sipm);
    
    double tResSum = 0.;
    double tResSumErr = 0.;
    
    for (int i = 0; i < npoints; i++)
    {

        int parNum = 0;

        if (collimator.Contains("Lead"))
        {

            if (fRatios[i]->GetFunction("fDGauss")->GetParameter(0) >
                fRatios[i]->GetFunction("fDGauss")->GetParameter(3))
                parNum = 1;
            else
                parNum = 4;

            mean  = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum);
            sigma = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum + 1);

            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 &&
            // ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",
            // mean-0.5*sigma, mean+0.5*sigma); cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590
            // && ch_1.fT0<590 && ch_0.fPE>15 && ch_1.fPE>15 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-0.5*sigma, mean+0.5*sigma);
            
            std::vector<double> customNum = {s, s, 15, 15, 
                                             /*mean - 0.5 * sigma, // lead collimator + Hamamatsu SiPMs, 1x1 mm fibers (Jan 2018)
                                             mean + 0.5 * sigma*/
                                             mean - 2*sigma, // lead collimator + Hamamatsu SiPMs, 2x3 mm fiber (March 2022)
                                             mean + 2*sigma};
            cut = SFDrawCommands::GetCut(SFCutType::kT0Diff, customNum);
            fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::kT0Difference,
                                                        cut, measIDs[i]));

            fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -100, 100));
            fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
            fun[i]->SetParLimits(0, 1, 1E6);
            fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
            fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
            fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()) / 4);
            fun[i]->SetParLimits(3, 1, 1E6);
            if (i < npoints / 2)
                fun[i]->SetParameter(4, fT0Diff[i]->GetMean() - 5);
            else
                fun[i]->SetParameter(4, fT0Diff[i]->GetMean() + 5);
            fun[i]->SetParameter(5, fT0Diff[i]->GetRMS() * 2);
            fT0Diff[i]->Fit(fun[i], "QR");
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("SensL"))
        {

            if (fRatios[i]->GetFunction("fDGauss")->GetParameter(0) >
                fRatios[i]->GetFunction("fDGauss")->GetParameter(3))
                parNum = 1;
            else
                parNum = 4;

            mean  = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum);
            sigma = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNum + 1);

            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 &&
            // ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",
            // mean-2*sigma,  mean+2*sigma); cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 &&
            // ch_1.fT0<590 && ch_0.fPE>20 && ch_1.fPE>20 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f", mean-2*sigma,  mean+2*sigma);
            
            std::vector<double> customNum = {s, s, 20, 20, 
                                             mean - 2 * sigma,
                                             mean + 2 * sigma};
            cut = SFDrawCommands::GetCut(SFCutType::kT0Diff, customNum);
            fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::kT0Difference,
                                                        cut, measIDs[i]));
            fT0Diff.back()->Rebin(2);

            fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -30, 30));
            fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
            fun[i]->SetParLimits(0, 1, 1E6);
            fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
            fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
            fun[i]->SetParLimits(2, 0, 10);
            fun[i]->SetParameter(3, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()) / 10.);
            fun[i]->SetParLimits(3, 1, 1E6);
            fun[i]->SetParameter(4, fT0Diff[i]->GetMean());
            fun[i]->SetParameter(5, fT0Diff[i]->GetRMS() * 10);
            fun[i]->SetParLimits(5, 0, 20);
            fT0Diff[i]->Fit(fun[i], "R");
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("Hamamatsu"))
        {

            mean  = fRatios[i]->GetFunction("fGauss")->GetParameter(1);
            sigma = fRatios[i]->GetFunction("fGauss")->GetParameter(2);
            
            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>0 &&
            // ch_1.fPE>0 && log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",
            // mean-3*sigma,  mean+3*sigma);
            
            std::vector<double> customNum = {s, s, 0, 0,
                                             mean - 3 * sigma,
                                             mean + 3 * sigma};
            cut = SFDrawCommands::GetCut(SFCutType::kT0Diff, customNum);
            fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::kT0Difference,
                                                        cut, measIDs[i]));
            fT0Diff.back()->Rebin(2);

            fun.push_back(new TF1("fun", "gaus(0)+gaus(3)", -50, 50));
            fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
            fun[i]->SetParLimits(0, 1, 1E6);
            fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
            fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
            fun[i]->SetParLimits(2, 0, 50);
            if (fiber.Contains("LuAG"))
            {
                fun[i]->FixParameter(3, 0);
                fun[i]->FixParameter(4, 0);
                fun[i]->FixParameter(5, 0);
            }
            else if (fiber.Contains("LYSO"))
            {
                fun[i]->SetParameter(3,
                                     fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()) / 10);
                fun[i]->SetParLimits(3, 1, 1E6);
                if (i < npoints / 2)
                    fun[i]->SetParameter(4, fT0Diff[i]->GetMean() - 2);
                else
                    fun[i]->SetParameter(4, fT0Diff[i]->GetMean() + 2);
                fun[i]->SetParameter(5, fT0Diff[i]->GetRMS() * 5);
                fun[i]->SetParLimits(5, 0, 50);
            }
            else if (fiber.Contains("GAGG"))
            {
                fun[i]->SetParameter(3,
                                     fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()) / 10);
                fun[i]->SetParLimits(3, 1, 1E6);
                if (i < npoints / 2)
                    fun[i]->SetParameter(4, fT0Diff[i]->GetMean() - 1);
                else
                    fun[i]->SetParameter(4, fT0Diff[i]->GetMean() + 1);
                fun[i]->SetParameter(5, fT0Diff[i]->GetRMS() * 2);
                fun[i]->SetParLimits(5, 0, 50);
            }
            fT0Diff[i]->Fit(fun[i], "QR");
        }
        else if (testBench == "PMI")
        {
            mean = fRatios[i]->GetFunction("fGauss")->GetParameter(1);
            sigma = fRatios[i]->GetFunction("fGauss")->GetParameter(2);
            
            std::vector<double> customNum = {mean - 2 * sigma,  //TODO optimize
                                             mean + 2 * sigma};
            cut = SFDrawCommands::GetCut(SFCutType::kPMIT0Diff, customNum);
            fT0Diff.push_back(fData->GetCustomHistogram(SFSelectionType::kPMIT0Difference,
                                                        cut, measIDs[i]));
            
            fun.push_back(new TF1("fun", "gaus", -50, 50)); //TODO optimize function
            fun[i]->SetParameter(0, fT0Diff[i]->GetBinContent(fT0Diff[i]->GetMaximumBin()));
            fun[i]->SetParameter(1, fT0Diff[i]->GetMean());
            fun[i]->SetParameter(2, fT0Diff[i]->GetRMS());
            fT0Diff[i]->Fit(fun[i], "QR");
        }

        parNum = 0;

        if (fun[i]->GetNpar() > 3)
        {
            if (fun[i]->GetParameter(0) > fun[i]->GetParameter(3))
                parNum = 2;
            else
                parNum = 5;
        }
        else
            parNum = 2;
        
        fTSigmaGraph->SetPoint(i, positions[i], fun[i]->GetParameter(parNum));
        fTSigmaGraph->SetPointError(i ,SFTools::GetPosError(collimator, testBench),
                                    fun[i]->GetParError(parNum));

        fTResGraph->SetPoint(i, positions[i], fun[i]->GetParameter(parNum - 1));
        fTResGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                             fun[i]->GetParameter(parNum));

        tResSum += fun[i]->GetParameter(parNum) * (1. / pow(fun[i]->GetParError(parNum), 2));
        tResSumErr += 1. / pow(fun[i]->GetParError(parNum), 2);
    }

    double tResAv    = tResSum / tResSumErr;
    double tResAvErr = sqrt(1. / tResSumErr); 
    
    fResults->AddResult(SFResultTypeNum::kTimeRes, tResAv, tResAvErr);
    fResults->AddObject(SFResultTypeObj::kTimeResGraph, fTResGraph);
    fResults->AddObject(SFResultTypeObj::kTimeSigGraph, fTSigmaGraph);

    std::cout << "Average timing resolution: " << tResAv << " +/- "
              << tResAvErr << " ns" << std::endl;

    return true;
}
//------------------------------------------------------------------
/// Method for timing resolution determination with an energy cut on 511 keV peak.
/// Timing resolution based on the Gaussian fit to the histogram - sigma of the fitted
/// function. Regardless the measurement type - always sigle Gauss fitted.
bool SFTimingRes::AnalyzeWithECut(void)
{

    std::cout << "----- Inside SFTimingRes::AnalyzeWithECut()" << std::endl;
    std::cout << "----- Series number " << fSeriesNo << std::endl;

    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    std::vector<int>    measIDs    = fData->GetMeasurementsIDs();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             sipm       = fData->GetSiPM();

    if (fRatios.empty()) LoadRatios();
    
    double s = SFTools::GetSigmaBL(sipm);
    std::vector<double> cut_sigma = {s};
     
    if (testBench == "PMI")
    {
        TString cutCh0    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh0, cut_sigma);
        TString cutCh1    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh1, cut_sigma);
        fSpecCh0  = fData->GetSpectra(0, SFSelectionType::kPMICharge, cutCh0);
        fSpecCh1  = fData->GetSpectra(1, SFSelectionType::kPMICharge, cutCh1);
    }
    else
    {
        TString cutCh0    = SFDrawCommands::GetCut(SFCutType::kSpecCh0, cut_sigma);
        TString cutCh1    = SFDrawCommands::GetCut(SFCutType::kSpecCh1, cut_sigma);
        fSpecCh0  = fData->GetSpectra(0, SFSelectionType::kPE, cutCh0);
        fSpecCh1  = fData->GetSpectra(1, SFSelectionType::kPE, cutCh1);
    }
    
    std::vector<SFPeakFinder*> peakFin_ch0;
    std::vector<SFPeakFinder*> peakFin_ch1;
    TF1*                       fun = new TF1("fun", "gaus", -200, 200);
    double                     xmin_ch0, xmax_ch0;
    double                     xmin_ch1, xmax_ch1;
    double                     mean_ratio, sigma_ratio;
    double                     mean, sigma;
    double                     const_1, const_2;
    // double f = 2*sqrt(2*log(2));   //to recalculate sigma into FWHM
    TString cut;

    double center_ch0, delta_ch0; // change here for smaller cut
    double center_ch1, delta_ch1; //

    double tResAv    = 0;
    double tResAvErr = 0;
    
    double tResSum    = 0;
    double tResSumErr = 0;

    int parNo = 0;

    TString  gname = Form("timeDiffECut_S%i", fSeriesNo);
    fTResECutGraph = new TGraphErrors(npoints);
    fTResECutGraph->GetXaxis()->SetTitle("source position [mm]");
    fTResECutGraph->GetYaxis()->SetTitle("#mu_{T_{D}} [ns]");
    fTResECutGraph->SetTitle(Form("#mu of T_{D} Distribution (with Energy Cut) S%i", fSeriesNo));
    fTResECutGraph->SetName(gname);
    fTResECutGraph->SetMarkerStyle(4);

    gname = Form("timeSigmaECut_S%i", fSeriesNo);
    fTSigmaECutGraph = new TGraphErrors(npoints);
    fTSigmaECutGraph->GetXaxis()->SetTitle("source position [mm]");
    fTSigmaECutGraph->GetYaxis()->SetTitle("#sigma_{T_{D}} [ns]");
    fTSigmaECutGraph->SetTitle(Form("#sigma of T_{D} Distribution (with Energy Cut) S%i", fSeriesNo));
    fTSigmaECutGraph->SetName(gname);
    fTSigmaECutGraph->SetMarkerStyle(4);
    
    for (int i = 0; i < npoints; i++)
    {
        peakFin_ch0.push_back(new SFPeakFinder(fSpecCh0[i], false));
        peakFin_ch1.push_back(new SFPeakFinder(fSpecCh1[i], false));
        peakFin_ch0[i]->FindPeakRange(xmin_ch0, xmax_ch0);
        peakFin_ch1[i]->FindPeakRange(xmin_ch1, xmax_ch1);

        center_ch0 = xmin_ch0 + (xmax_ch0 - xmin_ch0) / 2.; // change here for smaller cut
        delta_ch0  = (xmax_ch0 - xmin_ch0) / 6.;            //
        center_ch1 = xmin_ch1 + (xmax_ch1 - xmin_ch1) / 2.; //
        delta_ch1  = (xmax_ch1 - xmin_ch1) / 6.;            //

        if (collimator.Contains("Lead"))
        {
            const_1 = fRatios[i]->GetFunction("fDGauss")->GetParameter(0);
            const_2 = fRatios[i]->GetFunction("fDGauss")->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
            mean_ratio  = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo);
            sigma_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo + 1);
            
            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f
            // && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f", center_ch0-delta_ch0, center_ch0+delta_ch0,
            // center_ch1-delta_ch1, center_ch1+delta_ch1, mean_ratio-0.5*sigma_ratio,
            // mean_ratio+0.5*sigma_ratio);   //changed here for smaller cut
            
            std::vector<double> customNum = {s,
                                             s,
                                             center_ch0 - delta_ch0,
                                             center_ch0 + delta_ch0,
                                             center_ch1 - delta_ch1,
                                             center_ch1 + delta_ch1,
                                             /*mean_ratio - 0.5 * sigma_ratio, // lead collimator + Hamamatsu SiPMs, 1x1 mm fibers (Jan 2018)
                                             mean_ratio + 0.5 * sigma_ratio*/
                                             mean_ratio - 2 * sigma_ratio, // lead collimator + Hamamatsu SiPMs, 2x3 mm fibers (March 2018)
                                             mean_ratio + 2 * sigma_ratio}; 
     
            cut = SFDrawCommands::GetCut(SFCutType::kT0DiffECut, customNum);
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("Hamamatsu"))
        {
            parNo       = 1;
            mean_ratio  = fRatios[i]->GetFunction("fGauss")->GetParameter(parNo);
            sigma_ratio = fRatios[i]->GetFunction("fGauss")->GetParameter(parNo + 1);
            
            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f
            // && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0,
            // center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-3*sigma_ratio,
            // mean_ratio+3*sigma_ratio);   //changed here for smaller cut
            
            std::vector<double> customNum = {s,
                                             s,
                                             center_ch0 - 3 * delta_ch0,
                                             center_ch0 + 3 * delta_ch0,
                                             center_ch1 - 3 * delta_ch1,
                                             center_ch1 + 3 * delta_ch1,
                                             mean_ratio - 3 * sigma_ratio,
                                             mean_ratio + 3 * sigma_ratio};
            cut = SFDrawCommands::GetCut(SFCutType::kT0DiffECut, customNum);
        }
        else if (collimator.Contains("Electronic") && sipm.Contains("SensL"))
        {
            const_1 = fRatios[i]->GetFunction("fDGauss")->GetParameter(0);
            const_2 = fRatios[i]->GetFunction("fDGauss")->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
            mean_ratio  = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo);
            sigma_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo + 1);
            
            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f
            // && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0,
            // center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-2*sigma_ratio,
            // mean_ratio+2*sigma_ratio);
            
            std::vector<double> customNum = {s,
                                             s,
                                             center_ch0 - 3 * delta_ch0,
                                             center_ch0 + 3 * delta_ch0,
                                             center_ch1 - 3 * delta_ch1,
                                             center_ch1 + 3 * delta_ch1,
                                             mean_ratio - 2 * sigma_ratio,
                                             mean_ratio + 2 * sigma_ratio};
            cut = SFDrawCommands::GetCut(SFCutType::kT0DiffECut, customNum);
        }
        else if (testBench == "PMI")
        {
            const_1 = fRatios[i]->GetFunction("fDGauss")->GetParameter(0);
            const_2 = fRatios[i]->GetFunction("fDGauss")->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
            mean_ratio  = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo);
            sigma_ratio = fRatios[i]->GetFunction("fDGauss")->GetParameter(parNo + 1);
            
            // cut = Form("ch_0.fT0>0 && ch_1.fT0>0 && ch_0.fT0<590 && ch_1.fT0<590 && ch_0.fPE>%f
            // && ch_0.fPE<%f && ch_1.fPE>%f && ch_1.fPE<%f && log(sqrt(ch_1.fPE/ch_0.fPE))>%f &&
            // log(sqrt(ch_1.fPE/ch_0.fPE))<%f",  center_ch0-3*delta_ch0, center_ch0+3*delta_ch0,
            // center_ch1-3*delta_ch1, center_ch1+3*delta_ch1, mean_ratio-2*sigma_ratio,
            // mean_ratio+2*sigma_ratio);
            
            std::vector<double> customNum = {center_ch0 - 3 * delta_ch0,
                                             center_ch0 + 3 * delta_ch0,
                                             center_ch1 - 3 * delta_ch1,
                                             center_ch1 + 3 * delta_ch1,
                                             mean_ratio - 2 * sigma_ratio,
                                             mean_ratio + 2 * sigma_ratio};
            cut = SFDrawCommands::GetCut(SFCutType::kPMIT0DiffECut, customNum);
            
        }

        fT0DiffECut.push_back(fData->GetCustomHistogram(SFSelectionType::kT0Difference,
                                                        cut, measIDs[i]));
        mean  = fT0DiffECut[i]->GetMean();
        sigma = fT0DiffECut[i]->GetRMS();
        fT0DiffECut[i]->Fit(fun, "Q", "", mean - 5 * sigma, mean + 5 * sigma);
        
        fTSigmaECutGraph->SetPoint(i, positions[i], fun->GetParameter(2)); // Timing resolution as sigma, if FWHM needed multiply by f
        fTSigmaECutGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                                        fun->GetParError(2)); // FWHM - multiply by f

        fTResECutGraph->SetPoint(i, positions[i], fun->GetParameter(1));
        fTResECutGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench),
                             fun->GetParameter(2)); // FWHM - multiply by f

        tResSum    += fun->GetParameter(2) * (1. / pow(fun->GetParError(2), 2));
        tResSumErr += 1. / pow(fun->GetParError(2), 2);
    }

    tResAv    = tResSum / tResSumErr;
    tResAvErr = sqrt(1. / tResSumErr);
    
    fResultsECut->AddResult(SFResultTypeNum::kTimeRes, tResAv, tResAvErr);
    fResultsECut->AddObject(SFResultTypeObj::kTimeResGraph, fTResECutGraph);
    fResultsECut->AddObject(SFResultTypeObj::kTimeSigGraph, fTSigmaECutGraph);

    std::cout << "Average timing resolution: " << tResAv << " +/- "
              << tResAvErr << " ns" << std::endl;

    return true;
}
//------------------------------------------------------------------
/// Returns vector of histograms, which contains timing resolution
/// spectra for whole series.
/// \param type - type of analysis (1 - with cut, 0 - no cut)
std::vector<TH1D*> SFTimingRes::GetT0Diff(bool type)
{

    if ((type == 0 && fT0Diff.empty()) || (type == 1 && fT0DiffECut.empty()))
    {
        std::cerr << "##### Error in SFTimingRes::GetT0Diff()!" << std::endl;
        std::cerr << "No spectra available!" << std::endl;
    }

    std::vector<TH1D*> hist;
    
    if (type == 0)
        hist = fT0Diff;
    else if (type == 1)
        hist = fT0DiffECut;
    
    return hist;
}
//------------------------------------------------------------------
/// Returns vector of ratio histograms used for cut on scattered events.
std::vector<TH1D*> SFTimingRes::GetRatios(void)
{

    if (fRatios.empty())
    {
        std::cout << "##### Error in SFTimingRes::GetRatios()! Empty vector!" << std::endl;
        std::abort();
    }
    return fRatios;
}
//------------------------------------------------------------------
/// Returns vector of PE spectra used for 511 keV energy cut.
/// \param ch - channel number (0 or 1).
std::vector<TH1D*> SFTimingRes::GetSpectra(int ch)
{

    if ((ch == 0 && fSpecCh0.empty()) || (ch == 1 && fSpecCh1.empty()))
    {
        std::cerr << "##### Error in SFTimingRes::GetSpectra()!" << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    std::vector<TH1D*> spectra;
    
    if (ch == 0)
        spectra = fSpecCh0;
    else if (ch == 1)
        spectra = fSpecCh1;
    else
    {
        std::cerr << "##### Error in SFTimingRes::GetSpectra()" << std::endl;
        std::cerr << "Incorrect channel number!" << std::endl;
        std::abort();
    }
    
    return spectra;
}
//------------------------------------------------------------------
std::vector<SFResults*> SFTimingRes::GetResults(void)
{
    
    if(fResults == nullptr || 
       fResultsECut == nullptr)
    {
        std::cerr << "##### Error in SFTimingRes::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(2);
    
    results[0] = fResults;
    results[1] = fResultsECut;
    
    return results;
}
//------------------------------------------------------------------
/// Prints details of SFTimingRes class object.
void SFTimingRes::Print(void)
{
    std::cout << "\n--------------------------------------------" << std::endl;
    std::cout << "This is print out of SFTimingRes class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "--------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
