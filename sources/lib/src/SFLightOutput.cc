// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFLightOutput.cc           *
// *             Jonas Kasper              *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFLightOutput.hh"

ClassImp(SFLightOutput);

//------------------------------------------------------------------
SFLightOutput::SFLightOutput(int seriesNo) : fSeriesNo(seriesNo),
                                             fPDE(-1),
                                             fCrossTalk(-1),
                                             fLightOutGraph(nullptr),
                                             fLightOutCh0Graph(nullptr),
                                             fLightOutCh1Graph(nullptr),
                                             fLightColGraph(nullptr),
                                             fLightColCh0Graph(nullptr),
                                             fLightColCh1Graph(nullptr),
                                             fData(nullptr),
                                             fAtt(nullptr),
                                             fResultsLightOutCh0(nullptr),
                                             fResultsLightOutCh1(nullptr),
                                             fResultsLightOutSum(nullptr),
                                             fResultsLightColCh0(nullptr),
                                             fResultsLightColCh1(nullptr),
                                             fResultsLightColSum(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cout << message << std::endl;
        throw "##### Exception in SFLightOutput constructor!";
    }

    int     npoints     = fData->GetNpoints();
    TString testBench   = fData->GetTestBench();
    TString description = fData->GetDescription();
    
    if (!description.Contains("Regular series"))
    {
        std::cout << "##### Warning in SFLightOutput constructor!" << std::endl;
        std::cout << "Calculating light output for non-regular series!" << std::endl;
    }

    TString cutCh0 = " ";
    TString cutCh1 = " ";
    
    if (testBench == "PMI")
    {
        cutCh0 = SFDrawCommands::GetCut(SFCutType::kPMISpecCh0);
        cutCh1 = SFDrawCommands::GetCut(SFCutType::kPMISpecCh1);
        fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::kPMICharge, cutCh0);
        fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::kPMICharge, cutCh1);
    }
    else
    {
        double s = SFTools::GetSigmaBL(fData->GetSiPM());
        std::vector<double> sigma  = {s};
        cutCh0 = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
        cutCh1 = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);
        fSpectraCh0 = fData->GetSpectra(0, SFSelectionType::kPE, cutCh0);
        fSpectraCh1 = fData->GetSpectra(1, SFSelectionType::kPE, cutCh1);
    }

    for (int i = 0; i < npoints; i++)
    {
        fPFCh0.push_back(new SFPeakFinder(fSpectraCh0[i], 0));
        fPFCh1.push_back(new SFPeakFinder(fSpectraCh1[i], 0));
    }

    fCrossTalk = GetCrossTalk();
    fPDE       = GetPDE();

    fResultsLightOutCh0 = new SFResults(Form("LightOutResults_S%i_ch0", fSeriesNo));
    fResultsLightOutCh1 = new SFResults(Form("LightOutResults_S%i_ch1", fSeriesNo));
    fResultsLightOutSum = new SFResults(Form("LightOutResults_S%i_sum", fSeriesNo));
    
    fResultsLightColCh0 = new SFResults(Form("LightColResults_S%i_ch0", fSeriesNo));
    fResultsLightColCh1 = new SFResults(Form("LightColResults_S%i_ch1", fSeriesNo));
    fResultsLightColSum = new SFResults(Form("LightColResults_S%i_sum", fSeriesNo));
}
//------------------------------------------------------------------
double SFLightOutput::GetCrossTalk(void)
{

    double  crossTalk = -1;
    TString SiPM      = fData->GetSiPM();
    double  overvol   = fData->GetOvervoltage();

    TString fname = std::string(getenv("SFDATA")) + "/DB/SiPMs_Data.root";
    TFile*  file  = new TFile(fname, "READ");
    TGraph* gCrossTalk;

    if (SiPM == "Hamamatsu") { gCrossTalk = (TGraph*)file->Get("gHamamatsu_CT"); }
    else if (SiPM == "SensL")
    {
        crossTalk = 0.03;
        return crossTalk;
    }
    else
    {
        std::cerr << "##### Error in SFLightOutput::GetCrossTalk()!" << std::endl;
        std::cerr << "Unknown SiPM type! Please check!" << std::endl;
        std::abort();
    }

    if (gCrossTalk == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::GetCrossTalk()!" << std::endl;
        std::cerr << "Couldn't get requested cross talk graph!" << std::endl;
        std::abort();
    }

    crossTalk = gCrossTalk->Eval(overvol) / 100.;
    file->Close();

    std::cout << "\t----- Cross talk is: " << crossTalk << std::endl;

    return crossTalk;
}
//------------------------------------------------------------------
double SFLightOutput::GetPDE(void)
{

    double  PDE     = 0;
    double  overvol = fData->GetOvervoltage();
    TString fiber   = fData->GetFiber();
    TString SiPM    = fData->GetSiPM();

    TString fname = std::string(getenv("SFDATA")) + "/DB/SiPMs_Data.root";
    TFile*  file  = new TFile(fname, "READ");

    TGraph* gPDEvsVol;
    TH1F*   hPDEvsWavelen;
    TH1F*   hLightOutvsWavelen;

    if (fiber.Contains("LuAG:Ce")) { hLightOutvsWavelen = (TH1F*)file->Get("hLuAG_LOvsWave"); }
    else if (fiber.Contains("LYSO:Ce"))
    {
        hLightOutvsWavelen = (TH1F*)file->Get("hLYSO_LOvsWave");
    }
    else if (fiber.Contains("GAGG:Ce"))
    {
        hLightOutvsWavelen = (TH1F*)file->Get("hGAGG_LOvsWave");
    }
    else
    {
        std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
        std::cerr << "Unknown fiber type! Please check!" << std::endl;
    }

    if (SiPM == "Hamamatsu")
    {
        gPDEvsVol     = (TGraph*)file->Get("gHamamatsu_PDEvsVol");
        hPDEvsWavelen = (TH1F*)file->Get("hHamamatsu_PDEvsWave");
    }
    else if (SiPM == "SensL")
    {
        gPDEvsVol     = (TGraph*)file->Get("gSensL_PDEvsVol");
        hPDEvsWavelen = (TH1F*)file->Get("hSensL_PDEvsWave");
    }
    else
    {
        std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
        std::cerr << "Unknown SiPM type! Please check!" << std::endl;
    }

    if (gPDEvsVol == nullptr || hPDEvsWavelen == nullptr || hLightOutvsWavelen == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::GetPDE()!" << std::endl;
        std::cerr << "Couldn't get one of requested graphs/histograms!" << std::endl;
        std::abort();
    }

    int    nbins = hPDEvsWavelen->GetXaxis()->GetNbins();
    double f =
        gPDEvsVol->Eval(overvol) / hPDEvsWavelen->GetBinContent(hPDEvsWavelen->GetMaximumBin());

    for (int i = 1; i < nbins + 1; i++)
    {
        PDE += hLightOutvsWavelen->GetBinContent(i) * hPDEvsWavelen->GetBinContent(i) / 100.;
    }

    double PDE_final = PDE * f;

    //----- input data canvas
    TLatex  text;
    TString cname = Form("fInputData_S%i", fSeriesNo);

    TCanvas* can = new TCanvas(cname, cname, 1800, 1800);
    can->Divide(2, 2);

    can->cd(1);
    gPad->SetGrid(1, 1);
    gPDEvsVol->Draw("AP");

    can->cd(2);
    gPad->SetGrid(1, 1);
    hLightOutvsWavelen->Draw("hist");

    can->cd(3);
    text.DrawLatex(0.1, 0.9, Form("PDE scaling factor f = %.2f", f));
    text.DrawLatex(0.1, 0.8, Form("PDE from PDE vs voltage curve = %.2f", gPDEvsVol->Eval(overvol)));
    text.DrawLatex(0.1, 0.7, Form("SiPM maximum sensitivity = %.2f",
                   hPDEvsWavelen->GetBinContent(hPDEvsWavelen->GetMaximumBin())));
    text.DrawLatex(0.1, 0.6, Form("Nbins of emission spectrum = %i", hLightOutvsWavelen->GetNbinsX()));
    text.DrawLatex(0.1, 0.5, Form("Nbins of SiPM sensitivity curve = %i", hPDEvsWavelen->GetNbinsX()));
    text.DrawLatex(0.1, 0.4, Form("Emission and SiPMs sensitivity convolution = %.2f", PDE));
    text.DrawLatex(0.1, 0.3, Form("Final PDE value = %.2f", PDE_final));

    can->cd(4);
    gPad->SetGrid(1, 1);
    hPDEvsWavelen->Draw();

    fInputData     = (TCanvas*)can->Clone("fInputData");
    TCanvas* dummy = new TCanvas("dummy", "dummy", 100, 100);
    //-----

    file->Close();

    std::cout << "\t----- PDE is: " << PDE << std::endl;

    return PDE_final;
}
//------------------------------------------------------------------
/// Default destructor.
SFLightOutput::~SFLightOutput()
{
    if (fData != nullptr) delete fData;
    if (fAtt != nullptr) delete fAtt;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightOut(void)
{

    std::cout << "\n----- Inside SFLightOutput::CalculateLightOut()" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << std::endl;

    if (fLightOutCh0Graph == nullptr || fLightOutCh1Graph == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::CalculateLightOut()" << std::endl;
        std::cerr << "fLightOutCh0Graph and fLightOutCh1Graph don't exist!" << std::endl;
        std::cerr << "Run analysis for channels 0 and 1 first!" << std::endl;
        return false;
    }

    int                 npoints    = fData->GetNpoints();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             desc       = fData->GetDescription();
    std::vector<double> positions  = fData->GetPositions();

    TString       gname = Form("LightOutSum_S%i", fSeriesNo);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("light output [PE/MeV]");
    graph->SetName(gname);
    graph->SetTitle(Form("Summed Light Output S%i", fSeriesNo));
    graph->SetMarkerStyle(4);

    double lightOut      = 0;
    double lightOutErr   = 0;
    double lightOutAv    = 0;
    double lightOutAvErr = 0;
    double xCh0, xCh1, yCh0, yCh1;

    for (int i = 0; i < npoints; i++)
    {
        fLightOutCh0Graph->GetPoint(i, xCh0, yCh0);
        fLightOutCh1Graph->GetPoint(i, xCh1, yCh1);
        lightOut = yCh0 + yCh1;
        lightOutErr = sqrt(pow(fLightOutCh0Graph->GetErrorY(i), 2) + 
                           pow(fLightOutCh1Graph->GetErrorY(i), 2));
    
        graph->SetPoint(i, positions[i], lightOut);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightOutErr);
        lightOutAv += lightOut * (1. / pow(lightOutErr, 2));
        lightOutAvErr += 1. / pow(lightOutErr, 2);
    }

    lightOutAv     = lightOutAv / lightOutAvErr;
    lightOutAvErr  = sqrt(1. / lightOutAvErr);
    fLightOutGraph = graph;

    fResultsLightOutSum->AddResult(SFResultTypeNum::kLight, lightOutAv, lightOutAvErr);
    fResultsLightOutSum->AddObject(SFResultTypeObj::kLightGraph, fLightOutGraph);

    std::cout << "Averaged and summed light output: " << lightOutAv << " +/- " << lightOutAvErr
              << " ph/MeV" << std::endl;

    return true;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightOut(int ch)
{

    std::cout << "\n----- Inside SFLightOutput::CalculateLightOut()" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << "\t Channel: " << ch << std::endl;

    int                        npoints    = fData->GetNpoints();
    double                     fiberLen   = fData->GetFiberLength();
    TString                    collimator = fData->GetCollimator();
    TString                    testBench  = fData->GetTestBench();
    std::vector<double>        positions  = fData->GetPositions();
    std::vector<SFPeakFinder*> peakFin;


    SFAttenuation* att;

    try
    {
        att = new SFAttenuation(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cout << message << std::endl;
        throw "##### Exception in SFLightOutput::CalculateLightOut()!";
        return false;
    }

    att->AttCombinedCh();

    std::vector<SFResults*> results_tmp = att->GetResults();
    SFResults *results = results_tmp[2]; // combined channels/pol1

    for (int i = 0; i < npoints; i++)
    {
        if (ch == 0)
            peakFin = fPFCh0;
        else if (ch == 1)
            peakFin = fPFCh1;
        else
        {
            std::cerr << "##### Error in SFLightOutput::CalculateLightOut!" << std::endl;
            std::cerr << "Incorrect channel number" << std::endl;
            return false;
        }
    }

    TString gname = Form("LightOut_S%i_Ch%i", fSeriesNo, ch);

    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("light output [PE/MeV]");
    graph->SetTitle(Form("Light Output S%i ch%i", fSeriesNo, ch));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    SFResults* parameters;
    double     lightOut      = 0;
    double     lightOutErr   = 0;
    double     lightOutAv    = 0;
    double     lightOutAvErr = 0;
    double     distance      = 0;

    for (int i = 0; i < npoints; i++)
    {

        peakFin[i]->FindPeakFit();
        if (ch == 0) distance = positions[i];
        if (ch == 1) distance = fiberLen - positions[i];
        parameters = peakFin[i]->GetResults();

        lightOut = parameters->GetValue(SFResultTypeNum::kPeakPosition) * (1 - fCrossTalk) / fPDE /
                   0.511 / TMath::Exp(-distance / results->GetValue(SFResultTypeNum::kLambda));
        lightOutErr = sqrt((pow(lightOut, 2) * pow(parameters->GetValue(SFResultTypeNum::kPeakSigma), 2) /
                      pow(parameters->GetValue(SFResultTypeNum::kPeakPosition), 2)) +
                      (pow(lightOut, 2) * pow(results->GetUncertainty(SFResultTypeNum::kLambda), 2) /
                       pow(results->GetValue(SFResultTypeNum::kLambda), 4)));
        graph->SetPoint(i, positions[i], lightOut);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightOutErr);
        lightOutAv += lightOut * (1. / pow(lightOutErr, 2));
        lightOutAvErr += (1. / pow(lightOutErr, 2));
    }

    lightOutAv    = lightOutAv / lightOutAvErr;
    lightOutAvErr = sqrt(1. / lightOutAvErr);

    if (ch == 0)
    {
        fLightOutCh0Graph = graph;
        fResultsLightOutCh0->AddResult(SFResultTypeNum::kLight, lightOutAv, lightOutAvErr);
        fResultsLightOutCh0->AddObject(SFResultTypeObj::kLightGraph, fLightOutCh0Graph);
    }
    else if (ch == 1)
    {
        fLightOutCh1Graph = graph;
        fResultsLightOutCh1->AddResult(SFResultTypeNum::kLight, lightOutAv, lightOutAvErr);
        fResultsLightOutCh1->AddObject(SFResultTypeObj::kLightGraph, fLightOutCh1Graph);
    }

    std::cout << "Average light output for channel " << ch << ": " << lightOutAv << " +/- "
              << lightOutAvErr << " ph/MeV \n"
              << std::endl;

    return true;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightCol(void)
{

    std::cout << "\n----- Inside SFLightOutput::CalculateLCol()" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << std::endl;

    if (fLightColCh0Graph == nullptr || fLightColCh1Graph == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::CalculateLightCol()" << std::endl;
        std::cerr << "fLightColCh0Graph and fLightColCh1Graph don't exist!" << std::endl;
        std::cerr << "Run analysis for channels 0 and 1 first!" << std::endl;
        return false;
    }

    int                 npoints    = fData->GetNpoints();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    std::vector<double> positions  = fData->GetPositions();

    TString       gname = Form("LCol_s%i_sum", fSeriesNo);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("collected light [PE/MeV]");
    graph->SetTitle(Form("Summed Light Collection S%i", fSeriesNo));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    double lightCol      = 0;
    double lightColErr   = 0;
    double lightColAv    = 0;
    double lightColAvErr = 0;
    double xCh0, xCh1, yCh0, yCh1;
    
    for (int i = 0; i < npoints; i++)
    {
        fLightColCh0Graph->GetPoint(i, xCh0, yCh0);
        fLightColCh1Graph->GetPoint(i, xCh1, yCh1);
        lightCol = yCh0 + yCh1;
        lightColErr = sqrt(pow(fLightColCh0Graph->GetErrorY(i), 2) + 
                           pow(fLightColCh1Graph->GetErrorY(i), 2));
        
        graph->SetPoint(i, positions[i], lightCol);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightColErr);
        lightColAv += lightCol * (1. / pow(lightColErr, 2));
        lightColAvErr += 1. / pow(lightColErr, 2);
    }

    lightColAv     = lightColAv / lightColAvErr;
    lightColAvErr  = sqrt(1. / lightColAvErr);
    fLightColGraph = graph;

    fResultsLightColSum->AddResult(SFResultTypeNum::kLight, lightColAv, lightColAvErr);
    fResultsLightColSum->AddObject(SFResultTypeObj::kLightGraph, fLightColGraph);

    std::cout << "Average light collection for this series: " << lightColAv << " +/- "
              << lightColAvErr << " ph/MeV" << std::endl;

    return true;
}
//------------------------------------------------------------------
bool SFLightOutput::CalculateLightCol(int ch)
{

    std::cout << "\n----- Inside SFLightOutput::CalculateLightCol() for series " << fSeriesNo
              << std::endl;
    std::cout << "----- Analyzing channel " << ch << std::endl;

    int                        npoints    = fData->GetNpoints();
    TString                    collimator = fData->GetCollimator();
    TString                    testBench  = fData->GetTestBench();
    std::vector<double>        positions  = fData->GetPositions();
    std::vector<SFPeakFinder*> tempPF;

    TString       gname = Form("LCol_s%i_ch%i", fSeriesNo, ch);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("collected light [PE/MeV]");
    graph->SetTitle(Form("Light Collection S%i ch%i", fSeriesNo, ch));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);
        
    if (ch == 0)
        tempPF = fPFCh0;
    else if (ch == 1)
        tempPF = fPFCh1;
    else
    {
         std::cerr << "##### Error in SFLightOutput::CalculateLightCol()!" << std::endl;
         std::cerr << "Incorrect channel number! Please check!" << std::endl;
         return false;
    }

    SFResults* peak_par;
    double     lightCol      = 0;
    double     lightColErr   = 0;
    double     lightColAv    = 0;
    double     lightColAvErr = 0;
    
    for (int i = 0; i < npoints; i++)
    {

        peak_par    = tempPF[i]->GetResults();
        lightCol    = peak_par->GetValue(SFResultTypeNum::kPeakPosition) / 0.511;
        lightColErr = peak_par->GetValue(SFResultTypeNum::kPeakSigma) / 0.511;
        graph->SetPoint(i, positions[i], lightCol);
        graph->SetPointError(i, SFTools::GetPosError(collimator, testBench), lightColErr);

        lightColAv += lightCol * (1. / pow(lightColErr, 2));
        lightColAvErr += (1. / pow(lightColErr, 2));
    }

    lightColAv    = lightColAv / lightColAvErr;
    lightColAvErr = sqrt(1. / lightColAvErr);

    if (ch == 0)
    {
        fLightColCh0Graph = graph;
        fResultsLightColCh0->AddObject(SFResultTypeObj::kLightGraph, fLightColCh0Graph);
        fResultsLightColCh0->AddResult(SFResultTypeNum::kLight, lightColAv, lightColAvErr);
    }
    else if (ch == 1)
    {
        fLightColCh1Graph = graph;
        fResultsLightColCh1->AddObject(SFResultTypeObj::kLightGraph, fLightColCh1Graph);
        fResultsLightColCh1->AddResult(SFResultTypeNum::kLight, lightColAv, lightColAvErr);
    }

    std::cout << "Average light collection for this series and channel " << ch << ": " << lightColAv
              << " +/- " << lightColAvErr << " ph/MeV" << std::endl;

    return true;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFLightOutput::GetSpectra(int ch)
{

    if ((ch == 0 && fSpectraCh0.empty()) || (ch == 1 && fSpectraCh1.empty()))
    {
        std::cerr << "##### Error in SFLightOutput::GetSpectra() for ch" << ch << std::endl;
        std::cerr << "No spectra available!" << std::endl;
        std::abort();
    }

    if (ch == 0)
        return fSpectraCh0;
    else if (ch == 1)
        return fSpectraCh1;
    else
    {
        std::cerr << "##### Error in SFLightOutput::GetSpectra for ch" << ch << std::endl;
        std::cerr << "Incorrect channel number!" << std::endl;
        std::abort();
    }
}
//------------------------------------------------------------------
std::vector<SFResults*> SFLightOutput::GetLOResults()
{
    if (fResultsLightOutCh0 == nullptr ||
        fResultsLightOutCh1 == nullptr ||
        fResultsLightOutSum == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::GetLOResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(3);
    results[0] = fResultsLightOutCh0;
    results[1] = fResultsLightOutCh1;
    results[2] = fResultsLightOutSum;
    
    return results;
}
//------------------------------------------------------------------
std::vector<SFResults*> SFLightOutput::GetLCResults()
{
    if (fResultsLightColCh0 == nullptr ||
        fResultsLightColCh1 == nullptr ||
        fResultsLightColSum == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::GetLCResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }
    
    std::vector<SFResults*> results(3);
    results[0] = fResultsLightColCh0;
    results[1] = fResultsLightColCh1;
    results[2] = fResultsLightColSum;
    
    return results;
}
//------------------------------------------------------------------
/// Prints details of the SFLightOutput class object.
void SFLightOutput::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFLightOutput class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
