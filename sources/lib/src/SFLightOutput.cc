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

//------------------------------------------------------------------
double SFLightOutput::GetCrossTalk(TString SiPM, double overvol)
{

    double  crossTalk = -1;

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
auto SFLightOutput::GetPDE(TString fiber, TString SiPM, double overvol) -> std::tuple<double, TCanvas*>
{

    double  PDE     = 0;

    TString fname = std::string(getenv("SFDATA")) + "/DB/SiPMs_Data.root";
    TFile*  file  = new TFile(fname, "READ");

    TGraph* gPDEvsVol = nullptr;
    TH1F*   hPDEvsWavelen = nullptr;
    TH1F*   hLightOutvsWavelen = nullptr;

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

    TCanvas* can = new TCanvas("input_data", "input_data", 1800, 1800);
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
    //-----

    file->Close();

    std::cout << "\t----- PDE is: " << PDE << std::endl;

    return {PDE_final, can};
}
//------------------------------------------------------------------
SFResults* SFLightOutput::CalculateLightOut(TGraphErrors *gLightOutL, TGraphErrors *gLightOutR)
{

    if (gLightOutL == nullptr || gLightOutR == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "gLightOutL or gLightOutR don't exist!" << std::endl;
        std::abort();
    }

    int npoints = gLightOutL->GetN();
    double pos_uncert = gLightOutL->GetErrorX(0);
    double *positions = gLightOutL->GetX();

    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("light output [PE/MeV]");
    graph->SetName("LightOutSum");
    graph->SetTitle("Summed Light Output");
    graph->SetMarkerStyle(4);

    double lightOut      = 0;
    double lightOutErr   = 0;
    double lightOutAv    = 0;
    double lightOutAvErr = 0;
    double xCh0, xCh1, yCh0, yCh1;

    for (int i = 0; i < npoints; i++)
    {
        gLightOutL->GetPoint(i, xCh0, yCh0);
        gLightOutR->GetPoint(i, xCh1, yCh1);
        lightOut = yCh0 + yCh1;
        lightOutErr = sqrt(pow(gLightOutL->GetErrorY(i), 2) + 
                           pow(gLightOutR->GetErrorY(i), 2));
    
        graph->SetPoint(i, positions[i], lightOut);
        graph->SetPointError(i, pos_uncert, lightOutErr);
        lightOutAv += lightOut * (1. / pow(lightOutErr, 2));
        lightOutAvErr += 1. / pow(lightOutErr, 2);
    }

    lightOutAv     = lightOutAv / lightOutAvErr;
    lightOutAvErr  = sqrt(1. / lightOutAvErr);

    SFResults *results = new SFResults("LightOutSum");
    results->AddResult(SFResultTypeNum::kLight, lightOutAv, lightOutAvErr);
    results->AddObject(SFResultTypeObj::kLightGraph, graph);

    std::cout << "Averaged and summed light output: " << lightOutAv << " +/- " << lightOutAvErr
              << " ph/MeV" << std::endl;

    return results;
}
//------------------------------------------------------------------
SFResults* SFLightOutput::CalculateLightOut(TString side, std::vector<double> positions,
                                            std::vector<TH1D*> spectra, std::vector<TString> path,
                                            double pos_uncert, double fiberLen, double overvol,
                                            TString fiber, TString sipm, SFResults* attResults)
{

    if (!(side == "l" || side == "L" || side == "r" || side == "R"))
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "Incorrect side. Possible options are: l/L/r/R" << std::endl;
        std::abort();
    }
    
    if(spectra.empty() || spectra[0] == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "There's something wrong with std::vector<TH1D*> spectra! Please check!" << std::endl;
        std::abort();
    }

    int npoints = positions.size();

    TString gname = Form("LightOut_%s", side.Data());
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("light output [PE/MeV]");
    graph->SetTitle(Form("Light Output %s", side.Data()));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    double     lightOut      = 0;
    double     lightOutErr   = 0;
    double     lightOutAv    = 0;
    double     lightOutAvErr = 0;
    double     distance      = 0;

    double crossTalk = GetCrossTalk(sipm, overvol);
    auto pde_tuple = GetPDE(fiber, sipm, overvol);
    double pde = std::get<0>(pde_tuple);
    double lambda = attResults->GetValue(SFResultTypeNum::kLambda);
    double lambda_err = attResults->GetValue(SFResultTypeNum::kLambda);
    
    for (int i = 0; i < npoints; i++)
    {
        auto parameters = std::unique_ptr<SFResults>(SFPeakFinder::FindPeakFit(spectra[i], path[i], 0, 0));
        
        if (side == "l" || side == "L") distance = positions[i];
        if (side == "r" || side == "R") distance = fiberLen - positions[i];

        lightOut = parameters->GetValue(SFResultTypeNum::kPeakPosition) * (1 - crossTalk) / pde /
                   0.511 / TMath::Exp(-distance / lambda);
        lightOutErr = sqrt((pow(lightOut, 2) * pow(parameters->GetValue(SFResultTypeNum::kPeakSigma), 2) /
                      pow(parameters->GetValue(SFResultTypeNum::kPeakPosition), 2)) +
                      (pow(lightOut, 2) * pow(lambda_err, 2) / pow(lambda, 4)));
        graph->SetPoint(i, positions[i], lightOut);
        graph->SetPointError(i, pos_uncert, lightOutErr);
        lightOutAv += lightOut * (1. / pow(lightOutErr, 2));
        lightOutAvErr += (1. / pow(lightOutErr, 2));
    }

    lightOutAv    = lightOutAv / lightOutAvErr;
    lightOutAvErr = sqrt(1. / lightOutAvErr);

    
    SFResults *results = new SFResults(Form("LightOut_%s", side.Data()));
    results->AddResult(SFResultTypeNum::kLight, lightOutAv, lightOutAvErr);
    results->AddObject(SFResultTypeObj::kPDECan, std::get<1>(pde_tuple));
    results->AddObject(SFResultTypeObj::kLightGraph, graph);

    std::cout << "Average light output for this series and side " << side << ": "
              << lightOutAv << " +/- " << lightOutAvErr << " ph/MeV \n" << std::endl;

    return results;
}
//------------------------------------------------------------------
SFResults* SFLightOutput::CalculateLightCol(TGraphErrors* gLightColL, TGraphErrors* gLightColR)
{

    if (gLightColL == nullptr || gLightColR == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "gLightColL or gLightColR don't exist!" << std::endl;
        std::abort();
    }

    int npoints = gLightColL->GetN();
    double pos_uncert = gLightColL->GetErrorX(0);
    double *positions = gLightColL->GetX();

    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("collected light [PE/MeV]");
    graph->SetTitle("Summed Light Collection");
    graph->SetName("LCol_sum");
    graph->SetMarkerStyle(4);

    double lightCol      = 0;
    double lightColErr   = 0;   
    double lightColAv    = 0;
    double lightColAvErr = 0;
    double xCh0, xCh1, yCh0, yCh1;
    
    for (int i = 0; i < npoints; i++)
    {
        gLightColL->GetPoint(i, xCh0, yCh0);
        gLightColR->GetPoint(i, xCh1, yCh1);
        lightCol = yCh0 + yCh1;
        lightColErr = sqrt(pow(gLightColL->GetErrorY(i), 2) + 
                           pow(gLightColR->GetErrorY(i), 2));
        
        graph->SetPoint(i, positions[i], lightCol);
        graph->SetPointError(i, pos_uncert, lightColErr);
        lightColAv += lightCol * (1. / pow(lightColErr, 2));
        lightColAvErr += 1. / pow(lightColErr, 2);
    }

    lightColAv     = lightColAv / lightColAvErr;
    lightColAvErr  = sqrt(1. / lightColAvErr);

    SFResults* results = new SFResults("LightColSum");
    results->AddResult(SFResultTypeNum::kLight, lightColAv, lightColAvErr);
    results->AddObject(SFResultTypeObj::kLightGraph, graph);

    std::cout << "Average summed light collection for this series: " << lightColAv << " +/- "
              << lightColAvErr << " ph/MeV" << std::endl;

    return results;
}
//------------------------------------------------------------------
SFResults* SFLightOutput::CalculateLightCol(TString side, std::vector<double> positions,
                                            std::vector<TH1D*> spectra, std::vector<TString> path,
                                            double pos_uncert)
{
        
    if (!(side == "l" || side == "L" || side == "r" || side == "R"))
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "Incorrect side. Possible options are: l/L/r/R" << std::endl;
        std::abort();
    }
    
    if(spectra.empty() || spectra[0] == nullptr)
    {
        std::cerr << "##### Error in SFLightOutput::" << __func__ << std::endl;
        std::cerr << "There's something wrong with std::vector<TH1D*> spectra! Please check!" << std::endl;
        std::abort();
    }
    
    int npoints = positions.size();

    TString       gname = Form("LCol_%s", side.Data());
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("collected light [PE/MeV]");
    graph->SetTitle(Form("Light Collection %s", side.Data()));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    double     lightCol      = 0;
    double     lightColErr   = 0;
    double     lightColAv    = 0;
    double     lightColAvErr = 0;
    
    for (int i = 0; i < npoints; i++)
    {
        auto peak_par = std::unique_ptr<SFResults>(SFPeakFinder::FindPeakFit(spectra[i], path[i], 0, 0));
        lightCol    = peak_par->GetValue(SFResultTypeNum::kPeakPosition) / 0.511;
        lightColErr = peak_par->GetValue(SFResultTypeNum::kPeakSigma) / 0.511;
        graph->SetPoint(i, positions[i], lightCol);
        graph->SetPointError(i, pos_uncert, lightColErr);

        lightColAv += lightCol * (1. / pow(lightColErr, 2));
        lightColAvErr += (1. / pow(lightColErr, 2));
    }

    lightColAv    = lightColAv / lightColAvErr;
    lightColAvErr = sqrt(1. / lightColAvErr);

    SFResults* results = new SFResults(Form("LightCol_%s", side.Data()));
    results->AddObject(SFResultTypeObj::kLightGraph, graph);
    results->AddResult(SFResultTypeNum::kLight, lightColAv, lightColAvErr);
    
    std::cout << "Average light collection for this series and side " << side << ": " << lightColAv
              << " +/- " << lightColAvErr << " ph/MeV" << std::endl;

    return results;
}
//------------------------------------------------------------------
