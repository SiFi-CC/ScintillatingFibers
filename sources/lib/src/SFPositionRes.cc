// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionRes.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFPositionRes.hh"

ClassImp(SFPositionRes);

const double ampMax = 660;

//------------------------------------------------------------------
SFPositionRes::SFPositionRes(int seriesNo) : fSeriesNo(seriesNo),
                                             fData(nullptr),
                                             fAtt(nullptr),
                                             fPosVsMLRGraph(nullptr),
                                             fResultsPol3(nullptr),
                                             fResultsPol1(nullptr)
{

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFPositionRes constructor!";
    }

    TString desc = fData->GetDescription();

    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in SFPositionRes constructor! Non-regular series!" << std::endl;
        throw "##### Exception in SFPositionRes constructor!";
    }

    try
    {
        fAtt = new SFAttenuation(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Error in SFPositionRes constructor! Problem with SFAttenuation!"
                  << std::endl;
        throw "##### Exception in SFPositionRes constructor!";
    }

    double              s      = SFTools::GetSigmaBL(fData->GetSiPM());
    std::vector<double> sigmas = {s, s};
    TString             cut    = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);
    fSpecAv                    = fData->GetCustomHistograms(SFSelectionType::kPEAverage, cut);

    fResultsPol3 = new SFResults(Form("PositionResResultsPol3_S%i", fSeriesNo));
    fResultsPol1 = new SFResults(Form("PositionResResultsPol1_S%i", fSeriesNo));
}
//------------------------------------------------------------------
SFPositionRes::~SFPositionRes()
{
    if (fAtt != nullptr) delete fAtt;
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
bool SFPositionRes::AnalyzePositionRes(void)
{

    std::cout << "\n\n----- Position Resolution Analysis" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << std::endl;

    int                 npointsMax      = fData->GetNpoints();
    std::vector<double> positions       = fData->GetPositions();
    std::vector<int>    measurementsIDs = fData->GetMeasurementsIDs();
    TString             collimator      = fData->GetCollimator();
    TString             testBench       = fData->GetTestBench();
    TString             sipm            = fData->GetSiPM();
    TString             desc            = fData->GetDescription();

    TGraphErrors *gPosRecoVsPosPol3Graph = new TGraphErrors(npointsMax);
    gPosRecoVsPosPol3Graph->SetMarkerStyle(4);
    gPosRecoVsPosPol3Graph->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoVsPosPol3Graph->GetYaxis()->SetTitle("reconstructed source position [mm]");
    gPosRecoVsPosPol3Graph->SetName(Form("PosRecoVsPosPol3_S%i", fSeriesNo));
    gPosRecoVsPosPol3Graph->SetTitle(Form("Reconstructed Source Position Pol3 S%i", fSeriesNo));

    TGraphErrors *gPosRecoVsPosPol1Graph = new TGraphErrors(npointsMax);
    gPosRecoVsPosPol1Graph->SetMarkerStyle(4);
    gPosRecoVsPosPol1Graph->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoVsPosPol1Graph->GetYaxis()->SetTitle("reconstructed source position [mm]");
    gPosRecoVsPosPol1Graph->SetName(Form("PosRecoVsPosPol1_S%i", fSeriesNo));
    gPosRecoVsPosPol1Graph->SetTitle(Form("Reconstructed Source Position Pol1 S%i", fSeriesNo));

    TGraphErrors *gPosResVsPosPol3Graph = new TGraphErrors(npointsMax);
    gPosResVsPosPol3Graph->SetMarkerStyle(4);
    gPosResVsPosPol3Graph->GetXaxis()->SetTitle("source position [mm]");
    gPosResVsPosPol3Graph->GetYaxis()->SetTitle("position resolution [mm]");
    gPosResVsPosPol3Graph->SetName(Form("MLRPosResVsPosPol3_S%i", fSeriesNo));
    gPosResVsPosPol3Graph->SetTitle(Form("Position Resolution Pol3 S%i", fSeriesNo));

    TGraphErrors *gPosResVsPosPol1Graph = new TGraphErrors(npointsMax);
    gPosResVsPosPol1Graph->SetMarkerStyle(4);
    gPosResVsPosPol1Graph->GetXaxis()->SetTitle("source position [mm]");
    gPosResVsPosPol1Graph->GetYaxis()->SetTitle("position resolution [mm]");
    gPosResVsPosPol1Graph->SetName(Form("MLRPosResVsPosPol1_S%i", fSeriesNo));
    gPosResVsPosPol1Graph->SetTitle(Form("Position Resolution Pol1 S%i", fSeriesNo));

    TGraphErrors *gPosRecoDiffPol1 = new TGraphErrors(npointsMax);
    gPosRecoDiffPol1->SetMarkerStyle(4);
    gPosRecoDiffPol1->SetTitle(Form("Reconstructed Position Difference (Experimental) S%i", fSeriesNo));
    gPosRecoDiffPol1->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoDiffPol1->GetYaxis()->SetTitle("P_{reco} - P_{real} [mm]");
    gPosRecoDiffPol1->SetName("gPosRecoDiffPol1");
    
    double mean, sigma;
    double meanErr;
    double MLR, pos_pol1, pos_pol3;
    double posResAv_p1    = 0;
    double posResAvErr_p1 = 0;
    double posResAv_p3    = 0;
    double posResAvErr_p3 = 0;
    double xmin, xmax;

    //-----
    fAtt->AttCombinedCh();
    std::vector<SFResults*> results_tmp = fAtt->GetResults();
    TGraphErrors* tmp = (TGraphErrors*)results_tmp[2]->GetObject(SFResultTypeObj::kAttGraph);

    double* x  = tmp->GetX();
    double* ex = tmp->GetEX();
    double* y  = tmp->GetY();
    double* ey = tmp->GetEY();

    fPosVsMLRGraph = new TGraphErrors(npointsMax, y, x, ey, ex);
    fPosVsMLRGraph->SetName("PosVsMLR");
    fPosVsMLRGraph->SetTitle(Form("Source Position vs. M_{LR} S%i", fSeriesNo));
    fPosVsMLRGraph->GetXaxis()->SetTitle("M_{LR}");
    fPosVsMLRGraph->GetYaxis()->SetTitle("source position [mm]");
    fPosVsMLRGraph->SetMarkerStyle(4);
    fPosVsMLRGraph->GetXaxis()->SetRangeUser(-1, 1);

    TF1* funPol3 = new TF1("funpol3", "pol3", -1, 1);
    funPol3->SetParLimits(3, 0, 100000);
    fPosVsMLRGraph->Fit(funPol3, "QR+");

    TF1* funPol1 = new TF1("funpol1", "pol1", -1, 1);
    fPosVsMLRGraph->Fit(funPol1, "QR+");
    
    //-----

    if (funPol3 == nullptr)
    {
        std::cerr << "##### Error in SFPositionRes::AnalyzePositionRes()" << std::endl;
        std::cerr << "Attenuation function was not found!" << std::endl;
        return false;
    }

    std::vector<TF1*>   funGausPol3;
    std::vector<TF1*>   funGausPol1;
    std::vector<double> FWHM;

    double BL_sigma_cut = SFTools::GetSigmaBL(sipm);
    
    TString hname = Form("Reconstructed Position Distribution (Pol3 & Summed) S%i", fSeriesNo);
    TH1D *hPosRecoPol3All = new TH1D("hPosRecoAll", hname, 300, -100, 200);
    hPosRecoPol3All->GetXaxis()->SetTitle("reconstructed position - source position [mm]");
    hPosRecoPol3All->GetYaxis()->SetTitle("counts");
    
    hname = Form("Reconstructed Position Distribution (Pol1 & Summed) S%i", fSeriesNo);
    TH1D *hPosRecoPol1All = new TH1D("hPosRecoAll", hname, 300, -100, 200);
    hPosRecoPol1All->GetXaxis()->SetTitle("reconstructed position - source position [mm]");
    hPosRecoPol1All->GetYaxis()->SetTitle("counts");
    
    for (int npoint = 0; npoint < npointsMax; npoint++)
    {

        std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;

        //----- geting tree
        SLoop* loop = fData->GetTree(measurementsIDs[npoint]);
        int        nloopMax = loop->getEntries();
        SCategory* tSig     = SCategoryManager::getCategory(SCategory::CatDDSamples);

        //----- setting energy cut
        //peakFinAv.push_back(new SFPeakFinder(fSpecAv[npoint], false));
        //peakFinAv[npoint]->FindPeakRange(xmin, xmax);
        auto peakFinAv = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(fSpecAv[npoint], false));
        peakFinAv->FindPeakRange(xmin, xmax);
        
        //----- setting histogram
        hname = Form("hPosRecoPol3_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fPosRecoPol3Dist.push_back(new TH1D(hname, hname, 300, -100, 200));
        fPosRecoPol3Dist[npoint]->SetTitle(Form("Reconstructed Position Pol3 S%i %1f mm", fSeriesNo, positions[npoint]));
        
        hname = Form("hPosRecoPol1_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fPosRecoPol1Dist.push_back(new TH1D(hname, hname, 300, -100, 200));
        fPosRecoPol1Dist[npoint]->SetTitle(Form("Reconstructed Position Pol1 S%i %1f mm", fSeriesNo, positions[npoint]));

        //----- filling histogram
        for (int nloop = 0; nloop < nloopMax; ++nloop)
        {
            loop->getEvent(nloop);
            size_t tentriesMax = tSig->getEntries();

            for (int tentries = 0; tentries < tentriesMax; ++tentries)
            {

                int              m, l, f;
                SDDSamples*      samples = (SDDSamples*)tSig->getObject(tentries);
                samples->getAddress(m, l, f);

                if (m == 0)
                {
                    double t0Ch0  = samples->getSignalL()->GetT0();
                    double t0Ch1  = samples->getSignalR()->GetT0();
                    double peCh0  = samples->getSignalL()->GetPE();
                    double peCh1  = samples->getSignalR()->GetPE();
                    double blCh0  = samples->getSignalL()->GetBLSigma();
                    double blCh1  = samples->getSignalR()->GetBLSigma();
                    double totCh0 = samples->getSignalL()->GetTOT();
                    double totCh1 = samples->getSignalR()->GetTOT();
                    double ampCh0 = samples->getSignalL()->GetAmplitude();
                    double ampCh1 = samples->getSignalR()->GetAmplitude();
                    bool   vetoCh0 = samples->getSignalL()->GetVeto();
                    bool   vetoCh1 = samples->getSignalR()->GetVeto();
                    
                    if (t0Ch0 > 0 && t0Ch1 > 0 && 
                        totCh0 > 0 && totCh1 > 0 &&
                        blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
                        ampCh0 < ampMax && ampCh1 < ampMax && 
                        sqrt(peCh0 * peCh1) > xmin && sqrt(peCh0 * peCh1) < xmax &&
                        vetoCh0 == 0 && vetoCh1 == 0)
                    {
                        MLR = log(sqrt(peCh1 / peCh0));
                        pos_pol3 = funPol3->Eval(MLR);
                        fPosRecoPol3Dist[npoint]->Fill(pos_pol3);
                        hPosRecoPol3All->Fill(pos_pol3 - positions[npoint]);
                        
                        pos_pol1 = funPol1->Eval(MLR);
                        fPosRecoPol1Dist[npoint]->Fill(pos_pol1);
                        hPosRecoPol1All->Fill(pos_pol1 - positions[npoint]);
                    }
                }
            }
        }

        delete loop;

        //----- fitting histogram and calculating position resolution /pol3/
        mean        = fPosRecoPol3Dist[npoint]->GetMean();
        sigma       = fPosRecoPol3Dist[npoint]->GetRMS();
        double xmin = fPosRecoPol3Dist[npoint]->GetBinCenter(2);
        funGausPol3.push_back(new TF1("funGausPol3", "gaus", xmin, 200));

        // if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
        fPosRecoPol3Dist[npoint]->Fit(funGausPol3[npoint], "QR");
        mean    = funGausPol3[npoint]->GetParameter(1);
        meanErr = funGausPol3[npoint]->GetParError(1);
        
        if (npoint == 0) FWHM.resize(2);
        FWHM[0] = 2 * sqrt(2 * log(2)) * funGausPol3[npoint]->GetParameter(2);
        FWHM[1] = 2 * sqrt(2 * log(2)) * funGausPol3[npoint]->GetParError(2);
        // }
        // else
        // {
        //     fPosRecoDist[npoint]->Fit(funGausPol3[npoint], "Q", "", mean-5*sigma, mean+5*sigma);
        //     mean     = funGausPol3[npoint]->GetParameter(1);
        //     meanErr  = funGausPol3[npoint]->GetParError(1);
        //     FWHM = SFTools::GetFWHM(fPosRecoDist[npoint]);
        // }

        gPosRecoVsPosPol3Graph->SetPoint(npoint, positions[npoint], mean);
        gPosRecoVsPosPol3Graph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                          meanErr);

        gPosResVsPosPol3Graph->SetPoint(npoint, positions[npoint], FWHM[0]);
        gPosResVsPosPol3Graph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                         FWHM[1]);

        gPosRecoDiffPol1->SetPoint(npoint, positions[npoint], mean - positions[npoint]);
        gPosRecoDiffPol1->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), FWHM[0]);
        
        posResAv_p3 += FWHM[0] * (1. / pow(FWHM[1], 2));
        posResAvErr_p3 += (1. / pow(FWHM[1], 2));
        
        //----- fitting histogram and calculating position resolution /pol1/
        
        mean  = fPosRecoPol1Dist[npoint]->GetMean();
        sigma = fPosRecoPol1Dist[npoint]->GetRMS();
        xmin  = fPosRecoPol1Dist[npoint]->GetBinCenter(2);
        funGausPol1.push_back(new TF1("funGaus", "gaus", xmin, 200));

        // if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
        fPosRecoPol1Dist[npoint]->Fit(funGausPol1[npoint], "QR");
        mean    = funGausPol1[npoint]->GetParameter(1);
        meanErr = funGausPol1[npoint]->GetParError(1);
        
        FWHM[0] = 2 * sqrt(2 * log(2)) * funGausPol1[npoint]->GetParameter(2);
        FWHM[1] = 2 * sqrt(2 * log(2)) * funGausPol1[npoint]->GetParError(2);
        // }
        // else
        // {
        //     fPosRecoDist[npoint]->Fit(funGausPol1[npoint], "Q", "", mean-5*sigma, mean+5*sigma);
        //     mean     = funGausPol1[npoint]->GetParameter(1);
        //     meanErr  = funGausPol1[npoint]->GetParError(1);
        //     FWHM = SFTools::GetFWHM(fPosRecoDist[npoint]);
        // }

        gPosRecoVsPosPol1Graph->SetPoint(npoint, positions[npoint], mean);
        gPosRecoVsPosPol1Graph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                          meanErr);

        gPosResVsPosPol1Graph->SetPoint(npoint, positions[npoint], FWHM[0]);
        gPosResVsPosPol1Graph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                         FWHM[1]);

        posResAv_p1 += FWHM[0] * (1. / pow(FWHM[1], 2));
        posResAvErr_p1 += (1. / pow(FWHM[1], 2));
    }

    posResAv_p3    = posResAv_p3 / posResAvErr_p3;
    posResAvErr_p3 = sqrt(1. / posResAvErr_p3);

    posResAv_p1    = posResAv_p1 / posResAvErr_p1;
    posResAvErr_p1 = sqrt(1. / posResAvErr_p1);

    //----- fitting summed histogram
    
    double gconst = hPosRecoPol3All->GetBinContent(hPosRecoPol3All->GetMaximumBin());
    double fit_min = hPosRecoPol3All->GetBinCenter(2);
    mean    = hPosRecoPol3All->GetMean();
    sigma   = hPosRecoPol3All->GetRMS();
    
    TF1 *fun_gauss_pol3 = new TF1("fun_gauss_pol3", "gaus", fit_min, 200);
    fun_gauss_pol3->SetParameters(gconst, mean, sigma);
    hPosRecoPol3All->Fit(fun_gauss_pol3, "QR");
    
    gconst  = hPosRecoPol1All->GetBinContent(hPosRecoPol1All->GetMaximumBin());
    fit_min = hPosRecoPol1All->GetBinCenter(2);
    mean    = hPosRecoPol1All->GetMean();
    sigma   = hPosRecoPol1All->GetRMS();
    
    TF1 *fun_gauss_pol1 = new TF1("fun_gauss_pol1", "gaus", fit_min, 200);
    fun_gauss_pol1->SetParameters(gconst, mean, sigma);
    hPosRecoPol1All->Fit(fun_gauss_pol1, "QR");
    
    //-----
    
    TF1* funpol1_p3 = new TF1("funpol1_p3", "pol1", 0, 100);
    gPosRecoVsPosPol3Graph->Fit(funpol1_p3, "Q");
    
    TF1* funpol1_p1 = new TF1("funpol1_p1", "pol1", 0, 100);
    gPosRecoVsPosPol1Graph->Fit(funpol1_p1, "Q");

    //----- residuals
    
    TGraphErrors *gResidualPol3 = new TGraphErrors(npointsMax);
    gResidualPol3->SetName(Form("PosRecoResidualsPol3_S%i", fSeriesNo));
    gResidualPol3->SetTitle(Form("Reconstructed Position Residuals Pol3 S%i", fSeriesNo));
    gResidualPol3->GetXaxis()->SetTitle("source position [mm]");
    gResidualPol3->GetYaxis()->SetTitle("residual [mm]");
    gResidualPol3->SetMarkerStyle(4);
    
    TGraphErrors *gResidualPol1 = new TGraphErrors(npointsMax);
    gResidualPol1->SetName(Form("PosRecoResidualsPol1_S%i", fSeriesNo));
    gResidualPol1->SetTitle(Form("Reconstructed Position Residuals Pol1 S%i", fSeriesNo));
    gResidualPol1->GetXaxis()->SetTitle("source position [mm]");
    gResidualPol1->GetYaxis()->SetTitle("residual [mm]");
    gResidualPol1->SetMarkerStyle(4);

    double res;
    double point_x, point_y;

    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        gPosRecoVsPosPol3Graph->GetPoint(npoint, point_x, point_y);
        res = point_y - positions[npoint];
        gResidualPol3->SetPoint(npoint, positions[npoint], res);
        
        gPosRecoVsPosPol1Graph->GetPoint(npoint, point_x, point_y);
        res = point_y - positions[npoint];
        gResidualPol1->SetPoint(npoint, positions[npoint], res);
    }

    //----- results
    
    std::cout << "\nAverage position resolution for this series is (pol3): ";
    std::cout << posResAv_p3 << " +/- " << posResAvErr_p3 << " mm\n\n" << std::endl;
    
    std::cout << "\nAverage position resolution for this series is (pol1): ";
    std::cout << posResAv_p1 << " +/- " << posResAvErr_p1 << " mm\n\n" << std::endl;

    fResultsPol3->AddResult(SFResultTypeNum::kPositionRes, posResAv_p3, posResAvErr_p3);
    fResultsPol3->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, gPosRecoVsPosPol3Graph);
    fResultsPol3->AddObject(SFResultTypeObj::kPosResVsPosGraph, gPosResVsPosPol3Graph);
    fResultsPol3->AddObject(SFResultTypeObj::kPosVsMLRGraph, fPosVsMLRGraph);
    fResultsPol3->AddObject(SFResultTypeObj::kResidualGraph, gResidualPol3);
    fResultsPol3->AddObject(SFResultTypeObj::kPositionAllHist, hPosRecoPol3All);
    
    fResultsPol1->AddResult(SFResultTypeNum::kPositionRes, posResAv_p1, posResAvErr_p1);
    fResultsPol1->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, gPosRecoVsPosPol1Graph);
    fResultsPol1->AddObject(SFResultTypeObj::kPosResVsPosGraph, gPosResVsPosPol1Graph);
    fResultsPol1->AddObject(SFResultTypeObj::kResidualGraph, gResidualPol1);
    fResultsPol1->AddObject(SFResultTypeObj::kPositionAllHist, hPosRecoPol1All);
    //fResultsPol1->AddObject(SFResultTypeObj::kPositionDiff, gPosRecoDiffPol1);

    return true;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFPositionRes::GetPositionRecoDist(TString type)
{

    std::vector<TH1D*> tmp;
    
    if (type == "pol3")
    {
        if (fPosRecoPol3Dist.empty())
        {
            std::cerr << "##### Error in SFPositionRes::GetPositionsRecoDist()! Empty vector!"
                      << std::endl;
            std::abort();
        }
        tmp = fPosRecoPol3Dist;
    }
    else if (type == "pol1")
    {
        if (fPosRecoPol1Dist.empty())
        {
            std::cerr << "##### Error in SFPositionRes::GetPositionsRecoDist()! Empty vector!"
                    << std::endl;
            std::abort();
        }
        tmp = fPosRecoPol1Dist;
    }
    else
    {
        std::cerr << "##### Error in SFPositionResolution::GetPositionDist()" << std::endl;
        std::cerr << "Incorrect type! possible options are: pol1 and pol3" << std::endl;
        std::abort();
    }

    return tmp;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFPositionRes::GetSpectra(void)
{

    if (fSpecAv.empty())
    {
        std::cerr << "##### Error in SFPositionRes::GetSpectra()!" << std::endl;
        std::cerr << "Empty vector!" << std::endl;
        std::abort();
    }

    return fSpecAv;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFPositionRes::GetRatios(void)
{
    if (fQRatios.empty())
    {
        std::cerr << "##### Error in SFPositionRes::GetRatios()!" << std::endl;
        std::cerr << "Empty vector!" << std::endl;
        std::abort();
    }

    return fQRatios;
}
//------------------------------------------------------------------
std::vector<SFResults*> SFPositionRes::GetResults(void)
{
    if (fResultsPol3 == nullptr ||
        fResultsPol1 == nullptr)
    {
        std::cerr << "##### Error in SFPositionRes::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }

    std::vector<SFResults*> results(2);
    results[0] = fResultsPol1;
    results[1] = fResultsPol3;

    return results;
}
//------------------------------------------------------------------
void SFPositionRes::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFPositionRes class object" << std::endl;
    std::cout << "Experimental series number: " << fSeriesNo << std::endl;
    std::cout << "\n-------------------------------------------" << std::endl;
}
//------------------------------------------------------------------
