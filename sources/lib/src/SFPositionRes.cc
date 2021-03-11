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
                                             fPosRecoVsPosGraph(nullptr),
                                             fPosResVsPosGraph(nullptr), fMLRvsPosGraph(nullptr), fResidualGraph(nullptr),
                                             fResults(nullptr)
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

    fResults = new SFResults(Form("PositionResResults_S%i", fSeriesNo));
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

    fPosRecoVsPosGraph = new TGraphErrors(npointsMax);
    fPosRecoVsPosGraph->SetMarkerStyle(4);
    fPosRecoVsPosGraph->GetXaxis()->SetTitle("source position [mm]");
    fPosRecoVsPosGraph->GetYaxis()->SetTitle("reconstructed source position [mm]");
    fPosRecoVsPosGraph->SetName(Form("PosRecoVsPos_S%i", fSeriesNo));
    fPosRecoVsPosGraph->SetTitle(Form("Reconstructed source position S%i", fSeriesNo));

    fPosResVsPosGraph = new TGraphErrors(npointsMax);
    fPosResVsPosGraph->SetMarkerStyle(4);
    fPosResVsPosGraph->GetXaxis()->SetTitle("source position [mm]");
    fPosResVsPosGraph->GetYaxis()->SetTitle("position resolution [mm]");
    fPosResVsPosGraph->SetName(Form("MLRPosResVsPos_S%i", fSeriesNo));
    fPosResVsPosGraph->SetTitle(Form("Position resolution S%i", fSeriesNo));

    double mean, sigma;
    double meanErr;
    double MLR, pos;
    double posResAv    = 0;
    double posResAvErr = 0;
    double xmin, xmax;

    //-----
    fAtt->AttCombinedCh();
    std::vector<SFResults*> results_tmp = fAtt->GetResults();
    TGraphErrors* tmp = (TGraphErrors*)results_tmp[2]->GetObject(SFResultTypeObj::kAttGraph);

    double* x  = tmp->GetX();
    double* ex = tmp->GetEX();
    double* y  = tmp->GetY();
    double* ey = tmp->GetEY();

    fMLRvsPosGraph = new TGraphErrors(npointsMax, y, x, ey, ex);
    fMLRvsPosGraph->SetName("PosVsMLR");
    fMLRvsPosGraph->SetTitle(Form("Source position vs. ln(M_{LR}) S%i", fSeriesNo));
    fMLRvsPosGraph->GetXaxis()->SetTitle("ln(M_{LR})");
    fMLRvsPosGraph->GetYaxis()->SetTitle("source position [mm]");
    fMLRvsPosGraph->SetMarkerStyle(4);
    fMLRvsPosGraph->GetXaxis()->SetRangeUser(-1, 1);

    TF1* funPol3 = new TF1("funpol3", "pol3", -1, 1);
    funPol3->SetParLimits(3, 0, 100000);
    fMLRvsPosGraph->Fit(funPol3, "QR");

    //-----

    if (funPol3 == nullptr)
    {
        std::cerr << "##### Error in SFPositionRes::AnalyzePositionRes()" << std::endl;
        std::cerr << "Attenuation function was not found!" << std::endl;
        return false;
    }

    std::vector<TF1*>   funGaus;
    std::vector<double> FWHM;

    double BL_sigma_cut = SFTools::GetSigmaBL(sipm);

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
        TString hname = Form("hPosReco_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fPosRecoDist.push_back(new TH1D(hname, hname, 500, -50, 150));

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
                    if (t0Ch0 > 0 && t0Ch1 > 0 && totCh0 > 0 && totCh1 > 0 &&
                        blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut && ampCh0 < ampMax &&
                        ampCh1 < ampMax && sqrt(peCh0 * peCh1) > xmin && sqrt(peCh0 * peCh1) < xmax)
                    {
                        MLR = log(sqrt(peCh1 / peCh0));
                        pos = funPol3->Eval(MLR);
                        fPosRecoDist[npoint]->Fill(pos);
                    }
                }
            }
        }

        delete loop;

        //----- fitting histogram and calculating position resolution
        mean        = fPosRecoDist[npoint]->GetMean();
        sigma       = fPosRecoDist[npoint]->GetRMS();
        double xmin = fPosRecoDist[npoint]->GetBinCenter(2);
        funGaus.push_back(new TF1("funGaus", "gaus", xmin, 200));

        // if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
        fPosRecoDist[npoint]->Fit(funGaus[npoint], "QR");
        mean    = funGaus[npoint]->GetParameter(1);
        meanErr = funGaus[npoint]->GetParError(1);
        if (npoint == 0) FWHM.resize(2);
        FWHM[0] = 2.35 * funGaus[npoint]->GetParameter(2);
        FWHM[1] = 2.35 * funGaus[npoint]->GetParError(2);
        // }
        // else
        // {
        //     fPosRecoDist[npoint]->Fit(funGaus[npoint], "Q", "", mean-5*sigma, mean+5*sigma);
        //     mean     = funGaus[npoint]->GetParameter(1);
        //     meanErr  = funGaus[npoint]->GetParError(1);
        //     FWHM = SFTools::GetFWHM(fPosRecoDist[npoint]);
        // }

        //fResults.fPosReco.push_back(mean);
        //fResults.fPosRecoErr.push_back(meanErr);

        fPosRecoVsPosGraph->SetPoint(npoint, positions[npoint], mean);
        fPosRecoVsPosGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                          meanErr);

        fPosResVsPosGraph->SetPoint(npoint, positions[npoint], FWHM[0]);
        fPosResVsPosGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench),
                                         FWHM[1]);

        posResAv += FWHM[0] * (1. / pow(FWHM[1], 2));
        posResAvErr += (1. / pow(FWHM[1], 2));
    }

    posResAv    = posResAv / posResAvErr;
    posResAvErr = sqrt(1. / posResAvErr);

    TF1* funpol1 = new TF1("funpol1", "pol1", 0, 100);
    fPosRecoVsPosGraph->Fit(funpol1, "Q");

    fResidualGraph = new TGraphErrors(npointsMax);
    fResidualGraph->SetName(Form("PosRecoResiduals_S%i", fSeriesNo));
    fResidualGraph->SetTitle(Form("Reconstructed position residuals S%i", fSeriesNo));
    fResidualGraph->GetXaxis()->SetTitle("source position [mm]");
    fResidualGraph->GetYaxis()->SetTitle("residual [mm]");
    fResidualGraph->SetMarkerStyle(4);

    double res;
    double point_x, point_y;

    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        fPosRecoVsPosGraph->GetPoint(npoint, point_x, point_y);
        res = point_y - positions[npoint];
        fResidualGraph->SetPoint(npoint, positions[npoint], res);
    }

    std::cout << "Average position resolution for this series is: ";
    std::cout << posResAv << " +/- " << posResAvErr << " mm\n\n" << std::endl;

    fResults->AddResult(SFResultTypeNum::kPositionRes, posResAv, posResAvErr);
    fResults->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, fPosRecoVsPosGraph);
    fResults->AddObject(SFResultTypeObj::kPosResVsPosGraph, fPosResVsPosGraph);
    fResults->AddObject(SFResultTypeObj::kMLRvsPosGraph, fMLRvsPosGraph);
    fResults->AddObject(SFResultTypeObj::kResidualGraph, fResidualGraph);

    return true;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFPositionRes::GetPositionRecoDist(void)
{

    if (fPosRecoDist.empty())
    {
        std::cerr << "##### Error in SFPositionRes::GetPositionsRecoDist()! Empty vector!"
                  << std::endl;
        std::abort();
    }

    return fPosRecoDist;
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
void SFPositionRes::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFPositionRes class object" << std::endl;
    std::cout << "Experimental series number: " << fSeriesNo << std::endl;
    std::cout << "\n-------------------------------------------" << std::endl;
}
//------------------------------------------------------------------
