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

const double ampMax = 660;


auto SFTools::calculatePositionResolutions(TH1* spectrum) -> std::tuple<float, std::vector<float>>
{

    double mean, sigma;
    double meanErr;
    double MLR, pos_pol1, pos_pol3;
    double posResAv_p1    = 0;
    double posResAvErr_p1 = 0;
    double posResAv_p3    = 0;
    double posResAvErr_p3 = 0;
    double xmin, xmax;

    //-----
    TF1* funPol3 = new TF1("funpol3", "pol3", -1, 1);
    funPol3->SetParLimits(3, 0, 100000);
    fPosVsMLRGraph->Fit(funPol3, "QR+");

    TF1* funPol1 = new TF1("funpol1", "pol1", -1, 1);
    fPosVsMLRGraph->Fit(funPol1, "QR+");
    
    //-----

    if (mlr == nullptr)
    {
        std::cerr << "##### Error in SFTools::" << __func__ << std::endl;
        std::cerr << "Attenuation function was not found!" << std::endl;
        std::abort();
    }

    TF1*   funGaus;
    double FWHM;

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
        
        auto mean  = spectrum->GetMean();
        auto sigma = spectrum->GetRMS();
        auto xmin  = spectrum->GetBinCenter(2);
        auto fit_fun = new TF1("funGaus", "gaus", xmin, 200);
        std::vector<double> FWHM;

        spectrum->Fit(fit)fun, "QR");
        mean    = spectrm->GetParameter(1);
        auto meanErr = spectrum->GetParError(1);
        
        FWHM[0] = 2 * sqrt(2 * log(2)) * funGausPol1[npoint]->GetParameter(2);
        FWHM[1] = 2 * sqrt(2 * log(2)) * funGausPol1[npoint]->GetParError(2);

    return {mean, FWHM};
}

//------------------------------------------------------------------
double SFPositionRes::ReconstructPosition(SDDSamples* samples, TF1* fun_mlr, 
                                        float BL_sigma_cut, float qmin, float qmax)
                                        -> std::optional<double>
{
    double T0_l   =  samples->getSignalL()->GetT0();
    double T0_r   = samples->getSignalR()->GetT0();
    double PE_l   = samples->getSignalL()->GetPE();
    double PE_r   = samples->getSignalR()->GetPE();
    double BL_l   = samples->getSignalL()->GetBLSigma();
    double BL_r   = samples->getSignalR()->GetBLSigma();
    double TOT_l  = samples->getSignalL()->GetTOT();
    double TOT_r  = samples->getSignalR()->GetTOT();
    double amp_l  = samples->getSignalL()->GetAmplitude();
    double amp_r  = samples->getSignalR()->GetAmplitude();
    bool   veto_l = samples->getSignalL()->GetVeto();
    bool   veto_r = samples->getSignalR()->GetVeto();

    if (T0_l > 0 && T0_r > 0 &&
        TOT_l > 0 && TOT_r > 0 &&
        BL_l < BL_sigma_cut && BL_r < BL_sigma_cut &&
        amp_l < ampMax && amp_r < ampMax && 
        sqrt(PE_l * PE_r) > qmin && sqrt(PE_l * PE_r) < qmax &&
        veto_l == 0 && veto_r == 0)
    {
        double MLR = log(sqrt(PE_r / PE_l));
        return fun_mlr->Eval(MLR);
    }

    return std::nullopt;
}

//------------------------------------------------------------------
SFResults* SFPositionRes::ReconstructPositionDist(SFChAddr addr, SLoop *loop, TH1D* spectrum_av,
                                                  TF1* fun_mlr, double BL_sigma_cut)
{
    
    int nloopMax = loop->getEntries();
    auto tSig = std::unique_ptr<SCategory>(SCategoryManager::getCategory(SCategory::CatDDSamples));
    
    //----- setting energy cut
    double xmin = 0;
    double xmax = 0;
    
    auto peakFinAv = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(spectrum_av, false));
    peakFinAv->FindPeakRange(xmin, xmax);
    
    TH1D* hist = new TH1D("htemp", "htemp", 300, -100, 200);
    hist->GetXaxis()->SetTitle("reconstructed position [mm]");
    hist->GetYaxis()->SetTitle("counts");
    
    for (int nloop = 0; nloop < nloopMax; ++nloop)
    {
        loop->getEvent(nloop);
        size_t tentriesMax = tSig->getEntries();

        for (int tentries = 0; tentries < tentriesMax; ++tentries)
        {
            int m, l, f;
            SDDSamples* samples = (SDDSamples*)tSig->getObject(tentries);
            samples->getAddress(m, l, f);

                if (addr.fModule == m &&
                    addr.fLayer == l &&
                    addr.fFiber == f)
                {
                    hist->Fill(ReconstructPosition(samples, fun_mlr, BL_sigma_cut, xmin, xmax));
                }
            }
        }
    
    double fit_min = hist->GetBinCenter(2);
    TF1* fun_gauss = new TF1("fun_gaus", "gaus", fit_min, 200);
    fun_gauss->SetParameters(hist->GetBinContent(hist->GetMaximumBin()),
                             hist->GetMean(), hist->GetRMS());
    
    hist->Fit(fun_gauss, "QR");

    double sigma_to_fwhm = 2 * sqrt(2 * log(2));
    double pos_res = sigma_to_fwhm * fun_gauss->GetParameter(2);
    double pos_res_err = sigma_to_fwhm * fun_gauss->GetParError(2);
    double pos_reco = fun_gauss->GetParameter(1);
    double pos_reco_err = fun_gauss->GetParameter(1);
    
    SFResults* results = new SFResults("PositionReconstructionDistSinglePos");
    results->AddObject(SFResultTypeObj::kPositionDist, hist);
    results->AddResult(SFResultTypeNum::kPositionRes, pos_res, pos_res_err);
    results->AddResult(SFResultTypeNum::kPositionReco, pos_reco, pos_reco_err);
    
    return results;
}
//------------------------------------------------------------------
auto SFPositionRes::ReconstructPositionDistAll(SFChAddr addr, SLoop* loop, TH1D* spectrum_av, 
                                               TF1* fun_mlr, std::vector<double> positions,
                                               double BL_sigma_cut, double pos_uncert, TString suffix) 
                                               -> std::tuple<SFResults*,std::vector<TH1D*>>
{
    
    if (loop == nullptr ||
        spectrum_av == nullptr ||
        fun_mlr == nullptr)
    {
        std::cerr << "##### Error in SFPositionRes:: " << __func__ << std::endl;
        std::cerr << "One of the required pointers is invalid! Please check!" << std::endl;
        std::abort();
    }
    
    int npointsMax = positions.size();
    
    TGraphErrors *gPosRecoVsPos = new TGraphErrors(npointsMax);
    gPosRecoVsPos->SetMarkerStyle(4);
    gPosRecoVsPos->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoVsPos->GetYaxis()->SetTitle("reconstructed source position [mm]");
    gPosRecoVsPos->SetName(Form("PosRecoVsPos_%s", suffix.Data()));
    gPosRecoVsPos->SetTitle(Form("Reconstructed Source Position (%s)", suffix.Data()));
    
    TGraphErrors *gPosResVsPos = new TGraphErrors(npointsMax);
    gPosResVsPos->SetMarkerStyle(4);
    gPosResVsPos->GetXaxis()->SetTitle("source position [mm]");
    gPosResVsPos->GetYaxis()->SetTitle("position resolution [mm]");
    gPosResVsPos->SetName(Form("PosResVsPos_%s", suffix.Data()));
    gPosResVsPos->SetTitle(Form("Position Resolution (%s)", suffix.Data()));
    
    TGraphErrors *gPosRecoDiff = new TGraphErrors(npointsMax);
    gPosRecoDiff->SetMarkerStyle(4);
    gPosRecoDiff->SetTitle(Form("Reconstructed Position Difference (%s)", suffix.Data()));
    gPosRecoDiff->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoDiff->GetYaxis()->SetTitle("P_{reco} - P_{real} [mm]");
    gPosRecoDiff->SetName("gPosRecoDiff");
    
    std::vector<TH1D*> pos_distributions;
    
    double posResAv = 0;
    double posResAvErr = 0;
    
    for (int i = 0; i < npointsMax; i++)
    {
        auto results_single_pos = std::unique_ptr<SFResults>(ReconstructPositionDist(addr, loop, 
                                                             spectrum_av,fun_mlr, BL_sigma_cut));
        pos_distributions.push_back((TH1D*)results_single_pos->GetObject(SFResultTypeObj::kPositionDist));
        
        double pos_res = results_single_pos->GetValue(SFResultTypeNum::kPositionRes);
        double pos_res_err = results_single_pos->GetUncertainty(SFResultTypeNum::kPositionRes);
        
        double pos_reco = results_single_pos->GetValue(SFResultTypeNum::kPositionReco);
        double pos_reco_err = results_single_pos->GetUncertainty(SFResultTypeNum::kPositionReco);
        
        gPosRecoVsPos->SetPoint(i, positions[i], pos_reco);
        gPosRecoVsPos->SetPointError(i, pos_uncert, pos_reco_err);
        
        gPosResVsPos->SetPoint(i, positions[i], pos_res);
        gPosResVsPos->SetPointError(i, pos_uncert, pos_res_err);
        
        gPosRecoDiff->SetPoint(i, positions[i], (pos_reco - positions[i]));
        
        posResAv += pos_res * (1. / pow(pos_res_err, 2));
        posResAvErr += (1. / pow(pos_res_err, 2));
        
    }
    
    posResAv    = posResAv / posResAvErr;
    posResAvErr = sqrt(1. / posResAvErr);
    
    SFResults* results = new SFResults(Form("PositionResolution_%s", suffix.Data()));
    results->AddResult(SFResultTypeNum::kPositionRes, posResAv, posResAvErr);
    results->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, gPosRecoVsPos);
    results->AddObject(SFResultTypeObj::kPosResVsPosGraph, gPosResVsPos);
    results->AddObject(SFResultTypeObj::kResidualGraph, gPosRecoDiff);
    
    return {results, pos_distributions};
}



bool SFPositionRes::AnalyzePositionRes(void)
{

    int                 npointsMax      = fData->GetNpoints();
    std::vector<double> positions       = fData->GetPositions();
    std::vector<int>    measurementsIDs = fData->GetMeasurementsIDs();
    TString             collimator      = fData->GetCollimator();
    TString             testBench       = fData->GetTestBench();
    TString             sipm            = fData->GetSiPM();
    TString             desc            = fData->GetDescription();

//     TGraphErrors *gPosRecoVsPosPol3Graph = new TGraphErrors(npointsMax);
//     gPosRecoVsPosPol3Graph->SetMarkerStyle(4);
//     gPosRecoVsPosPol3Graph->GetXaxis()->SetTitle("source position [mm]");
//     gPosRecoVsPosPol3Graph->GetYaxis()->SetTitle("reconstructed source position [mm]");
//     gPosRecoVsPosPol3Graph->SetName(Form("PosRecoVsPosPol3_S%i", fSeriesNo));
//     gPosRecoVsPosPol3Graph->SetTitle(Form("Reconstructed Source Position Pol3 S%i", fSeriesNo));

//     TGraphErrors *gPosRecoVsPosPol1Graph = new TGraphErrors(npointsMax);
//     gPosRecoVsPosPol1Graph->SetMarkerStyle(4);
//     gPosRecoVsPosPol1Graph->GetXaxis()->SetTitle("source position [mm]");
//     gPosRecoVsPosPol1Graph->GetYaxis()->SetTitle("reconstructed source position [mm]");
//     gPosRecoVsPosPol1Graph->SetName(Form("PosRecoVsPosPol1_S%i", fSeriesNo));
//     gPosRecoVsPosPol1Graph->SetTitle(Form("Reconstructed Source Position Pol1 S%i", fSeriesNo));

//     TGraphErrors *gPosResVsPosPol3Graph = new TGraphErrors(npointsMax);
//     gPosResVsPosPol3Graph->SetMarkerStyle(4);
//     gPosResVsPosPol3Graph->GetXaxis()->SetTitle("source position [mm]");
//     gPosResVsPosPol3Graph->GetYaxis()->SetTitle("position resolution [mm]");
//     gPosResVsPosPol3Graph->SetName(Form("MLRPosResVsPosPol3_S%i", fSeriesNo));
//     gPosResVsPosPol3Graph->SetTitle(Form("Position Resolution Pol3 S%i", fSeriesNo));

//     TGraphErrors *gPosResVsPosPol1Graph = new TGraphErrors(npointsMax);
//     gPosResVsPosPol1Graph->SetMarkerStyle(4);
//     gPosResVsPosPol1Graph->GetXaxis()->SetTitle("source position [mm]");
//     gPosResVsPosPol1Graph->GetYaxis()->SetTitle("position resolution [mm]");
//     gPosResVsPosPol1Graph->SetName(Form("MLRPosResVsPosPol1_S%i", fSeriesNo));
//     gPosResVsPosPol1Graph->SetTitle(Form("Position Resolution Pol1 S%i", fSeriesNo));

//     TGraphErrors *gPosRecoDiffPol1 = new TGraphErrors(npointsMax);
//     gPosRecoDiffPol1->SetMarkerStyle(4);
//     gPosRecoDiffPol1->SetTitle(Form("Reconstructed Position Difference (Experimental) S%i", fSeriesNo));
//     gPosRecoDiffPol1->GetXaxis()->SetTitle("source position [mm]");
//     gPosRecoDiffPol1->GetYaxis()->SetTitle("P_{reco} - P_{real} [mm]");
//     gPosRecoDiffPol1->SetName("gPosRecoDiffPol1");
    
    double mean, sigma;
    double meanErr;
    double MLR, pos_pol1, pos_pol3;
    double posResAv_p1    = 0;
    double posResAvErr_p1 = 0;
    double posResAv_p3    = 0;
    double posResAvErr_p3 = 0;
    double xmin, xmax;

    //-----
//     fAtt->AttCombinedCh();
//     std::vector<SFResults*> results_tmp = fAtt->GetResults();
//     TGraphErrors* tmp = (TGraphErrors*)results_tmp[2]->GetObject(SFResultTypeObj::kAttGraph);
// 
//     double* x  = tmp->GetX();
//     double* ex = tmp->GetEX();
//     double* y  = tmp->GetY();
//     double* ey = tmp->GetEY();
// 
//     fPosVsMLRGraph = new TGraphErrors(npointsMax, y, x, ey, ex);
//     fPosVsMLRGraph->SetName("PosVsMLR");
//     fPosVsMLRGraph->SetTitle(Form("Source Position vs. M_{LR} S%i", fSeriesNo));
//     fPosVsMLRGraph->GetXaxis()->SetTitle("M_{LR}");
//     fPosVsMLRGraph->GetYaxis()->SetTitle("source position [mm]");
//     fPosVsMLRGraph->SetMarkerStyle(4);
//     fPosVsMLRGraph->GetXaxis()->SetRangeUser(-1, 1);
// 
//     TF1* funPol3 = new TF1("funpol3", "pol3", -1, 1);
//     funPol3->SetParLimits(3, 0, 100000);
//     fPosVsMLRGraph->Fit(funPol3, "QR+");
// 
//     TF1* funPol1 = new TF1("funpol1", "pol1", -1, 1);
//     fPosVsMLRGraph->Fit(funPol1, "QR+");
    
//     //-----
// 
//     if (funPol3 == nullptr)
//     {
//         std::cerr << "##### Error in SFPositionRes::AnalyzePositionRes()" << std::endl;
//         std::cerr << "Attenuation function was not found!" << std::endl;
//         return false;
//     }

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
        
        //----- setting histogram
        hname = Form("hPosRecoPol3_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fPosRecoPol3Dist.push_back(new TH1D(hname, hname, 300, -100, 200));
        fPosRecoPol3Dist[npoint]->SetTitle(Form("Reconstructed Position Pol3 S%i %1f mm", fSeriesNo, positions[npoint]));
        
        hname = Form("hPosRecoPol1_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fPosRecoPol1Dist.push_back(new TH1D(hname, hname, 300, -100, 200));
        fPosRecoPol1Dist[npoint]->SetTitle(Form("Reconstructed Position Pol1 S%i %1f mm", fSeriesNo, positions[npoint]));

        //----- filling histogram
//                if (t0Ch0 > 0 && t0Ch1 > 0 && 
//                         totCh0 > 0 && totCh1 > 0 &&
//                         blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
//                         ampCh0 < ampMax && ampCh1 < ampMax && 
//                         sqrt(peCh0 * peCh1) > xmin && sqrt(peCh0 * peCh1) < xmax &&
//                         vetoCh0 == 0 && vetoCh1 == 0)
//                     {
//                         MLR = log(sqrt(peCh1 / peCh0));
//                         pos_pol3 = funPol3->Eval(MLR);
//                         fPosRecoPol3Dist[npoint]->Fill(pos_pol3);
//                         hPosRecoPol3All->Fill(pos_pol3 - positions[npoint]);
//                         
//                         pos_pol1 = funPol1->Eval(MLR);
//                         fPosRecoPol1Dist[npoint]->Fill(pos_pol1);
//                         hPosRecoPol1All->Fill(pos_pol1 - positions[npoint]);
//                     }
        

        //----- fitting histogram and calculating position resolution /pol3/
//         mean        = fPosRecoPol3Dist[npoint]->GetMean();
//         sigma       = fPosRecoPol3Dist[npoint]->GetRMS();
//         double xmin = fPosRecoPol3Dist[npoint]->GetBinCenter(2);
//         funGausPol3.push_back(new TF1("funGausPol3", "gaus", xmin, 200));

        // if(collimator.Contains("Electronic") && sipm.Contains("SensL")){
        fPosRecoPol3Dist[npoint]->Fit(funGausPol3[npoint], "QR");
        mean    = funGausPol3[npoint]->GetParameter(1);
        meanErr = funGausPol3[npoint]->GetParError(1);
        
//         if (npoint == 0) FWHM.resize(2);
//         FWHM[0] = 2 * sqrt(2 * log(2)) * funGausPol3[npoint]->GetParameter(2);
//         FWHM[1] = 2 * sqrt(2 * log(2)) * funGausPol3[npoint]->GetParError(2);
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
