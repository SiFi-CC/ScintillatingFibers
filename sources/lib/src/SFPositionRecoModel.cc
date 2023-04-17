// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFPositionRecoModel.cc        *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFPositionRecoModel.hh"

const double ampMax = 660.;

//------------------------------------------------------------------
auto SFPositionRecoModel::CalculateMLR(TGraphErrors *gAttCorrL, TGraphErrors *gAttCorrR)
                                       -> std::tuple<TGraphErrors*, TGraphErrors*>
{    
    if (gAttCorrL == nullptr || gAttCorrR == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco::" << __func__ << std::endl;
        std::cerr << "Corrected attenuation graphs do not exist!" << std::endl;
        std::abort();
    }
    
    int npoints = gAttCorrL->GetN();
    double* positions = gAttCorrL->GetX();
    double* positions_err = gAttCorrL->GetEX();
    
    TGraphErrors *gMLRCorr = new TGraphErrors(npoints);
    gMLRCorr->SetName("gMLRCorr");
    gMLRCorr->SetTitle("M_{LR} vs. Position (Corrected)");
    gMLRCorr->GetXaxis()->SetTitle("source position [mm]");
    gMLRCorr->GetYaxis()->SetTitle("M_{LR} corrected");
    gMLRCorr->SetMarkerStyle(8);

    TGraphErrors *gRevMLRCorr = new TGraphErrors(npoints);
    gRevMLRCorr->SetName("gRevMLRCorr");
    gRevMLRCorr->SetTitle("Position vs. M_{LR} (Corrected)");
    gRevMLRCorr->GetXaxis()->SetTitle("M_{LR}");
    gRevMLRCorr->GetYaxis()->SetTitle("source position [mm]");
    gRevMLRCorr->SetMarkerStyle(8);
    
    for (int i = 0; i < npoints; i++)
    {
        double vl    = gAttCorrL->GetPointY(i);
        double vl_err = gAttCorrL->GetErrorY(i);
        double vr    = gAttCorrR->GetPointY(i);
        double vr_err = gAttCorrR->GetErrorY(i);
        double MLR     = log(sqrt(vr / vl));
        double MLRerr  = sqrt(pow(vr_err / (2 * vr), 2) +
                              pow(- vl_err / (2 * vl), 2));
        
        gMLRCorr->SetPoint(i, positions[i], MLR);
        gMLRCorr->SetPointError(i, positions_err[i], MLRerr);
        
        gRevMLRCorr->SetPoint(i, MLR, positions[i]);
        gRevMLRCorr->SetPointError(i, MLRerr, positions_err[i]);
    }
    
    TF1 *fpol1 = new TF1("fpol1", "pol1", positions[0], positions[npoints-1]); 
    TFitResultPtr ptr = gMLRCorr->Fit(fpol1, "SQR+");
    
    TF1 *fpol1_rev = new TF1("fpol1", "pol1", -1, 1);
    TFitResultPtr ptr_rev = gRevMLRCorr->Fit(fpol1_rev, "SQR+");
    
    return {gMLRCorr, gRevMLRCorr};
}
//------------------------------------------------------------------
auto SFPositionRecoModel::CalculateRecoCoefficients(TGraphErrors *gMLRCorr, double fiberLen)
                                                    -> std::tuple<TGraphErrors*, double, double>
{
    if (gMLRCorr == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco::" << __func__ << std::endl;
        std::cerr << "MLR graph do not exist!" << std::endl;
        std::abort();
    }
    
    int npoints = gMLRCorr->GetN();
    double* positions = gMLRCorr->GetX();
    double* positions_err = gMLRCorr->GetEX();
    
    double B = fiberLen / 2.;
    
    int n = 0;
    
    for (int i = 0; i < npoints; i++)
    {
        if (fabs(B - positions[i]) < 1E-2)
            continue;
        n++;
    }
    
    TGraphErrors* gACoeff = new TGraphErrors(n);
    gACoeff->SetName("gACoeff");
    gACoeff->SetTitle("A Coefficient vs. Position");
    gACoeff->GetXaxis()->SetTitle("source position [mm]");
    gACoeff->GetYaxis()->SetTitle("A [mm]");
    gACoeff->SetMarkerStyle(8);
    
    int j = -1;
    
    for (int i = 0; i < npoints; i++)
    {
        double val   = fMLRCorrGraph->GetPointY(i);
        double val_e = fMLRCorrGraph->GetErrorY(i);
        if (fabs(B - positions[i]) < 1E-2)
            continue;
        j++;
        double A    = (positions[i] - B) / (val);
        double Aerr = sqrt(pow((positions_err[i] / val), 2) +
                           pow((B - positions[i]) * val_e / (pow(val, 2)), 2));
        gACoeff->SetPoint(j, positions[i], A);
        gACoeff->SetPointError(j, positions_err[i], Aerr);
    }
    
    TF1* fpol0 = new TF1("fpol0", "pol0", positions[0], positions[npoints-1]);
    fAGraph->Fit(fpol0, "SQR+");
    
    std::cout << "A = " << fpol0->GetParameter(0) << " +/- " << fpol0->GetParError(0) << std::endl;
    std::cout << "B = " << B << std::endl;
    
    return {gACoeff, A, B};
}
//------------------------------------------------------------------
auto SFPositionRecoModel::ReconstructPosition(SDDSamples *samples, SFResults *results, TMatrixD covMatrix,
                                              TF1 *fPrReco, TF1 *fPlReco, std::vector<double> parsForErrors,
                                              float BL_sigma_cut, float xmin, float xmax)
                                              -> std::optional<std::tuple<double, double>>
{
    double T0_l   = samples->getSignalL()->GetT0();
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

    double A     = results->GetValue(SFResultTypeNum::kACoeff);
    double A_err = results->GetUncertainty(SFResultTypeNum::kACoeff);
    double B     = results->GetValue(SFResultTypeNum::kBCoeff);
    double B_err = results->GetUncertainty(SFResultTypeNum::kBCoeff);
    
    std::vector<double> vec(2);
    
    if (T0_l > 0 && T0_r > 0 &&
        TOT_l > 0 && TOT_r > 0 &&
        BL_l < BL_sigma_cut && BL_r < BL_sigma_cut &&
        amp_l < ampMax && amp_r < ampMax &&
        sqrt(PE_l * PE_r) > xmin &&
        sqrt(PE_l * PE_r) < xmax &&
        veto_l == 0 && veto_r == 0)
    {
        
    parsForErrors[5] = PE_l;
    parsForErrors[6] = PE_r;
                        
    double MLR = log(sqrt(fPrReco->Eval(PE_r, PE_l) /
                          fPlReco->Eval(PE_r, PE_l)));
    double pos = A * MLR + B;
    double Pr_err = SFAttenuationModel::CalculateUncertainty(parsForErrors, covMatrix, "R");
    double Pl_err = SFAttenuationModel::CalculateUncertainty(parsForErrors, covMatrix, "L");
    double pos_err = sqrt(pow(MLR * A_err, 2) +
                     pow(A / (2 * fPrReco->Eval(PE_r, PE_l)) * Pr_err, 2) +
                     pow(-A / (2 * fPlReco->Eval(PE_r, PE_l)) * Pl_err, 2) + 
                     pow(B_err, 2));

    return {pos, pos_err};
    }
    
    return std::nullopt;
}
//------------------------------------------------------------------
SFResults* SFPositionRecoModel::ReconstructPositionDist(SFChAddr addr, SLoop* loop,
                                                        TH1D* spectrum_av, TH1D *spectrum_l,
                                                        TH1D *spectrum_r, SFResults* results,
                                                        TMatrixD covMatrix, TF1* fPrReco,
                                                        std::vector<double> parsForErrors,
                                                        TF1* fPlReco, TString path,
                                                        double BL_sigma_cut, TString collimator)
{
    int nloopMax = loop->getEntries();
    auto tSig = std::unique_ptr<SCategory>(SCategoryManager::getCategory(SCategory::CatDDSamples));
    
    //----- setting energy cut
    auto peak_range_ave = SFPeakFinder::FindPeakRange(spectrum_av, path, collimator, 0, 0));
    double xmin_ave = std::get<0>(peak_range_ave);
    double xmax_ave = std::get<1>(peak_range_ave);
    
    TString hname = "hPosRecoDistCorr"
    TH1D *hPosReco = new TH1D(hname, hnam3, 300, -100, 200);
    hPosReco->GetXaxis()->SetTitle("reconstructed position [mm]");
    hPosReco->GetYaxis()->SetTitle("counts");
    
    hname = "hPosRecoCorrUncert";
    TH1D *hPosRecoUncert = new TH1D(hname, hname, 200, 0, 50);
    hPosRecoUncert->GetXaxis()->SetTitle("uncertainty of reconstructed position [mm]");
    hPosRecoUncert->GetYaxis()->SetTitle("counts");
        
    auto peak_fit_l = SFPeakFinder::FindPeakFit(spectrum_l, path, 0, 0);
    auto peak_fit_r = SFPeakFinder::FindPeakFit(spectrum_r, path, 0, 0);
        
    parsForErrors[7] = peak_fit_l->GetValue(SFResultTypeNum::kPeakSigma);    
    parsForErrors[8] = peak_fit_r->GetValue(SFResultTypeNum::kPeakSigma);
    
    for (int nloop = 0; nloop < nloopMax; nloop++)
    {
        loop->getEvent(nloop);
        size_t tentriesMax = tSig->getEntries();
            
        for (int tentries = 0; tentries < tentriesMax; ++tentries)
        {
            int m, l, f;
            auto samples = std::unique_ptr<SDDSamples>((SDDSamples*)tSig->getObject(tentries));
            samples->getAddress(m, l, f);
                
            if (addr.fModule == m &&
                addr.fLayer == l &&
                addr.fFiber == f)
            {
                auto pos_reco = ReconstructPosition(samples, results, covMatrix, fPrReco, fPlReco, 
                                                    parsForErrors, BL_sigma_cut, xmin_ave, xmax_ave);
                double pos = std::get<0>(pos_reco);
                double pos_err = std::get<1>(pos_reco);
                hPosReco->Fill(pos);
                hPosRecoUncert->Fill(pos_err);
            }
        }
    }
    
    double fit_min = hPosReco->GetBinCenter(2);
    TF1* fun_gauss = new TF1("fun_gaus", "gaus", fit_min, 200);
    fun_gauss->SetParameters(hPosReco->GetBinContent(hPosReco->GetMaximumBin()),
                             hPosReco->GetMean(), hPosReco->GetRMS());
    
    hPosReco->Fit(fun_gauss, "QR");

    double sigma_to_fwhm = 2 * sqrt(2 * log(2));
    double pos_res = sigma_to_fwhm * fun_gauss->GetParameter(2);
    double pos_res_err = sigma_to_fwhm * fun_gauss->GetParError(2);
    double pos_reco = fun_gauss->GetParameter(1);
    double pos_reco_err = fun_gauss->GetParameter(1);
    
    SFResults *single_position_results = new SFResults("CorrPositionReconstructionDistSinglePos");
    single_position_results->AddObject(SFResultTypeObj::kPositionDist, hPosReco);
    single_position_results->AddObject(SFResultTypeObj::kPositionDistUncert, hPosRecoUncert);
    single_position_results->AddResult(SFResultTypeNum::kPositionRes, pos_res, pos_res_err);
    single_position_results->AddResult(SFResultTypeNum::kPositionReco, pos_reco, pos_reco_err);
    
    return single_position_results;    
}
//------------------------------------------------------------------
auto SFPositionRecoModel::ReconstructPositionDistAll(SFChAddr addr, 
                                                     std::vector<SLoop*> loop,
                                                     std::vector<TH1D*> spectra_av,
                                                     std::vector<TH1D*> spectra_l,
                                                     std::vector<TH1D*> spectra_r,
                                                     std::vector<double> positions,
                                                     std::vector<TString> path,
                                                     TMatrixD covMatrix,
                                                     TF1 *fPrReco,
                                                     TF1 *fPlReco,
                                                     TGraphErrors *gAttCorrL,
                                                     TGraphErrors *gAttCorrR,
                                                     double BL_sigma_cut,
                                                     double pos_uncert,
                                                     double fiberLen,
                                                     TString collimator,
                                                     TString suffix)
                          -> std::tuple <SFResults*, std::vector<TH1D*>, std::vector <TH1D*>>
{
    
    int npointsMax = positions->GetN();
    
    //----- setting graphs
    TGraphErrors *gPosReco = new TGraphErrors(npointsMax);
    gPosReco->SetName("gPosReco");
    gPosReco->SetTitle("Reconstructed Source Position (Corrected)");
    gPosReco->GetXaxis()->SetTitle("source position [mm]");
    gPosReco->GetYaxis()->SetTitle("reconstructed position [mm]");
    gPosReco->SetMarkerStyle(8);
    
    TGraphErrors *gPosResolution = new TGraphErrors(npointsMax);
    gPosResolution->SetName("gPosResolution");
    gPosResolution->SetTitle("Position Resolution (Corrected) ");
    gPosResolution->GetXaxis()->SetTitle("source position [mm]");
    gPosResolution->GetYaxis()->SetTitle("position resolution [mm]");
    gPosResolution->SetMarkerStyle(8);
    
    TGraphErrors *gPosResiduals = new TGraphErrors(npointsMax);
    gPosResiduals->SetName("gPosResiduals");
    gPosResiduals->SetTitle("Reconstructed Position (Corrected) Residuals");
    gPosResiduals->GetXaxis()->SetTitle("source position [mm]");
    gPosResiduals->GetYaxis()->SetTitle("residual [mm]");
    gPosResiduals->SetMarkerStyle(8);
    
    TGraphErrors *gPosRecoDiff = new TGraphErrors(npointsMax);
    gPosRecoDiff->SetName("fPosRecoDiff");
    gPosRecoDiff->SetTitle("Reconstructed Position Difference (Corrected)");
    gPosRecoDiff->GetXaxis()->SetTitle("source position [mm]");
    gPosRecoDiff->GetYaxis()->SetTitle("P_{reco} - P_{real} [mm]");
    gPosRecoDiff->SetMarkerStyle(8);
    //-----
   
    //----- calculating reconstruction coefficients
    auto mlr = CalculateMLR(gAttCorrL, gAttCorrR);
    TGraphErrors *gMLRCorr = std::get<0>(mlr);
    TGraphErrors *gRevMLRCorr = std::get<1>(mlr);
    TF1 *fun_rev_mlr = gRevMLRCorr.GetFunction("fpol1");
    
    auto coeff = CalculateRecoCoefficients(gMLRCorr, fiberLen);
    
    double A = fun_rev_mlr.GetParameter(1);
    double A_err = fun_rev_mlr.GetParError(1);
    double B = std::get<2>(coeff);
    double B_err = 0.;
    TGraphErrors *gACoeff = std::get<0>(coeff);
    
    //----- position reconstruction event by event start
    
    double posResSum    = 0.;
    double posResSumErr = 0.;
    
    std::vector<double> parsForErrors(9);
    parsForErrors[0] = fModel->GetResults()->GetValue(SFResultTypeNum::kLambda);
    parsForErrors[1] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaR);
    parsForErrors[2] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaL);
    parsForErrors[3] = fModel->GetResults()->GetValue(SFResultTypeNum::kKsi);
    parsForErrors[4] = fModel->GetResults()->GetValue(SFResultTypeNum::kLength);
    
    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
//         std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;

//         SLoop* loop = fData->GetTree(measurementsIDs[npoint]);

//         int        nloopMax = loop->getEntries();
//         SCategory* tSig     = SCategoryManager::getCategory(SCategory::CatDDSamples);

//         //----- setting energy cut
//         auto specAv = std::unique_ptr<TH1D>(fData->GetCustomHistogram(SFSelectionType::kPEAverage, cut, measurementsIDs[npoint]));
//         
//         auto peakFinAv = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specAv.get(), false));
//         peakFinAv->FindPeakRange(xmin, xmax);
        
        //----- setting histograms
//         hname = Form("hRecoPositionsCorr_S%i_pos%.1f", fSeriesNo, positions[npoint]);
//         fRecoPositionsCorrHist.push_back(new TH1D(hname, hname, 300, -100, 200));
//         fRecoPositionsCorrHist[npoint]->SetTitle(Form("Reconstructed Position (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));
// 
//         hname = Form("hRecoPositionsUncertCorr_S%i_pos%.1f", fSeriesNo, positions[npoint]);
//         fRecoPositionsUncertCorrHist.push_back(new TH1D(hname, hname, 200, 0, 50));
//         fRecoPositionsUncertCorrHist[npoint]->SetTitle(Form("Uncertainties of Reconstructed Position (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));
        
//         TString cutCh0 = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
//         TString cutCh1 = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);
//         
//         auto specCh0 = std::unique_ptr<TH1D>(fData->GetSpectrum(0, SFSelectionType::kPE, cutCh0, measurementsIDs[npoint]));
//         auto specCh1 = std::unique_ptr<TH1D>(fData->GetSpectrum(1, SFSelectionType::kPE, cutCh1, measurementsIDs[npoint]));
//         
//         auto peakFinCh0 = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specCh0.get(), false));
//         peakFinCh0->FindPeakFit();
//         
//         auto peakFinCh1 = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specCh1.get(), false));
//         peakFinCh1->FindPeakFit();
//         
//         parsForErrors[7] = peakFinCh0->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
//         parsForErrors[8] = peakFinCh1->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
        
        //----- filling histograms
//         for (int nloop = 0; nloop < nloopMax; nloop++)
//         {
//             loop->getEvent(nloop);
//             size_t tentriesMax = tSig->getEntries();
//             
//             for (int tentries = 0; tentries < tentriesMax; tentries++)
//             {
//                 int m, l, f;
//                 SDDSamples* samples = (SDDSamples*)tSig->getObject(tentries);
//                 samples->getAddress(m, l, f);
//                 
//                 if (m == 0)
//                 {
//                     double t0Ch0 = samples->getSignalL()->GetT0();
//                     double t0Ch1 = samples->getSignalR()->GetT0();
//                     double peCh0  = samples->getSignalL()->GetPE();
//                     double peCh1  = samples->getSignalR()->GetPE();
//                     double blCh0  = samples->getSignalL()->GetBLSigma();
//                     double blCh1  = samples->getSignalR()->GetBLSigma();
//                     double totCh0 = samples->getSignalL()->GetTOT();
//                     double totCh1 = samples->getSignalR()->GetTOT();
//                     double ampCh0 = samples->getSignalL()->GetAmplitude();
//                     double ampCh1 = samples->getSignalR()->GetAmplitude();
//                     bool   vetoCh0 = samples->getSignalL()->GetVeto();
//                     bool   vetoCh1 = samples->getSignalR()->GetVeto();
                    
//                     if (t0Ch0 > 0 && t0Ch1 > 0 &&
//                         totCh0 > 0 && totCh1 > 0 &&
//                         blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
//                         ampCh0 < ampMax && ampCh1 < ampMax &&
//                         sqrt(peCh0 * peCh1) > xmin &&
//                         sqrt(peCh0 * peCh1) < xmax &&
//                         vetoCh0 == 0 && vetoCh1 == 0)
//                     {
//                         parsForErrors[5] = peCh0;
//                         parsForErrors[6] = peCh1;
//                         
//                         //for (int k=0; k<9; k++)
//                         //{
//                         //    std::cout << parsForErrors[k] << "\t";
//                         //}
//                         //std::cout << std::endl;
//                         
//                         double MLR = log(sqrt(fPrRecoFun->Eval(peCh1, peCh0) /
//                                               fPlRecoFun->Eval(peCh1, peCh0)));
//                         double pos = A * MLR + B;
//                         double Pr_err = fModel->CalculateUncertainty(parsForErrors, "R");
//                         double Pl_err = fModel->CalculateUncertainty(parsForErrors, "L");
//                         double pos_err = sqrt(pow(MLR * A_err, 2) +
//                                               pow(A / (2 * fPrRecoFun->Eval(peCh1, peCh0)) * Pr_err, 2) +
//                                               pow(-A / (2 * fPlRecoFun->Eval(peCh1, peCh0)) * Pl_err, 2) + 
//                                               pow(B_err, 2));
//                          fRecoPositionsCorrHist[npoint]->Fill(pos);
//                          fRecoPositionsUncertCorrHist[npoint]->Fill(pos_err);
                          hPosRecoAll->Fill(pos - positions[npoint]);
//                     }
/*                    
                }
            }

        }*/

//         delete loop;

        double mean, mean_err, fwhm, fwhm_err;
        
        SFTools::FitGaussSingle(fRecoPositionsCorrHist[npoint], 5);
        mean     = fRecoPositionsCorrHist[npoint]->GetFunction("fGauss")->GetParameter(1);
        mean_err = fRecoPositionsCorrHist[npoint]->GetFunction("fGauss")->GetParError(1);
        fwhm     = 2 * sqrt(2 * log(2)) *fRecoPositionsCorrHist[npoint]->GetFunction("fGauss")->GetParameter(2);
        fwhm_err = 2 * sqrt(2 * log(2)) * fRecoPositionsCorrHist[npoint]->GetFunction("fGauss")->GetParError(2);

        fPosRecoCorrGraph->SetPoint(npoint, positions[npoint], mean);
        fPosRecoCorrGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), mean_err);
        
        fPosResCorrGraph->SetPoint(npoint, positions[npoint], fwhm);
        fPosResCorrGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), fwhm_err);
        
        fPosRecoDiffCorr->SetPoint(npoint, positions[npoint], mean - positions[npoint]);
        fPosRecoDiffCorr->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), fwhm);
        
        posResSum    += fwhm * (1. / pow(fwhm_err, 2));
        posResSumErr += 1. / pow(fwhm_err, 2);
    }

    double posResAv    = posResSum / posResSumErr;
    double posResAvErr = sqrt(1. / posResSumErr);

    std::cout << "Average position resolution for this series is: ";
    std::cout << posResAv << " +/- " << posResAvErr << " mm\n\n" << std::endl;

    double gconst  = hPosRecoAll->GetBinContent(hPosRecoAll->GetMaximumBin());
    double fit_min = hPosRecoAll->GetBinCenter(2);
    double mean    = hPosRecoAll->GetMean();
    double rms     = hPosRecoAll->GetRMS();
    
    TF1 *fun_gauss = new TF1("fun_gauss", "gaus", fit_min, 200);
    fun_gauss->SetParameters(gconst, mean, rms);
    hPosRecoAll->Fit(fun_gauss, "QR");
    
    TF1* funpol1 = new TF1("funpol1", "pol1", 0, 100);
    fPosRecoCorrGraph->Fit(funpol1, "Q");
    
    //----- position reconstruction end
    double res;
    double point_x, point_y;

    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        fPosRecoCorrGraph->GetPoint(npoint, point_x, point_y);
        res = point_y - positions[npoint];
        fPosResidualsCorr->SetPoint(npoint, positions[npoint], res);
    }
    
    //----- setting results
    fResultsCorr->AddResult(SFResultTypeNum::kPositionRes, posResAv, posResAvErr);
    fResultsCorr->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, fPosRecoCorrGraph);
    fResultsCorr->AddObject(SFResultTypeObj::kPosResVsPosGraph, fPosResCorrGraph);
    fResultsCorr->AddObject(SFResultTypeObj::kResidualGraph, fPosResidualsCorr);
    fResultsCorr->AddObject(SFResultTypeObj::kPositionAllHist, hPosRecoAll);
//     fResultsCorr->AddObject(SFResultTypeObj::kPositionDiff, fPosRecoDiffCorr);
    
    fResultsExp->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, fPosRecoGraph);
    fResultsExp->AddObject(SFResultTypeObj::kPosResVsPosGraph, fPosResGraph);
    fResultsExp->AddObject(SFResultTypeObj::kResidualGraph, fPosResiduals);
    fResultsExp->AddObject(SFResultTypeObj::kPositionAllHist, fPosRecoAll);
//     fResultsExp->AddObject(SFResultTypeObj::kPositionDiff, fPosRecoDiff);
    //-----
    
    return true;
}
//------------------------------------------------------------------
SFResults* SFPositionRecoModel::ReconstructPositionDistAll(SFChAddr addr, 
                                                           std::vector<SLoop*> loop,
                                                           std::vector<TH1D*> spectra_av,
                                                           std::vector<TH1D*> spectra_l,
                                                           std::vector<TH1D*> spectra_r,
                                                           std::vector<double> positions,
                                                           std::vector<TString> path,
                                                           TMatrixD covMatrix,
                                                           TF1 *fPrReco,
                                                           TF1 *fPlReco,
                                                           TGraphErrors *gAttCorrL,
                                                           TGraphErrors *gAttCorrR,
                                                           double BL_sigma_cut,
                                                           double pos_uncert,
                                                           TString collimator)
{
    int npointsMax = positions->GetN();
    
    TString hname = Form("Summed Corrected Reconstructed Position Distribution (%s)", suffix.Data());
    TH1D *hPosRecoSum = new TH1D("hPosRecoSum", hname, 300, -100, 200);
    hPosRecoSum->GetXaxis()->SetTitle("reconstructed position - source position [mm]");
    hPosRecoSum->GetYaxis()->SetTitle("counts");
    
}
//------------------------------------------------------------------
