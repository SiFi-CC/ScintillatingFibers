// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFPositionReco.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFPositionReco.hh"

const double ampMax = 660.;

//------------------------------------------------------------------
SFPositionReco::SFPositionReco(int seriesNo) : fSeriesNo(seriesNo), 
                                               fData(nullptr),
                                               fModel(nullptr),
                                               fMAttCh0CorrGraph(nullptr),
                                               fMAttCh1CorrGraph(nullptr),
                                               fMLRGraph(nullptr),
                                               fMLRCorrGraph(nullptr),
                                               fPosRecoGraph(nullptr),
                                               fPosRecoCorrGraph(nullptr),
                                               fPosResiduals(nullptr),
                                               fPosResidualsCorr(nullptr),
                                               fPosRecoDiff(nullptr),
                                               fPosRecoDiffCorr(nullptr),
                                               fPosResGraph(nullptr),
                                               fPosResCorrGraph(nullptr),
                                               fPlRecoFun(nullptr),
                                               fPrRecoFun(nullptr),
                                               fPosRecoAll(nullptr),
                                               fResultsExp(nullptr),
                                               fResultsCorr(nullptr)
{
    try
    {
       fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFPositionReco constructor!";
    }

    TString desc = fData->GetDescription();

    if (!desc.Contains("Regular series"))
    {
        std::cout << "##### Error in SFPositionReco constructor! Non-regular series!"
                  << std::endl;
        throw "##### Exception in SFPositionReco constructor!";
    }

    fResultsExp  = new SFResults(Form("PositionRecoResults_S%i_Exp", fSeriesNo));
    fResultsCorr = new SFResults(Form("PositionRecoResults_S%i_Corr", fSeriesNo));
    
    //----- accessing attenuation analysis results
    SFAttenuation* att;
     
    try
    {
        att = new SFAttenuation(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFPositionReco constructor!";
    }

    att->AttCombinedCh();

    std::vector<SFResults*> att_results = att->GetResults();
    
    fMLRGraph = (TGraphErrors*)att_results[2]->GetObject(SFResultTypeObj::kAttGraph);
    
    delete att;
    delete att_results[0];
    delete att_results[1];
    delete att_results[3];
    delete att_results[4];
    //-----
    
    //----- accessing attenuation model results
    //SFAttenuationModel* fModel;
    
    try
    {
        fModel = new SFAttenuationModel(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFPositionReco constructor!";
    }
    
    fModel->FitModel();
    
    SFResults *model_results = fModel->GetResults();
    
    fPlRecoFun = (TF2*)model_results->GetObject(SFResultTypeObj::kPlRecoFun);
    fPrRecoFun = (TF2*)model_results->GetObject(SFResultTypeObj::kPrRecoFun);
    fMAttCh0CorrGraph = (TGraphErrors*)model_results->GetObject(SFResultTypeObj::kPlVsPosGraph);
    fMAttCh1CorrGraph = (TGraphErrors*)model_results->GetObject(SFResultTypeObj::kPrVsPosGraph);
    
    //delete model_results;
    //-----
    
    //----- accessing position resolution analysis results 
    SFPositionRes *posres;
    
    try
    {
        posres = new SFPositionRes(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFPositionReco constructor!";
    }
    
    posres->AnalyzePositionRes();
    
    std::vector<SFResults*> posres_results = posres->GetResults();
    
    // using pol1 results! change for [1] for pol3
    fPosResGraph  = (TGraphErrors*)posres_results[0]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    fPosRecoGraph = (TGraphErrors*)posres_results[0]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    fPosResiduals = (TGraphErrors*)posres_results[0]->GetObject(SFResultTypeObj::kResidualGraph);
//     fPosRecoDiff = (TGraphErrors*)posres_results[0]->GetObject(SFResultTypeObj::kPositionDiff);
    fPosRecoAll   = (TH1D*)posres_results[0]->GetObject(SFResultTypeObj::kPositionAllHist);
    
    fRecoPositionsHist = posres->GetPositionRecoDist("pol3");
    
    fResultsExp->AddResult(SFResultTypeNum::kPositionRes, posres_results[0]->GetValue(SFResultTypeNum::kPositionRes),
                           posres_results[0]->GetUncertainty(SFResultTypeNum::kPositionRes));
    //delete posres;
    //delete posres_results;
    //-----
}
//------------------------------------------------------------------
SFPositionReco::~SFPositionReco()
{
    if (fData != nullptr) delete fData;
};
//------------------------------------------------------------------
bool SFPositionReco::CalculateMLR(void)
{    
    std::cout << "\n----- Calculating Corrected MLR curve for series " << fSeriesNo << std::endl;
    
    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();

    fMLRGraph->SetName("fMLRGraph");
    fMLRGraph->SetTitle("M_{LR} vs. Position");
    fMLRGraph->GetXaxis()->SetTitle("source position [mm]");
    fMLRGraph->GetYaxis()->SetTitle("M_{LR}");
    fMLRGraph->SetMarkerStyle(8);
    
    fMLRCorrGraph = new TGraphErrors(npoints);
    fMLRCorrGraph->SetName("fMAttCorrGraph");
    fMLRCorrGraph->SetTitle("M_{LR} vs. Position (Corrected)");
    fMLRCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fMLRCorrGraph->GetYaxis()->SetTitle("M_{LR} corrected");
    fMLRCorrGraph->SetMarkerStyle(8);

    TGraphErrors *gRevMLRCorr = new TGraphErrors(npoints);
    gRevMLRCorr->SetName("gRevMLRCorr");
    gRevMLRCorr->SetTitle("Position vs. M_{LR} (Corrected)");
    gRevMLRCorr->GetXaxis()->SetTitle("M_{LR}");
    gRevMLRCorr->GetYaxis()->SetTitle("source position [mm]");
    gRevMLRCorr->SetMarkerStyle(8);
    
    if (fMAttCh0CorrGraph == nullptr || fMAttCh1CorrGraph == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco::CalculateMLR()" << std::endl;
        std::cerr << "Corrected attenuation graphs not existing!" << std::endl;
        std::abort();
    }
    
    for (int i = 0; i < npoints; i++)
    {
        double vCh0    = fMAttCh0CorrGraph->GetPointY(i);
        double vCh0err = fMAttCh0CorrGraph->GetErrorY(i);
        double vCh1    = fMAttCh1CorrGraph->GetPointY(i);
        double vCh1err = fMAttCh1CorrGraph->GetErrorY(i);
        double MLR     = log(sqrt(vCh1 / vCh0));
        double MLRerr  = sqrt(pow(vCh1err / (2 * vCh1), 2) +
                              pow(- vCh0err / (2 * vCh0), 2));
        fMLRCorrGraph->SetPoint(i, positions[i], MLR);
        fMLRCorrGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), MLRerr);
        
        gRevMLRCorr->SetPoint(i, MLR, positions[i]);
        gRevMLRCorr->SetPointError(i, MLRerr, SFTools::GetPosError(collimator, testBench));
    }
    
    TF1 *fpol1 = new TF1("fpol1", "pol1", positions[0], positions[npoints-1]); 
    TFitResultPtr ptr = fMLRCorrGraph->Fit(fpol1, "SQR+");
    
    TF1 *fpol1_rev = new TF1("fpol1", "pol1", -1, 1);
    TFitResultPtr ptr_rev = gRevMLRCorr->Fit(fpol1_rev, "SQR+");
    
    fResultsCorr->AddObject(SFResultTypeObj::kAttGraph, fMLRCorrGraph);
    fResultsCorr->AddResult(SFResultTypeNum::kMLRSlope, fpol1->GetParameter(1), fpol1->GetParError(1));
    fResultsCorr->AddResult(SFResultTypeNum::kMLROffset, fpol1->GetParameter(0), fpol1->GetParError(0));
    fResultsCorr->AddObject(SFResultTypeObj::kPosVsMLRGraph, gRevMLRCorr);
    //fResultsCorr->AddResult(SFResultTypeNum::kACoeff, fpol1_rev->GetParameter(0), fpol1_rev->GetParError(0));
    fResultsCorr->AddResult(SFResultTypeNum::kACoeff, fpol1_rev->GetParameter(1), fpol1_rev->GetParError(1));
    
    fResultsExp->AddObject(SFResultTypeObj::kAttGraph, fMLRGraph);
    fResultsExp->AddResult(SFResultTypeNum::kMLRSlope, fMLRGraph->GetFunction("fpol1")->GetParameter(1),
                           fMLRGraph->GetFunction("fpol1")->GetParError(1));
    fResultsExp->AddResult(SFResultTypeNum::kMLROffset, fMLRGraph->GetFunction("fpol1")->GetParameter(0),
                           fMLRGraph->GetFunction("fpol1")->GetParError(0));

    return true;
}
//------------------------------------------------------------------
bool SFPositionReco::CalculatePosRecoCoefficients(void)
{
    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    
    double B = fData->GetFiberLength() / 2.;
    
    int n = 0;
    
    for (int i = 0; i < npoints; i++)
    {
        if (fabs(B - positions[i]) < 1E-2)
            continue;
        n++;
    }
    
    std::cout << "\n----- Calculating Position Reconstruction Coefficients for series " << fSeriesNo << std::endl;
    std::cout << "Graph has " << n << " points" << std::endl; 
    
    fAGraph = new TGraphErrors(n);
    fAGraph->SetName("fAGraph");
    fAGraph->SetTitle("A Coefficient vs. Position");
    fAGraph->GetXaxis()->SetTitle("source position [mm]");
    fAGraph->GetYaxis()->SetTitle("A [mm]");
    fAGraph->SetMarkerStyle(8);
    
    if (fMLRCorrGraph == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco::CalculatePosRecoCoefficients()" << std::endl;
        std::cerr << "MLR graph does not exist!" << std::endl;
        std::abort();
    }
    
    int j = -1;
    
    for (int i = 0; i < npoints; i++)
    {
        double val   = fMLRCorrGraph->GetPointY(i);
        double val_e = fMLRCorrGraph->GetErrorY(i);
        if (fabs(B - positions[i]) < 1E-2)
            continue;
        j++;
        double A    = (positions[i] - B) / (val);
        double Aerr = sqrt(pow((SFTools::GetPosError(collimator, testBench) / val), 2) +
                           pow((B - positions[i]) * val_e / (pow(val, 2)), 2));
        fAGraph->SetPoint(j, positions[i], A);
        fAGraph->SetPointError(j, SFTools::GetPosError(collimator, testBench), Aerr);
    }
    
    TF1* fpol0 = new TF1("fpol0", "pol0", positions[0], positions[npoints-1]);
    fAGraph->Fit(fpol0, "SQR+");
    
    fResultsCorr->AddObject(SFResultTypeObj::kAGraph, fAGraph);
    //fResultsCorr->AddResult(SFResultTypeNum::kACoeff, fpol0->GetParameter(0), fpol0->GetParError(0));
    fResultsCorr->AddResult(SFResultTypeNum::kBCoeff, B, 0);
    //fResultsCorr->AddResult(SFResultTypeNum::kACoeff, fAGraph->GetPointY(0), fAGraph->GetErrorY(0));
    
    std::cout << "A = " << fpol0->GetParameter(0) << " +/- " << fpol0->GetParError(0) << std::endl;
    std::cout << "B = " << B << std::endl;
    
    return true;
}
//------------------------------------------------------------------
bool SFPositionReco::PositionReco(void)
{
    std::cout << "\n\n----- Position Resolution (Corrected) Analysis" << std::endl;
    std::cout << "----- Series: " << fSeriesNo << std::endl;
    
    int                 npointsMax = fData->GetNpoints();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             sipm       = fData->GetSiPM();
    std::vector<double> positions  = fData->GetPositions();
    std::vector<int>    measurementsIDs = fData->GetMeasurementsIDs();
    
    //----- setting graphs
    fPosRecoCorrGraph = new TGraphErrors(npointsMax);
    fPosRecoCorrGraph->SetName("fPosRecoCorrGraph");
    fPosRecoCorrGraph->SetTitle(Form("Reconstructed Source Position (Corrected) S%i", fSeriesNo));
    fPosRecoCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fPosRecoCorrGraph->GetYaxis()->SetTitle("reconstructed position [mm]");
    fPosRecoCorrGraph->SetMarkerStyle(8);
    
    fPosResCorrGraph = new TGraphErrors(npointsMax);
    fPosResCorrGraph->SetName("fPosResCorrGraph");
    fPosResCorrGraph->SetTitle(Form("Position Resolution (Corrected) S%i", fSeriesNo));
    fPosResCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fPosResCorrGraph->GetYaxis()->SetTitle("position resolution [mm]");
    fPosResCorrGraph->SetMarkerStyle(8);
    
    fPosResidualsCorr = new TGraphErrors(npointsMax);
    fPosResidualsCorr->SetName("fPosResidualsCorr");
    fPosResidualsCorr->SetTitle(Form("Reconstructed Position (Corrected) Residuals S%i", fSeriesNo));
    fPosResidualsCorr->GetXaxis()->SetTitle("source position [mm]");
    fPosResidualsCorr->GetYaxis()->SetTitle("residual [mm]");
    fPosResidualsCorr->SetMarkerStyle(8);
    
    fPosRecoDiffCorr = new TGraphErrors(npointsMax);
    fPosRecoDiffCorr->SetName("fPosRecoDiffCorr");
    fPosRecoDiffCorr->SetTitle(Form("Reconstructed Position Difference (Corrected) S%i", fSeriesNo));
    fPosRecoDiffCorr->GetXaxis()->SetTitle("source position [mm]");
    fPosRecoDiffCorr->GetYaxis()->SetTitle("P_{reco} - P_{real} [mm]");
    fPosRecoDiffCorr->SetMarkerStyle(8);
    //-----
    
    //----- position reconstruction event by event start
    double BL_sigma_cut = SFTools::GetSigmaBL(sipm);  
    double s = SFTools::GetSigmaBL(fData->GetSiPM());
    std::vector<double> sigmas = {s, s};
    std::vector<double> sigma = {s};
    TString cut = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);
   
    //-----
    double xmin, xmax;
    double A     = fResultsCorr->GetValue(SFResultTypeNum::kACoeff);
    double A_err = fResultsCorr->GetUncertainty(SFResultTypeNum::kACoeff);
    double B     = fResultsCorr->GetValue(SFResultTypeNum::kBCoeff);
    double B_err = fResultsCorr->GetUncertainty(SFResultTypeNum::kBCoeff);
    
    double posResSum    = 0.;
    double posResSumErr = 0.;
    
    TString hname = Form("Reconstructed Position Distribution (Summed) S%i", fSeriesNo);
    TH1D *hPosRecoAll = new TH1D("hPosRecoAll", hname, 300, -100, 200);
    hPosRecoAll->GetXaxis()->SetTitle("reconstructed position - source position [mm]");
    hPosRecoAll->GetYaxis()->SetTitle("counts");
    
    std::vector<double> parsForErrors(9);
    parsForErrors[0] = fModel->GetResults()->GetValue(SFResultTypeNum::kLambda);
    parsForErrors[1] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaR);
    parsForErrors[2] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaL);
    parsForErrors[3] = fModel->GetResults()->GetValue(SFResultTypeNum::kKsi);
    parsForErrors[4] = fModel->GetResults()->GetValue(SFResultTypeNum::kLength);
    
    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;

        SLoop* loop = fData->GetTree(measurementsIDs[npoint]);

        int        nloopMax = loop->getEntries();
        SCategory* tSig     = SCategoryManager::getCategory(SCategory::CatDDSamples);

        //----- setting energy cut
        auto specAv = std::unique_ptr<TH1D>(fData->GetCustomHistogram(SFSelectionType::kPEAverage, cut, measurementsIDs[npoint]));
        
        auto peakFinAv = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specAv.get(), false));
        peakFinAv->FindPeakRange(xmin, xmax);
        
        //----- setting histograms
        hname = Form("hRecoPositionsCorr_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fRecoPositionsCorrHist.push_back(new TH1D(hname, hname, 300, -100, 200));
        fRecoPositionsCorrHist[npoint]->SetTitle(Form("Reconstructed Position (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));

        hname = Form("hRecoPositionsUncertCorr_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        fRecoPositionsUncertCorrHist.push_back(new TH1D(hname, hname, 200, 0, 50));
        fRecoPositionsUncertCorrHist[npoint]->SetTitle(Form("Uncertainties of Reconstructed Position (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));
        
        TString cutCh0 = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
        TString cutCh1 = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);
        
        auto specCh0 = std::unique_ptr<TH1D>(fData->GetSpectrum(0, SFSelectionType::kPE, cutCh0, measurementsIDs[npoint]));
        auto specCh1 = std::unique_ptr<TH1D>(fData->GetSpectrum(1, SFSelectionType::kPE, cutCh1, measurementsIDs[npoint]));
        
        auto peakFinCh0 = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specCh0.get(), false));
        peakFinCh0->FindPeakFit();
        
        auto peakFinCh1 = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(specCh1.get(), false));
        peakFinCh1->FindPeakFit();
        
        parsForErrors[7] = peakFinCh0->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
        parsForErrors[8] = peakFinCh1->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
        
        //----- filling histograms
        for (int nloop = 0; nloop < nloopMax; nloop++)
        {
            loop->getEvent(nloop);
            size_t tentriesMax = tSig->getEntries();
            
            for (int tentries = 0; tentries < tentriesMax; tentries++)
            {
                int m, l, f;
                SDDSamples* samples = (SDDSamples*)tSig->getObject(tentries);
                samples->getAddress(m, l, f);
                
                if (m == 0)
                {
                    double t0Ch0 = samples->getSignalL()->GetT0();
                    double t0Ch1 = samples->getSignalR()->GetT0();
                    double peCh0  = samples->getSignalL()->GetPE();
                    double peCh1  = samples->getSignalR()->GetPE();
                    //double blCh0  = samples->getSignalL()->GetBLSigma();
                    //double blCh1  = samples->getSignalR()->GetBLSigma();
                    double totCh0 = samples->getSignalL()->GetTOT();
                    double totCh1 = samples->getSignalR()->GetTOT();
                    double ampCh0 = samples->getSignalL()->GetAmplitude();
                    double ampCh1 = samples->getSignalR()->GetAmplitude();
                    //bool   vetoCh0 = samples->getSignalL()->GetVeto();
                    //bool   vetoCh1 = samples->getSignalR()->GetVeto();
                    
                    if (t0Ch0 > 0 && t0Ch1 > 0 &&
                        totCh0 > 0 && totCh1 > 0 &&
                        //blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
                        ampCh0 < ampMax && ampCh1 < ampMax &&
                        sqrt(peCh0 * peCh1) > xmin &&
                        sqrt(peCh0 * peCh1) < xmax)// &&
                        //vetoCh0 == 0 && vetoCh1 == 0)
                    {
                        parsForErrors[5] = peCh0;
                        parsForErrors[6] = peCh1;
                        
                        //for (int k=0; k<9; k++)
                        //{
                        //    std::cout << parsForErrors[k] << "\t";
                        //}
                        //std::cout << std::endl;
                        
                        double MLR = log(sqrt(fPrRecoFun->Eval(peCh1, peCh0) /
                                              fPlRecoFun->Eval(peCh1, peCh0)));
                        double pos = A * MLR + B;
                        double Pr_err = fModel->CalculateUncertainty(parsForErrors, "R");
                        double Pl_err = fModel->CalculateUncertainty(parsForErrors, "L");
                        double pos_err = sqrt(pow(MLR * A_err, 2) +
                                              pow(A / (2 * fPrRecoFun->Eval(peCh1, peCh0)) * Pr_err, 2) +
                                              pow(-A / (2 * fPlRecoFun->Eval(peCh1, peCh0)) * Pl_err, 2) + 
                                              pow(B_err, 2));
                        fRecoPositionsCorrHist[npoint]->Fill(pos);
                        fRecoPositionsUncertCorrHist[npoint]->Fill(pos_err);
                        hPosRecoAll->Fill(pos - positions[npoint]);
                    }
                    
                }
            }

        }

        delete loop;

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
std::vector<SFResults*> SFPositionReco::GetResults(void)
{
    if (fResultsExp == nullptr ||
        fResultsCorr == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco::GetResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }

    std::vector<SFResults*> results(2);
    results[0] = fResultsExp;
    results[1] = fResultsCorr;

    return results;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFPositionReco::GetPositionDistributions(TString type)
{
    std::vector<TH1D*> histograms;
    
    if (type == "experimental" && !fRecoPositionsHist.empty())
    {
        histograms = fRecoPositionsHist;
    }
    else if (type == "corrected" && !fRecoPositionsCorrHist.empty())
    {
        histograms = fRecoPositionsCorrHist;
    }
    else
    {
        std::cerr << "##### Error in SFPositionReco::GetPositionDistributions()!" << std::endl;
        std::cerr << "Incorrect type. Possible options are: experimental and corrected." << std::endl;
        std::abort();
    }
    
    return histograms;
}
//------------------------------------------------------------------
void SFPositionReco::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFPositionReco class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
