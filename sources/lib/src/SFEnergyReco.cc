 // *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFEnergyReco.c  c          *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFEnergyReco.hh"

const double ampMax = 660;

//------------------------------------------------------------------
SFEnergyReco::SFEnergyReco(int seriesNo) : fSeriesNo(seriesNo), 
                                           fData(nullptr),
                                           fModel(nullptr),
                                           fMAttCh0Graph(nullptr),
                                           fMAttCh1Graph(nullptr),
                                           fMAttCh0CorrGraph(nullptr),
                                           fMAttCh1CorrGraph(nullptr),
                                           fPlRecoFun(nullptr),
                                           fPrRecoFun(nullptr),
                                           fResultsExp(nullptr),
                                           fResultsCorr(nullptr),
                                           fEref(511.0)
{
    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFEnergyReco constructor!";
    }

    TString desc = fData->GetDescription();

    if (!desc.Contains("Regular series"))
    {
        std::cout << "##### Error in SFEnergyReco constructor! Non-regular series!"
                  << std::endl;
        throw "##### Exception in SFEnergyReco constructor!";
    }

    try
    {
        fModel = new SFAttenuationModel(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFEnergyReco constructor!";
    }

    fModel->FitModel();
    
    SFResults* model_res = fModel->GetResults();

    fMAttCh0Graph     = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kSlVsPosGraph);
    fMAttCh1Graph     = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kSrVsPosGraph);
    fMAttCh0CorrGraph = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kPlVsPosGraph);
    fMAttCh1CorrGraph = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kPrVsPosGraph);

    fPlRecoFun = (TF2*)model_res->GetObject(SFResultTypeObj::kPlRecoFun);
    fPrRecoFun = (TF2*)model_res->GetObject(SFResultTypeObj::kPrRecoFun);
    
    fResultsExp  = new SFResults(Form("EnergyRecoResults_S%i_Exp", fSeriesNo));
    fResultsCorr = new SFResults(Form("EnergyRecoResults_S%i_Corr", fSeriesNo));
}
//------------------------------------------------------------------
SFEnergyReco::~SFEnergyReco()
{
    if (fData != nullptr) 
        delete fData;
};
//------------------------------------------------------------------
bool SFEnergyReco::CalculateAlpha(void)
{
    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    
    //----- setting up graphs start
    TGraphErrors* gEnergyAlpha = new TGraphErrors(npoints);
    gEnergyAlpha->SetName("gEnergyAlpha");
    gEnergyAlpha->SetTitle("Alpha Factor for Energy Reconstruction");
    gEnergyAlpha->GetXaxis()->SetTitle("source position [mm]");
    gEnergyAlpha->GetYaxis()->SetTitle("#alpha [keV/PE]");
    gEnergyAlpha->SetMarkerStyle(8);

    TGraphErrors* gEnergyAlphaCorr = new TGraphErrors(npoints);
    gEnergyAlphaCorr->SetName("gEnergyAlphaCorr");
    gEnergyAlphaCorr->SetTitle("Alpha Factor for Energy Reconstruction (Corrected)");
    gEnergyAlphaCorr->GetXaxis()->SetTitle("source position [mm]");
    gEnergyAlphaCorr->GetYaxis()->SetTitle("#alpha [keV/PE]");
    gEnergyAlphaCorr->SetMarkerStyle(8);
    //----- setting up graphs end
    
    //----- calculating alpha start
    double alpha, alpha_err;
    double alpha_corr, alpha_corr_err;

    for (int i = 0; i < npoints; i++)
    {
        alpha = fEref / sqrt(fMAttCh0Graph->GetPointY(i) * fMAttCh1Graph->GetPointY(i));
        alpha_err = sqrt(pow(- (fEref * fMAttCh1Graph->GetPointY(i)) /
                    (2 * pow(fMAttCh1Graph->GetPointY(i) *
                    fMAttCh0Graph->GetPointY(i), 1.5)), 2) *
                    pow(fMAttCh0Graph->GetErrorY(i), 2) +
                    pow(- (fEref * fMAttCh0Graph->GetPointY(i)) /
                    (2 * pow(fMAttCh1Graph->GetPointY(i) *
                    fMAttCh0Graph->GetPointY(i), 1.5)), 2) *
                    pow(fMAttCh1Graph->GetErrorY(i), 2));
        gEnergyAlpha->SetPoint(i, positions[i], alpha);
        gEnergyAlpha->SetPointError(i, SFTools::GetPosError(collimator, testBench), alpha_err);
        
        alpha_corr = fEref / sqrt(fMAttCh0CorrGraph->GetPointY(i) * fMAttCh1CorrGraph->GetPointY(i));
        alpha_corr_err = sqrt(pow(- (fEref * fMAttCh1CorrGraph->GetPointY(i)) /
                         (2 * pow(fMAttCh1CorrGraph->GetPointY(i) *
                         fMAttCh0CorrGraph->GetPointY(i), 1.5)), 2) *
                         pow(fMAttCh0CorrGraph->GetErrorY(i), 2) + 
                         pow(- (fEref * fMAttCh0CorrGraph->GetPointY(i)) /
                         (2 * pow(fMAttCh1CorrGraph->GetPointY(i) *
                         fMAttCh0CorrGraph->GetPointY(i), 1.5)), 2) *
                         pow(fMAttCh1CorrGraph->GetErrorY(i), 2));
        gEnergyAlphaCorr->SetPoint(i, positions[i], alpha_corr);
        gEnergyAlphaCorr->SetPointError(i, SFTools::GetPosError(collimator, testBench), alpha_corr_err);
    }

    TF1 *fun_alpha = new TF1("fun_alpha", "pol0", positions[0], positions[npoints-1]);
    gEnergyAlpha->Fit(fun_alpha, "QR+");
    alpha = fun_alpha->GetParameter(0);
    alpha_err = fun_alpha->GetParError(0);
    
    TF1 *fun_alpha_corr = new TF1("fun_alpha_corr", "pol0", positions[0], positions[npoints-1]);
    gEnergyAlphaCorr->Fit(fun_alpha_corr, "QR+");
    alpha_corr = fun_alpha_corr->GetParameter(0);
    alpha_corr_err = fun_alpha_corr->GetParError(0);
    //----- calculating alpha end
    
    //----- setting results start
    fResultsCorr->AddResult(SFResultTypeNum::kAlpha, alpha_corr, alpha_corr_err);
    fResultsCorr->AddObject(SFResultTypeObj::kAlphaGraph, gEnergyAlphaCorr);
    
    fResultsExp->AddResult(SFResultTypeNum::kAlpha, alpha, alpha_err);
    fResultsExp->AddObject(SFResultTypeObj::kAlphaGraph, gEnergyAlpha);
    //----- setting results end
    
    return true;
}
//------------------------------------------------------------------
bool SFEnergyReco::EnergyReco(void)
{
    int                 npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();

    //----- setting up graphs start
    TGraphErrors* gEnergyReco = new TGraphErrors(npoints);
    gEnergyReco->SetName("gEnergyReco");
    gEnergyReco->SetTitle("Energy Reconstruction");
    gEnergyReco->GetXaxis()->SetTitle("source position [mm]");
    gEnergyReco->GetYaxis()->SetTitle("E_{reco}/E_{real}");
    gEnergyReco->SetMarkerStyle(8);

    TGraphErrors* gEnergyRecoCorr = new TGraphErrors(npoints);
    gEnergyRecoCorr->SetName("gEnergyRecoCorr");
    gEnergyRecoCorr->SetTitle("Energy Reconstruction (Corrected)");
    gEnergyRecoCorr->GetXaxis()->SetTitle("source position [mm]");
    gEnergyRecoCorr->GetYaxis()->SetTitle("E_{reco}/E_{real}");
    gEnergyRecoCorr->SetMarkerStyle(8);
    //----- setting up graphs end

    //----- energy reconstruction start
    double       e, e_err;
    double       ecorr, ecorr_err;
    double       fiberLen = fData->GetFiberLength();

    double alpha          = fResultsExp->GetValue(SFResultTypeNum::kAlpha);
    double alpha_err      = fResultsExp->GetUncertainty(SFResultTypeNum::kAlpha);
    double alpha_corr     = fResultsCorr->GetValue(SFResultTypeNum::kAlpha);
    double alpha_corr_err = fResultsCorr->GetUncertainty(SFResultTypeNum::kAlpha);

    for (int i = 0; i < npoints; i++)
    {
        double q_av = sqrt(fMAttCh0Graph->GetPointY(i) * fMAttCh1Graph->GetPointY(i));
        e = alpha * q_av;
        e_err = sqrt(pow((q_av * alpha_err), 2) +
                pow((alpha * fMAttCh1Graph->GetPointY(i)) / 
                (2 * q_av), 2) *
                pow(fMAttCh0Graph->GetErrorY(i), 2) +
                pow((alpha * fMAttCh0Graph->GetPointY(i)) /
                (2 * q_av), 2) *
                pow(fMAttCh1Graph->GetErrorY(i), 2));
        gEnergyReco->SetPoint(i, positions[i], e / fEref);
        gEnergyReco->SetPointError(i, SFTools::GetPosError(collimator, testBench), e_err / fEref);

        q_av = sqrt(fMAttCh0CorrGraph->GetPointY(i) * fMAttCh1CorrGraph->GetPointY(i));
        ecorr = alpha_corr * q_av;
        ecorr_err = sqrt(pow(q_av * alpha_corr_err, 2) +
                    pow((alpha_corr * fMAttCh1CorrGraph->GetPointY(i)) / 
                    (2 * q_av), 2) *
                    pow(fMAttCh0CorrGraph->GetErrorY(i), 2) +
                    pow((alpha_corr * fMAttCh0CorrGraph->GetPointY(i)) /
                    (2 * q_av), 2) *
                    pow(fMAttCh1CorrGraph->GetErrorY(i), 2));
        gEnergyRecoCorr->SetPoint(i, positions[i], ecorr / fEref);
        gEnergyRecoCorr->SetPointError(i, SFTools::GetPosError(collimator, testBench), ecorr_err / fEref);
    }
    
    TF1* funPol0 = new TF1("funPol0", "pol0", 0, fiberLen);
    gEnergyRecoCorr->Fit(funPol0, "QR+");

    TF1* funEReco = new TF1("funEReco", "sqrt(fun_Sl(x)*fun_Sr(x))*[6]/[7]", 0, fiberLen); 
    funEReco->FixParameter(6, alpha);
    funEReco->FixParameter(7, fEref);
    //----- energy reconstruction end

    //----- setting results start
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyRecoFun, funPol0);
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyRecoGraph, gEnergyRecoCorr);
    
    fResultsExp->AddObject(SFResultTypeObj::kEnergyRecoFun, funEReco);
    fResultsExp->AddObject(SFResultTypeObj::kEnergyRecoGraph, gEnergyReco);
    //----- setting results end

    return true;
}
//------------------------------------------------------------------
bool SFEnergyReco::EnergyRecoByEvent(void)
{    
    
    int                 npointsMax = fData->GetNpoints();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();
    TString             sipm       = fData->GetSiPM();
    std::vector<double> positions  = fData->GetPositions();
    std::vector<int>    measurementsIDs = fData->GetMeasurementsIDs();
    
    double              s     = SFTools::GetSigmaBL(sipm);
    std::vector<double> sigma = {s};
    
    //----- setting up graphs start
    TGraphErrors* gEnergyReco = new TGraphErrors(npointsMax);
    gEnergyReco->SetName("gEnergyReco_fromSpectrum");
    gEnergyReco->SetTitle("(E_{reco} - E_{ref}) vs. Position");
    gEnergyReco->GetXaxis()->SetTitle("source position [mm]");
    gEnergyReco->GetYaxis()->SetTitle("E_{reco} - E_{ref}");
    gEnergyReco->SetMarkerStyle(8);
    
    TGraphErrors* gEnergyRecoCorr = new TGraphErrors(npointsMax);
    gEnergyRecoCorr->SetName("gEnergyRecoCorr_fromSpectrum");
    gEnergyRecoCorr->SetTitle("(E_{reco corr} - E_{ref}) vs. Position");
    gEnergyRecoCorr->GetXaxis()->SetTitle("source position [mm]");
    gEnergyRecoCorr->GetYaxis()->SetTitle("E_{reco corr} - E_{ref}");
    gEnergyRecoCorr->SetMarkerStyle(8);
    
    TGraphErrors* gEnergyRes = new TGraphErrors(npointsMax);
    gEnergyRes->SetName("fEnergyRes");
    gEnergyRes->SetTitle("Energy Resolution");
    gEnergyRes->GetXaxis()->SetTitle("source position [mm]");
    gEnergyRes->GetYaxis()->SetTitle("energy resolution [%]");
    gEnergyRes->SetMarkerStyle(8);
    
    TGraphErrors* gEnergyResCorr = new TGraphErrors(npointsMax);
    gEnergyResCorr->SetName("fEnergyResCorr");
    gEnergyResCorr->SetTitle("Energy Resolution (Corrected)");
    gEnergyResCorr->GetXaxis()->SetTitle("source position [mm]");
    gEnergyResCorr->GetYaxis()->SetTitle("energy resolution [%]");
    gEnergyResCorr->SetMarkerStyle(8);
    
    TString hname_e = Form("S%i_hEnergyRecoAllExp", fSeriesNo);
    TH1D *hEnergySpecAll = new TH1D(hname_e , hname_e, 500, 0, 1300);
    hEnergySpecAll->GetXaxis()->SetTitle("energy [keV]");
    hEnergySpecAll->GetYaxis()->SetTitle("counts");
    hEnergySpecAll->SetTitle(Form("Energy Reconstruction Spectrum (Summed) S%i", fSeriesNo));
    
    TString hname_c = Form("S%i_hEnergyRecoAllCorr", fSeriesNo);
    TH1D *hEnergySpecAllCorr = new TH1D(hname_c, hname_c, 500, 0, 1300);
    hEnergySpecAllCorr->GetXaxis()->SetTitle("energy [keV]");
    hEnergySpecAllCorr->GetYaxis()->SetTitle("counts");
    hEnergySpecAllCorr->SetTitle(Form("Energy Reconstruction Spectrum (Summed & Corrected) S%i", fSeriesNo));
    //----- setting up graphs end
    
    //----- energy reconstruction event by event start
    double BL_sigma_cut = SFTools::GetSigmaBL(sipm);

    double alpha          = fResultsExp->GetValue(SFResultTypeNum::kAlpha);
    double alpha_corr     = fResultsCorr->GetValue(SFResultTypeNum::kAlpha);
    double alpha_corr_err = fResultsCorr->GetUncertainty(SFResultTypeNum::kAlpha);
    
    double eres_exp_sum    = 0;
    double eres_exp_sumerr = 0;
    double eres_cor_sum    = 0;
    double eres_cor_sumerr = 0;
    
    std::vector<double> parsForErrors(9);
    parsForErrors[0] = fModel->GetResults()->GetValue(SFResultTypeNum::kLambda);
    parsForErrors[1] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaR);
    parsForErrors[2] = fModel->GetResults()->GetValue(SFResultTypeNum::kEtaL);
    parsForErrors[3] = fModel->GetResults()->GetValue(SFResultTypeNum::kKsi);
    parsForErrors[4] = fModel->GetResults()->GetValue(SFResultTypeNum::kLength);
    
    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;
        
        //----- getting tree
        SLoop*     loop     = fData->GetTree(measurementsIDs[npoint]);
        int        nloopMax = loop->getEntries();
        SCategory* tSig     = SCategoryManager::getCategory(SCategory::CatDDSamples);
        
        //----- setting histograms
        hname_e = Form("S%i_hEnergyRecoExp_pos%.1f", fSeriesNo, positions[npoint]);
        hname_c = Form("S%i_hEnergyRecoCorr_pos%.1f", fSeriesNo, positions[npoint]);
        
        fEnergySpectra.push_back(new TH1D(hname_e, hname_e, 750, 0, 1300));
        fEnergySpectraCorr.push_back(new TH1D(hname_c, hname_c, 750, 0, 1300));
        
        fEnergySpectra[npoint]->SetTitle(Form("Energy Reconstruction Spectrum S%i %.1f mm", fSeriesNo, positions[npoint]));
        fEnergySpectraCorr[npoint]->SetTitle(Form("Energy Reconstruction Spectrum (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));
        
        hname_c = Form("S%i_hEnergyRecoCorrUncert_pos%.1f", fSeriesNo, positions[npoint]);
        fEnergyUncertDistCorr.push_back(new TH1D(hname_c, hname_c, 750, 0, 500));
        fEnergyUncertDistCorr[npoint]->SetTitle(Form("Reconstructed Energy Uncertainty Distribution (Corrected) S%i %.1f mm", fSeriesNo, positions[npoint]));

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
                        //blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
                        ampCh0 < ampMax && ampCh1 < ampMax &&
                        peCh0 > 0 && peCh1 > 0)// &&
                        //vetoCh0 == 0 && vetoCh1 == 0)
                    {
                        double e_reco = alpha * sqrt(peCh0 * peCh1);
                        
                        parsForErrors[5] = peCh0;
                        parsForErrors[6] = peCh1;
                        
//                         for (int k=0; k<9; k++)
//                         {
//                             std::cout << parsForErrors[k] << "\t";
//                         }
//                         std::cout << std::endl;
                        
                        double Pl_err = fModel->CalculateUncertainty(parsForErrors, "L");
                        double Pr_err = fModel->CalculateUncertainty(parsForErrors, "R");
                        
                        double q_av = sqrt(fPlRecoFun->Eval(peCh1, peCh0) *
                                           fPrRecoFun->Eval(peCh1, peCh0));
                        double e_corr_reco = alpha_corr * q_av;
                        double e_reco_corr_err = sqrt(pow(q_av * alpha_corr_err, 2) + 
                                                      pow((alpha_corr * fPrRecoFun->Eval(peCh1, peCh0)) / (2 * q_av) * Pl_err, 2) + 
                                                      pow((alpha_corr * fPlRecoFun->Eval(peCh1, peCh0)) / (2 * q_av) * Pr_err, 2));
                        
                        fEnergySpectra[npoint]->Fill(e_reco);
                        fEnergySpectraCorr[npoint]->Fill(e_corr_reco);
                        fEnergyUncertDistCorr[npoint]->Fill(e_reco_corr_err);
                        
                        hEnergySpecAll->Fill(e_reco);
                        hEnergySpecAllCorr->Fill(e_corr_reco);
                    }
                    
                }
            }
        }
        
        double mean, mean_err;
        double sigma, sigma_err;
        
        fEnergySpectra[npoint]->Rebin(2);
        
        auto pf_exp = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(fEnergySpectra[npoint], measurementsIDs[npoint], false, true));
        pf_exp->FindPeakFit();
        
        mean      = pf_exp->GetResults()->GetValue(SFResultTypeNum::kPeakPosition);
        mean_err  = pf_exp->GetResults()->GetUncertainty(SFResultTypeNum::kPeakPosition);
        sigma     = pf_exp->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
        sigma_err = pf_exp->GetResults()->GetUncertainty(SFResultTypeNum::kPeakSigma);
        double eres_exp = sigma / mean;
        double eres_exp_err = eres_exp * sqrt(pow((mean_err / mean), 2) 
                                       + pow((sigma_err / sigma), 2));
        
        gEnergyReco->SetPoint(npoint, positions[npoint], mean - fEref);
        gEnergyReco->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), sigma);
        gEnergyRes->SetPoint(npoint, positions[npoint], eres_exp * 100);
        gEnergyRes->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), eres_exp_err * 100);
        
        eres_exp_sum += eres_exp * (1. / pow(eres_exp_err, 2));
        eres_exp_sumerr += (1. / pow(eres_exp_err, 2));
        
        fEnergySpectraCorr[npoint]->Rebin(2);
        
        auto pf_corr = std::unique_ptr<SFPeakFinder>(new SFPeakFinder(fEnergySpectraCorr[npoint],  measurementsIDs[npoint], false, true));
        pf_corr->FindPeakFit();
        
        mean      = pf_corr->GetResults()->GetValue(SFResultTypeNum::kPeakPosition);
        mean_err  = pf_corr->GetResults()->GetUncertainty(SFResultTypeNum::kPeakPosition);
        sigma     = pf_corr->GetResults()->GetValue(SFResultTypeNum::kPeakSigma);
        sigma_err = pf_corr->GetResults()->GetUncertainty(SFResultTypeNum::kPeakSigma);
        double eres_cor = sigma / mean;
        double eres_cor_err = eres_cor * sqrt(pow((mean_err / mean), 2) 
                                       + pow((sigma_err / sigma), 2));
        
        gEnergyRecoCorr->SetPoint(npoint, positions[npoint], mean - fEref);
        gEnergyRecoCorr->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), sigma);
        gEnergyResCorr->SetPoint(npoint, positions[npoint], eres_cor * 100);
        gEnergyResCorr->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), eres_cor_err * 100);
        
        eres_cor_sum += eres_cor * (1. / pow(eres_cor_err, 2));
        eres_cor_sumerr += (1. / pow(eres_cor_err, 2));
    }
    //----- energy reconstruction end
    
    double eres_exp_av    = (eres_exp_sum / eres_exp_sumerr) * 100; 
    double eres_exp_averr = (sqrt(1. / eres_exp_sumerr)) * 100;
    
    double eres_cor_av    = (eres_cor_sum / eres_cor_sumerr) * 100;
    double eres_cor_averr = (sqrt(1. / eres_cor_sumerr)) * 100;
    
    //----- setting results start
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyRecoSpecGraph, gEnergyRecoCorr);
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyResGraph, gEnergyResCorr);
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyAllHist, hEnergySpecAllCorr);
    fResultsCorr->AddResult(SFResultTypeNum::kEnergyRes, eres_cor_av, eres_cor_averr);
    
    fResultsExp->AddObject(SFResultTypeObj::kEnergyRecoSpecGraph, gEnergyReco);
    fResultsExp->AddObject(SFResultTypeObj::kEnergyResGraph, gEnergyRes);
    fResultsExp->AddObject(SFResultTypeObj::kEnergyAllHist, hEnergySpecAll);
    fResultsExp->AddResult(SFResultTypeNum::kEnergyRes, eres_exp_av, eres_exp_averr);
    //----- setting results end

    return true;
}
//------------------------------------------------------------------
void SFEnergyReco::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFEnergyReco class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
std::vector<SFResults*> SFEnergyReco::GetResults(void)
{
    if (fResultsExp == nullptr ||
        fResultsCorr == nullptr)
    {
        std::cerr << "##### Error in SFReconstruction::GetRecoResults()!" << std::endl;
        std::cerr << "Empty SFResults object pointer!" << std::endl;
        std::abort();
    }

    std::vector<SFResults*> results(2);
    results[0] = fResultsExp;
    results[1] = fResultsCorr;

    return results;
}
//------------------------------------------------------------------
std::vector<TH1D*> SFEnergyReco::GetEnergySpectra(TString type)
{
    if ((type == "experimental" && fEnergySpectra.empty()) || 
        (type == "corrected" && fEnergySpectraCorr.empty()))
    {
        std::cerr << "##### Error in SFEnergyReco::GetEnergyRecoHistograms(). Empty vector!" << std::endl;
        std::abort();
    }

    std::vector<TH1D*> hist;
    
    if (type == "experimental")
        hist = fEnergySpectra;
    else if (type == "corrected")
        hist = fEnergySpectraCorr;
    else
    {
        std::cerr << "##### Error in SFEnergyReco::GetEnergyRecoHistograms()!" << std::endl;
        std::cerr << "Incorrect type. Possible options are: experimental and corrected" << std::endl;
        std::abort();
    }
    
    return hist;
}
//------------------------------------------------------------------
