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
                                           fMAttCh0Graph(nullptr),
                                           fMAttCh1Graph(nullptr),
                                           fMAttCh0CorrGraph(nullptr),
                                           fMAttCh1CorrGraph(nullptr),
                                           fSlFun(nullptr),
                                           fSrFun(nullptr),
                                           fPlRecoFun(nullptr),
                                           fPrRecoFun(nullptr),
                                           fEnergyRecoGraph(nullptr),
                                           fEnergyRecoCorrGraph(nullptr),
                                           fEnergyAlphaGraph(nullptr),
                                           fEnergyAlphaCorrGraph(nullptr),
                                           fEnergyDiffGraph(nullptr),
                                           fEnergyDiffCorrGraph(nullptr),
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

    SFAttenuationModel* model = nullptr;
    
    try
    {
        model = new SFAttenuationModel(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFEnergyReco constructor!";
    }

    model->FitModel();
    
    SFResults* model_res = model->GetResults();
    
    fMAttCh0Graph     = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kSlVsPosGraph);
    fMAttCh1Graph     = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kSrVsPosGraph);
    fMAttCh0CorrGraph = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kPlVsPosGraph);
    fMAttCh1CorrGraph = (TGraphErrors*)model_res->GetObject(SFResultTypeObj::kPrVsPosGraph);
    
    fSlFun = (TF1*)model_res->GetObject(SFResultTypeObj::kSlFun);
    fSrFun = (TF1*)model_res->GetObject(SFResultTypeObj::kSrFun);

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
    fEnergyAlphaGraph = new TGraphErrors(npoints);
    fEnergyAlphaGraph->SetName("fEnergyAlphaGraph");
    fEnergyAlphaGraph->SetTitle("Alpha Factor for Energy Reconstruction");
    fEnergyAlphaGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyAlphaGraph->GetYaxis()->SetTitle("#alpha [keV/PE]");
    fEnergyAlphaGraph->SetMarkerStyle(4);

    fEnergyAlphaCorrGraph = new TGraphErrors(npoints);
    fEnergyAlphaCorrGraph->SetName("fEnergyAlphaCorrGraph");
    fEnergyAlphaCorrGraph->SetTitle("Alpha Factor for Energy Reconstruction (Corrected)");
    fEnergyAlphaCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyAlphaCorrGraph->GetYaxis()->SetTitle("#alpha [keV/PE]");
    fEnergyAlphaCorrGraph->SetMarkerStyle(8);
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
                    fMAttCh0Graph->GetErrorY(i) +
                    pow(- (fEref * fMAttCh0Graph->GetPointY(i)) /
                    (2 * pow(fMAttCh1Graph->GetPointY(i) *
                    fMAttCh0Graph->GetPointY(i), 1.5)), 2) *
                    fMAttCh1Graph->GetErrorY(i));
        fEnergyAlphaGraph->SetPoint(i, positions[i], alpha);
        fEnergyAlphaGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), alpha_err);
        
        alpha_corr = fEref / sqrt(fMAttCh0CorrGraph->GetPointY(i) * fMAttCh1CorrGraph->GetPointY(i));
        alpha_corr_err = sqrt(pow(- (fEref * fMAttCh1CorrGraph->GetPointY(i)) /
                         (2 * pow(fMAttCh1CorrGraph->GetPointY(i) *
                         fMAttCh0CorrGraph->GetPointY(i), 1.5)), 2) *
                         fMAttCh0CorrGraph->GetErrorY(i) +
                         pow(- (fEref * fMAttCh0CorrGraph->GetPointY(i)) /
                         (2 * pow(fMAttCh1CorrGraph->GetPointY(i) *
                         fMAttCh0CorrGraph->GetPointY(i), 1.5)), 2) *
                         fMAttCh1CorrGraph->GetErrorY(i));
        fEnergyAlphaCorrGraph->SetPoint(i, positions[i], alpha_corr);
        fEnergyAlphaCorrGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), alpha_corr_err);
    }

    TF1 *fun_alpha = new TF1("fun_alpha", "pol0", positions[0], positions[npoints-1]);
    fEnergyAlphaGraph->Fit(fun_alpha, "QR+");
    alpha = fun_alpha->GetParameter(0);
    alpha_err = fun_alpha->GetParError(0);
    
    TF1 *fun_alpha_corr = new TF1("fun_alpha_corr", "pol0", positions[0], positions[npoints-1]);
    fEnergyAlphaCorrGraph->Fit(fun_alpha_corr, "QR+");
    alpha_corr = fun_alpha_corr->GetParameter(0);
    alpha_corr_err = fun_alpha_corr->GetParError(0);
    //----- calculating alpha end
    
    //----- setting results start
    fResultsCorr->AddResult(SFResultTypeNum::kAlpha, alpha_corr, alpha_corr_err);
    fResultsCorr->AddObject(SFResultTypeObj::kAlphaGraph, fEnergyAlphaCorrGraph);
    
    fResultsExp->AddResult(SFResultTypeNum::kAlpha, alpha, alpha_err);
    fResultsExp->AddObject(SFResultTypeObj::kAlphaGraph, fEnergyAlphaGraph);
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
    fEnergyRecoGraph = new TGraphErrors(npoints);
    fEnergyRecoGraph->SetName("fEnergyRecoGraph");
    fEnergyRecoGraph->SetTitle("Energy Reconstruction");
    fEnergyRecoGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyRecoGraph->GetYaxis()->SetTitle("E_{reco}/E_{real}");
    fEnergyRecoGraph->SetMarkerStyle(4);

    fEnergyRecoCorrGraph = new TGraphErrors(npoints);
    fEnergyRecoCorrGraph->SetName("fEnergyRecoCorrGraph");
    fEnergyRecoCorrGraph->SetTitle("Energy Reconstruction (Corrected)");
    fEnergyRecoCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyRecoCorrGraph->GetYaxis()->SetTitle("E_{reco}/E_{real}");
    fEnergyRecoCorrGraph->SetMarkerStyle(8);
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
        e = alpha * sqrt(fMAttCh0Graph->GetPointY(i) * fMAttCh1Graph->GetPointY(i));
        e_err = sqrt(pow(sqrt(fMAttCh0Graph->GetPointY(i) *
                fMAttCh1Graph->GetPointY(i)), 2) * alpha_err +
                pow((alpha_corr * fMAttCh1Graph->GetPointY(i)) / 
                (2 * sqrt(fMAttCh0Graph->GetPointY(i) *
                fMAttCh1Graph->GetPointY(i))), 2) *
                fMAttCh0Graph->GetErrorY(i) +
                pow((alpha_corr * fMAttCh0Graph->GetPointY(i)) /
                (2 * sqrt(fMAttCh0Graph->GetPointY(i) *
                fMAttCh1Graph->GetPointY(i))), 2) *
                fMAttCh1Graph->GetErrorY(i));
        fEnergyRecoGraph->SetPoint(i, positions[i], e / fEref);
        fEnergyRecoGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), e_err / fEref);

        ecorr = alpha_corr * sqrt(fMAttCh0CorrGraph->GetPointY(i) * fMAttCh1CorrGraph->GetPointY(i));
        ecorr_err = sqrt(pow(sqrt(fMAttCh0CorrGraph->GetPointY(i) *
                    fMAttCh1CorrGraph->GetPointY(i)), 2) * alpha_corr_err +
                    pow((alpha_corr * fMAttCh1CorrGraph->GetPointY(i)) / 
                    (2 * sqrt(fMAttCh0CorrGraph->GetPointY(i) *
                    fMAttCh1CorrGraph->GetPointY(i))), 2) *
                    fMAttCh0CorrGraph->GetErrorY(i) +
                    pow((alpha_corr * fMAttCh0CorrGraph->GetPointY(i)) /
                    (2 * sqrt(fMAttCh0CorrGraph->GetPointY(i) *
                    fMAttCh1CorrGraph->GetPointY(i))), 2) *
                    fMAttCh1CorrGraph->GetErrorY(i));
        fEnergyRecoCorrGraph->SetPoint(i, positions[i], ecorr / fEref);
        fEnergyRecoCorrGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), ecorr_err / fEref);
    }
    
    TF1* funPol0 = new TF1("funPol0", "pol0", 0, fiberLen);
    fEnergyRecoCorrGraph->Fit(funPol0, "QR+");

    TF1* funEReco = new TF1("funEReco", "sqrt(fun_Sl(x)*fun_Sr(x))*[6]/[7]", 0, fiberLen); 
    funEReco->FixParameter(6, alpha);
    funEReco->FixParameter(7, fEref);
    //----- energy reconstruction end

    //----- setting results start
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyRecoFun, funPol0);
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyRecoGraph, fEnergyRecoCorrGraph);
    
    fResultsExp->AddObject(SFResultTypeObj::kEnergyRecoFun, funEReco);
    fResultsExp->AddObject(SFResultTypeObj::kEnergyRecoGraph, fEnergyRecoGraph);
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

    //----- setting up graphs start
    fEnergyDiffGraph = new TGraphErrors(npointsMax);
    fEnergyDiffGraph->SetName("fEnergyDiffGraph");
    fEnergyDiffGraph->SetTitle("(E_{reco} - E_{ref}) vs. position");
    fEnergyDiffGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyDiffGraph->GetYaxis()->SetTitle("E_{reco}-E_{ref} [keV]");
    fEnergyDiffGraph->SetMarkerStyle(4);
    
    fEnergyDiffCorrGraph = new TGraphErrors(npointsMax);
    fEnergyDiffCorrGraph->SetName("fEnergyDiffCorrGraph");
    fEnergyDiffCorrGraph->SetTitle("(E_{reco corr} - E_{ref}) vs. position");
    fEnergyDiffCorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fEnergyDiffCorrGraph->GetYaxis()->SetTitle("E_{reco corr}-E_{ref} [keV]");
    fEnergyDiffCorrGraph->SetMarkerStyle(8);
    //----- setting up graphs end
    
    //----- energy reconstruction event by event start
    double BL_sigma_cut = SFTools::GetSigmaBL(sipm);
    
    std::vector<SLoop*> trees; 

    double alpha          = fResultsExp->GetValue(SFResultTypeNum::kAlpha);
    double alpha_corr     = fResultsCorr->GetValue(SFResultTypeNum::kAlpha);
    
    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        std::cout << "\t Analyzing position " << positions[npoint] << " mm..." << std::endl;
        
        //----- getting tree
        trees.push_back(fData->GetTree(measurementsIDs[npoint]));
        int        nloopMax = trees[npoint]->getEntries();
        SCategory* tSig     = SCategoryManager::getCategory(SCategory::CatDDSamples);
        
        //----- setting histograms
        TString hname_e = Form("hEnergyRecoExp_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        TString hname_c = Form("hEnergyRecoCorr_S%i_pos%.1f", fSeriesNo, positions[npoint]);
        
        fEnergyRecoHist.push_back(new TH1D(hname_e, hname_e, 500, 0, 1000));
        fEnergyRecoCorrHist.push_back(new TH1D(hname_c, hname_c, 500, 0, 1000));
        
        //----- filling histograms
        for (int nloop = 0; nloop < nloopMax; nloop++)
        {
            trees[npoint]->nextEvent();
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
                    
                    if (t0Ch0 > 0 && t0Ch1 > 0 &&
                        totCh0 > 0 && totCh1 > 0 &&
                        blCh0 < BL_sigma_cut && blCh1 < BL_sigma_cut &&
                        ampCh0 < ampMax && ampCh1 < ampMax &&
                        peCh0 > 0 && peCh1 >0)
                    {
                        double e_reco      = alpha * sqrt(peCh0 * peCh1);
                        double e_corr_reco = alpha_corr * sqrt(fPlRecoFun->Eval(peCh1, peCh0) *
                                                               fPrRecoFun->Eval(peCh1, peCh0));
                        
                        fEnergyRecoHist[npoint]->Fill(e_reco);
                        fEnergyRecoCorrHist[npoint]->Fill(e_corr_reco);
                    }
                    
                }
            }
        }
        
        double mean, sigma;
        
        //TODO fitting of energy-reconstructed spectra
        
        //SFTools::FitGaussSingle(fEnergyRecoHist[npoint], 5);
        //mean  = fEnergyRecoHist[npoint]->GetFunction("fGauss")->GetParameter(1);
        //sigma = fEnergyRecoHist[npoint]->GetFunction("fGauss")->GetParameter(2);
        mean = 1;
        sigma = 0;
        fEnergyDiffGraph->SetPoint(npoint, positions[npoint], mean / fEref);
        fEnergyDiffGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), sigma / fEref);
        
        //SFTools::FitGaussSingle(fEnergyRecoCorrHist[npoint], 5);
        //mean  = fEnergyRecoCorrHist[npoint]->GetFunction("fGauss")->GetParameter(1);
        //sigma = fEnergyRecoCorrHist[npoint]->GetFunction("fGauss")->GetParameter(2);
        fEnergyDiffCorrGraph->SetPoint(npoint, positions[npoint], mean / fEref);
        fEnergyDiffCorrGraph->SetPointError(npoint, SFTools::GetPosError(collimator, testBench), sigma / fEref);
    }
    //----- energy reconstruction end
    
    //----- setting results start
    fResultsCorr->AddObject(SFResultTypeObj::kEnergyDiffGraph, fEnergyDiffCorrGraph);
    
    fResultsExp->AddObject(SFResultTypeObj::kEnergyDiffGraph, fEnergyDiffGraph);
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
std::vector<TH1D*> SFEnergyReco::GetEnergyRecoHistograms(TString type)
{
    if ((type == "experimental" && fEnergyRecoHist.empty()) || 
        (type == "corrected" && fEnergyRecoCorrHist.empty()))
    {
        std::cerr << "##### Error in SFEnergyReco::GetEnergyRecoHistograms(). Empty vector!" << std::endl;
        std::abort();
    }

    std::vector<TH1D*> hist;
    
    if (type == "experimental")
        hist = fEnergyRecoHist;
    else if (type == "corrected")
        hist = fEnergyRecoCorrHist;
    else
    {
        std::cerr << "##### Error in SFEnergyReco::GetEnergyRecoHistograms()!" << std::endl;
        std::cerr << "Incorrect type. Possible options are: experimental and corrected" << std::endl;
        std::abort();
    }
    
    return hist;
}
//------------------------------------------------------------------
