// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPositionReco.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFPositionReco.hh"

const double ampMax = 660;

//------------------------------------------------------------------
auto SFPositionReco::ReconstructPosition(SDDSamples* samples, TF1* fun_mlr, 
                                        float BL_sigma_cut, float qmin, float qmax)
                                        -> std::optional<double>
{
    double T0_l   = samples->getSignalL()->GetT0();
    double T0_r   = samples->getSignalR()->GetT0();
    double PE_l   = samples->getSignalL()->GetPE();
    double PE_r   = samples->getSignalR()->GetPE();
//     double BL_l   = samples->getSignalL()->GetBLSigma();
//     double BL_r   = samples->getSignalR()->GetBLSigma();
    double TOT_l  = samples->getSignalL()->GetTOT();
    double TOT_r  = samples->getSignalR()->GetTOT();
    double amp_l  = samples->getSignalL()->GetAmplitude();
    double amp_r  = samples->getSignalR()->GetAmplitude();
//     bool   veto_l = samples->getSignalL()->GetVeto();
//     bool   veto_r = samples->getSignalR()->GetVeto();

    if (T0_l > 0 && T0_r > 0 &&
        TOT_l > 0 && TOT_r > 0 &&
        /*BL_l < BL_sigma_cut && BL_r < BL_sigma_cut &&*/
        amp_l < ampMax && amp_r < ampMax && 
        sqrt(PE_l * PE_r) > qmin && sqrt(PE_l * PE_r) < qmax /*&&
        veto_l == 0 && veto_r == 0*/)
    {
        double MLR = log(sqrt(PE_r / PE_l));
        return fun_mlr->Eval(MLR);
    }

    return std::nullopt;
}

//------------------------------------------------------------------
SFResults* SFPositionReco::ReconstructPositionDist(SFChAddr addr,
                                                   TH1D* spectrum_av,
                                                   TF1* fun_mlr,
                                                   TString path,
                                                   double BL_sigma_cut,
                                                   TString collimator)
{
 
    std::cout << "Reconstructing: " << path << std::endl;
    
    TString full_tree_path = path + "/sifi_results.root";
    std::string full_tree_path_str = std::string(full_tree_path);
    std::cout << "\nOpening file: " << full_tree_path << std::endl;
    
    auto tfile = std::make_unique<TFile>(full_tree_path, "READ");
        
    if (!tfile.get()->IsOpen())
    {
        std::cerr << "##### Error in SFPositionReco:: " << __func__ << std::endl;
        std::cerr << "Could not open file: " << std::endl;
        std::cerr << full_tree_path << std::endl;
        std::abort();
    }
    tfile->Close();
        
    SLoop *loop = new SLoop();
    loop->addFile(full_tree_path_str);
    loop->setInput({});
    
    int nloopMax = loop->getEntries();
    auto tSig = SCategoryManager::getCategory(SCategory::CatDDSamples);
    
    //----- setting energy cut
    auto peak_range = SFPeakFinder::FindPeakRange(spectrum_av, path, collimator, 0, 0);
    double xmin = std::get<0>(peak_range);
    double xmax = std::get<1>(peak_range);
    
    TH1D* hist = new TH1D("htemp", "htemp", 200, -200, 200);
    hist->GetXaxis()->SetTitle("reconstructed position [mm]");
    hist->GetYaxis()->SetTitle("counts");

    for (int nloop = 0; nloop < nloopMax; ++nloop)
    {
        loop->getEvent(nloop);
        size_t tentriesMax = tSig->getEntries();

        for (int tentries = 0; tentries < tentriesMax; ++tentries)
        {
            int m, l, f;
            auto samples = (SDDSamples*)tSig->getObject(tentries);
            samples->getAddress(m, l, f);
            
            if (addr.fModule == m &&
                addr.fLayer == l &&
                addr.fFiber == f)
            {
                if (auto val = ReconstructPosition(samples, fun_mlr, BL_sigma_cut, xmin, xmax); val)
                {
                    hist->Fill(*val);
                }
                else
                    continue;
             }
        }
    }
    
    double fit_min = hist->GetBinCenter(2);
    TF1* fun_gauss = new TF1("fun_gaus", "gaus", fit_min, 200);
    fun_gauss->SetParameters(hist->GetBinContent(hist->GetMaximumBin()),
                             hist->GetMean(), hist->GetRMS());
    
    std::cout << "Fit starting parameters: " << hist->GetBinContent(hist->GetMaximumBin()) 
              << "\t" << hist->GetMean() << "\t" << hist->GetRMS() << std::endl;
    
    hist->Fit(fun_gauss, "R");
    
    double sigma_to_fwhm = 2 * sqrt(2 * log(2));
    double pos_res = sigma_to_fwhm * fun_gauss->GetParameter(2);
    double pos_res_err = sigma_to_fwhm * fun_gauss->GetParError(2);
    double pos_reco = fun_gauss->GetParameter(1);
    double pos_reco_err = fun_gauss->GetParError(1);
    
    std::cout << "Position resolution: " << pos_res << " +/- " << pos_res_err  << std::endl;
    std::cout << "Reconstructed position: " << pos_reco << " +/- " << pos_reco_err << std::endl; 
    
    SFResults* results = new SFResults("PositionReconstructionDistSinglePos");
    results->AddObject(SFResultTypeObj::kPositionDist, hist);
    results->AddResult(SFResultTypeNum::kPositionRes, pos_res, pos_res_err);
    results->AddResult(SFResultTypeNum::kPositionReco, pos_reco, pos_reco_err);
    
    return results;
}
//------------------------------------------------------------------
auto SFPositionReco::ReconstructPositionDistAll(SFChAddr addr, 
                                                std::vector<TH1D*> spectrum_av, 
                                                std::vector<double> positions,
                                                std::vector<TString> path,
                                                TF1* fun_mlr, 
                                                double BL_sigma_cut, double pos_uncert, 
                                                TString collimator) 
                                                -> std::tuple<SFResults*,std::vector<TH1D*>>
{
    
    if (fun_mlr == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco:: " << __func__ << std::endl;
        std::cerr << "One of the required pointers is invalid! Please check!" << std::endl;
        std::abort();
    }
    
    std::cout << std::endl;
    fun_mlr->Print();
    std::cout << "fMLR params: " << fun_mlr->GetParameter(0) << " +/- " << fun_mlr->GetParError(0) << "\n" << fun_mlr->GetParameter(1) << " +/- " << fun_mlr->GetParError(1) << std::endl;
    
    if (spectrum_av.empty() || spectrum_av[0] == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco:: " << __func__ << std::endl;
        std::cerr << "Something wrong with std::vector<TH1D*> spectrum_av! Please check!"<< std::endl;
        std::abort();
    }
    
    TString suffix = Form("M%iL%iF%i", addr.fModule, addr.fLayer, addr.fFiber);
    
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
    
    TGraphErrors *gResiduals = new TGraphErrors(npointsMax);
    gResiduals->SetName(Form("PosRecoResiduals_%s", suffix.Data()));
    gResiduals->SetTitle(Form("Reconstructed Position Residuals (%s)", suffix.Data()));
    gResiduals->GetXaxis()->SetTitle("source position [mm]");
    gResiduals->GetYaxis()->SetTitle("residual [mm]");
    gResiduals->SetMarkerStyle(4);
    
    std::vector<TH1D*> pos_distributions;
    
    double posResAv = 0;
    double posResAvErr = 0;
    
    for (int i = 0; i < npointsMax; i++)
    {
        std::cout << "\n\nMeasurement: " << i << std::endl;
        auto results_single_pos = std::unique_ptr<SFResults>(ReconstructPositionDist(addr, 
                                                             spectrum_av[i], fun_mlr, path[i],
                                                             BL_sigma_cut, collimator));
        pos_distributions.push_back((TH1D*)results_single_pos->GetObject(SFResultTypeObj::kPositionDist));
        
        pos_distributions[i]->SetName(Form("hPosReco_%i", i));
        pos_distributions[i]->SetTitle(Form("Reconstructed Position Distribution (M%iL%iF%i) position %.1f mm", addr.fModule, addr.fLayer, addr.fFiber, positions[i]));
        
        double pos_res = results_single_pos->GetValue(SFResultTypeNum::kPositionRes);
        double pos_res_err = results_single_pos->GetUncertainty(SFResultTypeNum::kPositionRes);
        
        double pos_reco = results_single_pos->GetValue(SFResultTypeNum::kPositionReco);
        double pos_reco_err = results_single_pos->GetUncertainty(SFResultTypeNum::kPositionReco);
        
        gPosRecoVsPos->SetPoint(i, positions[i], pos_reco);
        gPosRecoVsPos->SetPointError(i, pos_uncert, pos_reco_err);
        
        gPosResVsPos->SetPoint(i, positions[i], pos_res);
        gPosResVsPos->SetPointError(i, pos_uncert, pos_res_err);
        
        gPosRecoDiff->SetPoint(i, positions[i], (pos_reco - positions[i]));
        gPosRecoDiff->SetPointError(i, pos_uncert, pos_res_err);
        
        posResAv += pos_res * (1. / pow(pos_res_err, 2));
        posResAvErr += (1. / pow(pos_res_err, 2));
    }
    
    posResAv    = posResAv / posResAvErr;
    posResAvErr = sqrt(1. / posResAvErr);
    
    TF1* fun_pol1 = new TF1("fun_pol1", "pol1", 0, 100);
    gPosRecoVsPos->Fit(fun_pol1, "Q");
    
    double residual;
    double point_x, point_y;

    for (int npoint = 0; npoint < npointsMax; npoint++)
    {
        gPosRecoVsPos->GetPoint(npoint, point_x, point_y);
        residual = point_y - positions[npoint];
        gResiduals->SetPoint(npoint, positions[npoint], residual);
    }
    
    std::cout << "\nAverage position resolution for this series is (" << suffix << ") : ";
    std::cout << posResAv << " +/- " << posResAvErr << " mm\n\n" << std::endl;
    
    SFResults* results = new SFResults(Form("PositionResolution_%s", suffix.Data()));
    results->AddResult(SFResultTypeNum::kPositionRes, posResAv, posResAvErr);
    results->AddObject(SFResultTypeObj::kPosRecoVsPosGraph, gPosRecoVsPos);
    results->AddObject(SFResultTypeObj::kPosResVsPosGraph, gPosResVsPos);
    results->AddObject(SFResultTypeObj::kPosDiffVsPosGraph, gPosRecoDiff);
    results->AddObject(SFResultTypeObj::kResidualGraph, gResiduals);
    
    return {results, pos_distributions};
}

//------------------------------------------------------------------
SFResults* SFPositionReco::ReconstructPositionDistSum(SFChAddr addr, 
                                                      std::vector<TH1D*> spectrum_av,
                                                      std::vector<double> positions,
                                                      std::vector<TString> path,
                                                      TF1* fun_mlr,
                                                      double BL_sigma_cut, 
                                                      TString collimator)
{
    
    if (fun_mlr == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco:: " << __func__ << std::endl;
        std::cerr << "One of the required pointers is invalid! Please check!" << std::endl;
        std::abort();
    }
    
    if (spectrum_av.empty() || spectrum_av[0] == nullptr)
    {
        std::cerr << "##### Error in SFPositionReco:: " << __func__ << std::endl;
        std::cerr << "Something wrong with std::vector<TH1D*> spectrum_av! Please check!"<< std::endl;
        std::abort();
    }
    
    int npointsMax = positions.size();
    
    TString suffix = Form("M%iL%iF%i", addr.fModule, addr.fLayer, addr.fFiber);
    
    TString hname = Form("Summed Reconstructed Position Distribution (%s)", suffix.Data());
    TH1D *hPosRecoSum = new TH1D("hPosRecoSum", hname, 300, -300, 300);
    hPosRecoSum->GetXaxis()->SetTitle("reconstructed position - source position [mm]");
    hPosRecoSum->GetYaxis()->SetTitle("counts");
    
    for (int i = 0; i < npointsMax; i++)
    {
        TString full_tree_path = path[i] + "/sifi_results.root";
        std::string full_tree_path_str = std::string(full_tree_path);
        std::cout << "\nOpening file: " << full_tree_path << std::endl;
    
        auto tfile = std::make_unique<TFile>(full_tree_path, "READ");
        
        if (!tfile.get()->IsOpen())
        {
            std::cerr << "##### Error in SFPositionReco::" << __func__ << std::endl;
            std::cerr << "Could not open file:" << std::endl;
            std::cerr << full_tree_path << std::endl;
            std::abort();
        }
        
        tfile->Close();
        
        SLoop *loop = new SLoop();
        loop->addFile(full_tree_path_str);
        loop->setInput({});
        
        int nloopMax = loop->getEntries();
        auto tSig = SCategoryManager::getCategory(SCategory::CatDDSamples);
    
        auto peak_range = SFPeakFinder::FindPeakRange(spectrum_av[i], path[i], collimator, 0, 0);
        double xmin = std::get<0>(peak_range);
        double xmax = std::get<1>(peak_range);
        
        for (int nloop = 0; nloop < nloopMax; ++nloop)
        {
            loop->getEvent(nloop);
            size_t tentriesMax = tSig->getEntries();

            for (int tentries = 0; tentries < tentriesMax; ++tentries)
            {
                int m, l, f;
                auto samples = (SDDSamples*)tSig->getObject(tentries);
                samples->getAddress(m, l, f);

                if (addr.fModule == m &&
                    addr.fLayer == l &&
                    addr.fFiber == f)
                {
                    if (auto val = ReconstructPosition(samples, fun_mlr, BL_sigma_cut, xmin, xmax); val)
                        hPosRecoSum->Fill(*val-positions[i]);
                    else
                        continue;
                }
            }
        }
    }
    
    double fit_min = hPosRecoSum->GetBinCenter(2);
    TF1* fun_gauss = new TF1("fun_gaus", "gaus", fit_min, 200);
    fun_gauss->SetParameters(hPosRecoSum->GetBinContent(hPosRecoSum->GetMaximumBin()),
                             hPosRecoSum->GetMean(), hPosRecoSum->GetRMS());
    
    hPosRecoSum->Fit(fun_gauss, "QR");

    double sigma_to_fwhm = 2 * sqrt(2 * log(2));
    double pos_res = sigma_to_fwhm * fun_gauss->GetParameter(2);
    double pos_res_err = sigma_to_fwhm * fun_gauss->GetParError(2);
    
    SFResults* results = new SFResults(Form("PositionResolutionSum_%s", suffix.Data()));
    results->AddObject(SFResultTypeObj::kPositionDist, hPosRecoSum);
    results->AddResult(SFResultTypeNum::kPositionRes, pos_res, pos_res_err);
    
    return results;
}
//------------------------------------------------------------------
