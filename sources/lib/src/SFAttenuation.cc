// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFAttenuation.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFAttenuation.hh"

#include <TGraphErrors.h>

#include <memory>

int iparSlAtt[] = {0, 1, 2, 3};
int iparSrAtt[] = {0, 1, 2, 3};

//------------------------------------------------------------------
/// Custom pol3 formula for fitting of MLR dependence.
double myPol3(double* x, double* par)
{
    double f = par[0] + par[1] * (x[0] + par[2] * pow(x[0], 3));
    return f;
}
//------------------------------------------------------------------
/// Method to determine attenuation length used in Pauwels et al., JINST 8 (2013) P09019.
/// For both ends of the fiber one MLR value is calculated, since combined signal from both 
/// channels is taken into account. Obtained MLR dependence is fitted with pol1 or pol3 
/// function.
SFResults* SFAttenuation::AttCombinedCh(TString fun, SFInfo* info,
                                        double pos_uncert,
                                        std::vector<double> positions,
                                        std::vector<TH1D*> spectra)
{
    int npoints = spectra.size();
    TString collimator = info->GetCollimator();
    TString sipm = info->GetSiPM();
    TString testBench = info->GetTestBench();
    
    TString gname = "attenuation_mlr";
    auto mlr = new TGraphErrors(npoints);
    mlr->GetXaxis()->SetTitle("source position [mm]");
    mlr->GetYaxis()->SetTitle("M_{LR}");
    mlr->SetTitle(gname);
    mlr->SetName(gname);
    mlr->SetMarkerStyle(4);

    gname = "sigma_mlr";
    auto sigma = new TGraphErrors(npoints);
    sigma->GetXaxis()->SetTitle("source position [mm]");
    sigma->GetYaxis()->SetTitle("#sigma M_{LR}");
    sigma->SetTitle(gname);
    sigma->SetName(gname);
    sigma->SetMarkerStyle(4);

    int     parNo    = 1;
    double  const_1  = 0;
    double  const_2  = 0;
    TString fun_name = "";

    for (int i = 0; i < npoints; i++)
    {
        if(testBench == "PL")
        {
            if (collimator.Contains("Lead"))
            {
                SFTools::RatiosFitDoubleGauss(spectra, 5);
                fun_name = "fDGauss";
                const_1  = spectra[i]->GetFunction(fun_name)->GetParameter(0);
                const_2  = spectra[i]->GetFunction(fun_name)->GetParameter(3);
                if (const_1 > const_2)
                    parNo = 1;
                else
                    parNo = 4;
            }
            else if (collimator.Contains("Electronic") && sipm.Contains("SensL"))
            {
                SFTools::RatiosFitDoubleGauss(spectra, 2);
                fun_name = "fDGauss";
                const_1  = spectra[i]->GetFunction(fun_name)->GetParameter(0);
                const_2  = spectra[i]->GetFunction(fun_name)->GetParameter(3);
                if (const_1 > const_2)
                    parNo = 1;
                else
                    parNo = 4;
            }
            else if (collimator.Contains("Electronic") && sipm.Contains("Hamamatsu"))
            {
                SFTools::RatiosFitGauss(spectra, 2);
                fun_name = "fGauss";
                parNo    = 1;
            }
        }
        else if (testBench == "PRO KRK")
        {
                SFTools::RatiosFitDoubleGauss(spectra, 5);
                fun_name = "fDGauss";
                const_1  = spectra[i]->GetFunction(fun_name)->GetParameter(0);
                const_2  = spectra[i]->GetFunction(fun_name)->GetParameter(3);
                if (const_1 > const_2)
                    parNo = 1;
                else
                    parNo = 4;
        }
        else if (testBench == "PMI")
        {
            SFTools::RatiosFitDoubleGauss(spectra, 5);   /// TODO tune!
            fun_name = "fDGauss";
            const_1  = spectra[i]->GetFunction(fun_name)->GetParameter(0);
            const_2  = spectra[i]->GetFunction(fun_name)->GetParameter(3);
            if (const_1 > const_2)
                parNo = 1;
            else
                parNo = 4;
        }

        mlr->SetPoint(i, positions[i], spectra[i]->GetFunction(fun_name)->GetParameter(parNo));
        mlr->SetPointError(i, pos_uncert, spectra[i]->GetFunction(fun_name)->GetParError(parNo));
        sigma->SetPoint(i, positions[i], spectra[i]->GetFunction(fun_name)->GetParameter(parNo + 1));
        sigma->SetPointError(i, pos_uncert, spectra[i]->GetFunction(fun_name)->GetParError(parNo + 1));
    }

    double att     = 0;
    double err     = 0;
    double chi2ndf = 0;
    
    if(fun.Contains("pol1"))
    {
        TF1* fpol1 = new TF1("fpol1", "pol1", -50, 150);
        fpol1->SetParameters(-1E-3, -9E-3);
        TFitResultPtr ptr = mlr->Fit(fpol1, "SQR+");

        att = fabs(1. / fpol1->GetParameter(1));
        err = fpol1->GetParError(1) / pow(fpol1->GetParameter(1), 2);
        chi2ndf = ptr->Chi2() / ptr->Ndf();
        
        std::cout << "Attenuation lenght (MLR) from pol1 fit is: " << att << " +/- " 
                  << err << " mm" << std::endl;
        std::cout << "Chi2/NDF = " << chi2ndf << std::endl << std::endl;
    }
    else if(fun.Contains("pol3"))
    {
        TF1* fpol3 = new TF1("fpol3", myPol3, -50, 150, 3);
        fpol3->SetParameters(-1E-3, -9E-3, 50);
        TFitResultPtr ptr = mlr->Fit(fpol3, "SQR+");
        
        att     = fabs(1. / fpol3->GetParameter(1));
        err     = fpol3->GetParError(1) / pow(fpol3->GetParameter(1), 2);
        chi2ndf = ptr->Chi2() / ptr->Ndf();   
        
        std::cout << "Attenuation lenght (MLR) from pol3 fit is: " << att << " +/- " 
                  << err << " mm" << std::endl;
        std::cout << "Chi2/NDF = " << chi2ndf << std::endl << std::endl;
    }
    else
    {
        std::cerr << "##### Error in SFAttenuation::AttCombinedCh()!" << std::endl;
        std::cerr << "Unknown function type. Possible options are: pol1, pol3" << std::endl;
        std::abort();
    }    

    SFResults* results = new SFResults("AttenuationMLR");

    results->AddResult(SFResultTypeNum::kLambda, att, err);
    results->AddResult(SFResultTypeNum::kChi2NDF, chi2ndf, -1);
    results->AddObject(SFResultTypeObj::kAttGraph, mlr);
    results->AddObject(SFResultTypeObj::kMLRSigmaGraph, sigma);

    return results;
}
//------------------------------------------------------------------
/// Method to determine attenuation length for both channels independently.
/// If series was measured with lead collimator peak position is determied
/// with the FindPeakNoBackground() method of the SFPeakFinder class. If
/// series was measured with electronic collimator - FindPeakFit() method
/// of the SFPeakFinder class is used.
/// \param ch - channel number
SFResults* SFAttenuation::AttSeparateCh(char side, double pos_uncert, SFInfo *info,
                                        std::vector<double> positions,
                                        std::vector<TH1D*> spectra, 
                                        std::vector<TString> path)
{
    int npoints = spectra.size();
    double fiberLen = info->GetFiberLength();
    
    TString gname = Form("attenuation_separate_%c", side);
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    graph->SetTitle(gname);
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    for (int i = 0; i < npoints; i++)
    {
        auto peakParams = std::unique_ptr<SFResults>(SFPeakFinder::FindPeakFit(spectra[i], path[i], 0, 0));
        graph->SetPoint(i, positions[i], peakParams->GetValue(SFResultTypeNum::kPeakPosition));
        graph->SetPointError(i, pos_uncert, peakParams->GetUncertainty(SFResultTypeNum::kPeakPosition));
    }

    //----- fitting
    TF1* fexp;

    if (side == 'L' || side == 'l')
    {
        fexp = new TF1("fun_l", "[0]*exp(-x/[1])", positions[0], positions[npoints - 1]);
        fexp->SetParameters(10, 100);
    }
    else if (side == 'R' || side == 'r')
    {
        fexp = new TF1("fun_r", "[0]*exp(-([2]-x)/[1])", positions[0], positions[npoints - 1]);
        fexp->SetParameter(0, 10);
        fexp->SetParameter(1, 100);
        fexp->FixParameter(2, fiberLen);
    }
    else
    {
        std::cerr << "##### Error in SFAttenuation::AttSeparateCh()!" << std::endl;
        std::cerr << "Incorrect side. Possible options: R/r/L/l" << std::endl;
        fexp = nullptr;
        std::abort();
    }

    TFitResultPtr ptr     = graph->Fit(fexp, "QRS");
    double        Chi2NDF = ptr->Chi2() / ptr->Ndf();

    std::cout << "\n\tAttenuation for side " << side << ": " << fexp->GetParameter(1) << " +/- "
              << fexp->GetParError(1) << " mm" << std::endl;
    std::cout << "Chi2/NDF = " << Chi2NDF << std::endl << std::endl;

    SFResults* results = new SFResults(Form("AttenuationSeparateCh_%c", side));

    results->AddResult(SFResultTypeNum::kLambda, fexp->GetParameter(1), fexp->GetParError(1));
    results->AddResult(SFResultTypeNum::kChi2NDF, Chi2NDF, -1);
    results->AddObject(SFResultTypeObj::kAttGraph, graph);

    return results;
}
//------------------------------------------------------------------
/// Method for simultaneous fitting of data from both channels. Simple 
/// exponential attenuation model is fitted.
SFResults* SFAttenuation::AttSimultaneousFit(TGraphErrors *attL, TGraphErrors *attR)
{
    double fiberLen = 100;
    
    TFormula* form_Sl = new TFormula("Sl", "[0]*exp(-x/[2])");
    TFormula* form_Sr = new TFormula("Sr", "[1]*exp(-([3]-x)/[2])");

    TF1* fun_Sl = new TF1("fun_Sl", "Sl(x)", 0, fiberLen);
    TF1* fun_Sr = new TF1("fun_Sr", "Sr(x)", 0, fiberLen);
    
    //----- fitting start -----
    ROOT::Math::WrappedMultiTF1 wfSl(*fun_Sl, 1);
    ROOT::Math::WrappedMultiTF1 wfSr(*fun_Sr, 1);

    ROOT::Fit::DataOptions opt;

    ROOT::Fit::DataRange rangeSl;
    rangeSl.SetRange(0, fiberLen);
    ROOT::Fit::BinData dataSl(opt, rangeSl);
    ROOT::Fit::FillData(dataSl, attL);
    
    ROOT::Fit::DataRange rangeSr;
    rangeSr.SetRange(0, fiberLen);
    ROOT::Fit::BinData dataSr(opt, rangeSr);
    ROOT::Fit::FillData(dataSr, attR);
    
    ROOT::Fit::Chi2Function chi2_Sl(dataSl, wfSl);
    ROOT::Fit::Chi2Function chi2_Sr(dataSr, wfSr);

    GlobalChi2Att globalChi2(chi2_Sl, chi2_Sr);

    ROOT::Fit::Fitter fitter;

    const int npar       = 4;
    double    par0[npar] = {1000, 1000, 500, fiberLen};

    fitter.Config().SetParamsSettings(npar, par0);

    fitter.Config().ParSettings(0).SetLimits(0, 1E6);
    fitter.Config().ParSettings(1).SetLimits(0, 1E6);
    fitter.Config().ParSettings(2).SetLimits(0, 1E4);
    fitter.Config().ParSettings(3).Fix();

    fitter.Config().ParSettings(0).SetName("S0");
    fitter.Config().ParSettings(1).SetName("S1");
    fitter.Config().ParSettings(2).SetName("Latt");
    fitter.Config().ParSettings(3).SetName("L");

    fitter.Config().MinimizerOptions().SetPrintLevel(1);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");
    fitter.Config().MinimizerOptions().SetMaxIterations(1e6);
    fitter.Config().MinimizerOptions().SetMaxFunctionCalls(1e6);

    fitter.FitFCN(npar, globalChi2, 0, dataSl.Size() + dataSr.Size(), true);
    
    ROOT::Fit::FitResult* fitterResults = new ROOT::Fit::FitResult(fitter.Result());
    const double* params = fitterResults->GetParams();
    const double* errors = fitterResults->GetErrors();
    
    fun_Sl->SetFitResult(*fitterResults, iparSlAtt);
    fun_Sl->SetRange(rangeSl().first, rangeSl().second);

    fun_Sr->SetFitResult(*fitterResults, iparSrAtt);
    fun_Sr->SetRange(rangeSr().first, rangeSr().second);
    
    std::cout << "\n Attenuation length from simultaneuus fit: " 
              << params[2] << " +/- " << errors[2] << "\n" << std::endl;
    //----- fitting end -----
    
    SFResults *results = new SFResults("AttenuationSimultaneusFit"); 
              
    results->AddResult(SFResultTypeNum::kLambda, params[2], errors[2]); 
    results->AddResult(SFResultTypeNum::kChi2NDF, fitterResults->Chi2() / fitterResults->Ndf(), -1);
    results->AddObject(SFResultTypeObj::kSlFun, fun_Sl);
    results->AddObject(SFResultTypeObj::kSrFun, fun_Sr);
     
    return results;
}
