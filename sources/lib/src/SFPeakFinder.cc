// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFPeakFinder.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFPeakFinder.hh"

//------------------------------------------------------------------
/// Fits function describing spectrum in the area where 511 keV peak
/// is visible. Fitted function: expo + pol0 + gaus. Results can be
/// accessed via SFPeakFinder::GetParameters() function and other
/// defined getters. fPeak histogram is not filled.
SFResults* SFPeakFinder::FindPeakFit(TH1D* spectrum, TString path, bool verbose, bool tests)
{
    
    FitterFactory fitter;
    fitter.initFactoryFromFile((path + "/fitconfig.txt").Data(),
                               (path + "/fitparams.out").Data());
    
    auto name = spectrum->GetName();
    auto* hfp = fitter.findParams(name);
    
    if(!hfp)
    {
        std::cerr << "No fit parameters for this histogram in fitparams.out!" << std::endl;
        std::cerr << "Please do initial fitting first!" << std::endl;
        std::abort();
    }

    auto histFP = fitter.findParams(name);
    histFP->print();
    fitter.fit(histFP, spectrum);
    fitter.exportFactoryToFile();
    
    TF1* fitted = (TF1*)histFP->function_sum.Clone();
    
    if (fitted == nullptr)
    {
        std::cerr << "##### Error in SFPeakFinder::" << __func__ 
                  << " Function is null pointer" << std::endl;
        std::abort();
    }

    TF1 *tmpfun = spectrum->GetFunction(fitted->GetName());
    double chi2NDF = tmpfun->GetChisquare() / tmpfun->GetNDF();
    
    SFResults* results = new SFResults(Form("PeakFinder_%s", TString(spectrum->GetName()).Data()));
    
    results->AddResult(SFResultTypeNum::kPeakConst, 
                        fitted->GetParameter(0),
                        fitted->GetParError(0));
    results->AddResult(SFResultTypeNum::kPeakPosition, 
                        fitted->GetParameter(1),
                        fitted->GetParError(1));
    results->AddResult(SFResultTypeNum::kPeakSigma, 
                        fitted->GetParameter(2),
                        fitted->GetParError(2));
    results->AddResult(SFResultTypeNum::kChi2NDF, chi2NDF, -1);

    // for tests
    if (tests)
    {
        double fit_min = 0;
        double fit_max = 0;
        fitted->GetRange(fit_min, fit_max);

        TString opt;
        if (verbose)
            opt = "BS+";
        else
            opt = "BSQ+";

        TF1* fun_bg_clone = new TF1("fun_bg_clone", "pol0(0)+[1]*TMath::Exp((x-[2])*[3])", fit_min, fit_max);
        fun_bg_clone->SetLineColor(kGreen + 3);
        fun_bg_clone->FixParameter(0, fitted->GetParameter(3));
        fun_bg_clone->FixParameter(1, fitted->GetParameter(4));
        fun_bg_clone->FixParameter(2, fitted->GetParameter(5));
        fun_bg_clone->FixParameter(3, fitted->GetParameter(6));
        spectrum->Fit("fun_bg_clone", opt);

        if (verbose)
        {
            std::cout << "Fitted functions parameters: \n"
                      << fitted->GetParameter(3) << "\t" 
                      << fitted->GetParameter(4) << "\t"
                      << fitted->GetParameter(5) << "\t" 
                      << fitted->GetParameter(6)
                      << std::endl;

            std::cout << "Expo clone parameters: \n"
                      << fun_bg_clone->GetParameter(0) << "\t" 
                      << fun_bg_clone->GetParameter(1) << "\t" 
                      << fun_bg_clone->GetParameter(2) << "\t"
                      << fun_bg_clone->GetParameter(3) << std::endl;
        }

        TF1* fun_gaus_clone = new TF1("fun_gaus_clone", "gaus", fit_min, fit_max);
        fun_gaus_clone->SetLineColor(kMagenta);
        fun_gaus_clone->FixParameter(0, fitted->GetParameter(0));
        fun_gaus_clone->FixParameter(1, fitted->GetParameter(1));
        fun_gaus_clone->FixParameter(2, fitted->GetParameter(2));
        spectrum->Fit("fun_gaus_clone", opt);

        if (verbose)
        {
            std::cout << "Fitted function parameters: \n"
                      << fitted->GetParameter(0) << "\t" 
                      << fitted->GetParameter(1) << "\t"
                      << fitted->GetParameter(2) << std::endl;

            std::cout << "Gaus clone parameters: \n"
                      << fun_gaus_clone->GetParameter(0) << "\t" 
                      << fun_gaus_clone->GetParameter(1) << "\t" 
                      << fun_gaus_clone->GetParameter(2) << std::endl;
        }
    }

    if (results->GetValue(SFResultTypeNum::kPeakPosition) < 0)
    {
        std::cerr << "##### Error in SFPeakFinder::" << __func__ << 
                     " Position cannot be negative!" << std::endl;
        std::cerr << "fPosition = " << results->GetValue(SFResultTypeNum::kPeakPosition) << std::endl;
        std::abort();
    }

    results->AddObject(SFResultTypeObj::kSpectrum, spectrum);

    return results;
}
//------------------------------------------------------------------
/// Finds range of the 511 keV peak. Range is returned as references.
/// For measurements with lead collimator range is defined as position
/// +/- 1.5 sigma and for measurements with electronic collimator range 
/// is defined as position +/- 2 sigma.
auto SFPeakFinder::FindPeakRange(TH1D* spectrum, TString path, TString colimator, 
                                 bool verbose, bool tests)
                                 -> std::tuple<double, double>
{

    TF1* fitted = spectrum->GetFunction(Form("f_%s" ,TString(spectrum->GetName()).Data()));
    
    double pos = -1;
    double sig = -1;
    
    if (fitted == nullptr)
    {
        std::cout << "##### Warning in SFPeakFinder::" << __func__ << ": peak_pos = "
                  << pos << "\t peak_sigma = " << sig << std::endl; 
        std::cout << "Fitting peak... " << std::endl;
        
        SFResults *fitting_results = FindPeakFit(spectrum, path, verbose, tests);
        fitted = spectrum->GetFunction(Form("f_%s" ,TString(spectrum->GetName()).Data()));
        pos = fitting_results->GetValue(SFResultTypeNum::kPeakPosition);
        sig = fabs(fitting_results->GetValue(SFResultTypeNum::kPeakSigma));
        
        std::cout << "After fitting: peak_pos = " << pos 
                  << "\t peak_sig = " << sig << std::endl;
    }
    else
    {
        pos = fitted->GetParameter(1);
        sig = fabs(fitted->GetParameter(2));
        std::cout << "In SFPeakFinder::" << __func__ << ": spectrum already fitted." << std::endl;
        std::cout << "peak_pos = " << pos << "\t peak_sigma = " << sig << std::endl;
    }
    
    // Calculating peak range
    double min = -1;
    double max = -1;
    
    if (colimator == "Lead")
    {
        min = pos - sig;
        max = pos + 1.5 * sig;
    }
    else if (colimator == "Electronic" )
    {
        min = pos - 2 * sig;
        max = pos + 2 * sig;
    }
    else if (colimator == "Electronic long")
    {
        min = pos - 2 * sig;
        max = pos + 2 * sig;
    }
    else
    {
        std::cerr << "##### Error in SFPeakFinder::" << __func__ << ". Unknown colimator type" << std::endl;
        std::cerr << "Colimator: " << colimator << ". Posiible options are: Lead, Electronic" << std::endl;
        std::abort();
    }
    
    const double delta = 1E-8;
    
    if (max < min || fabs(min + 1) < delta || fabs(max + 1) < delta)
    {
        std::cerr << "##### Error in SFPeakFinder::" << __func__ << ". Incorrect range." << std::endl;
        std::cerr << "min = " << min << "\t max = " << max << std::endl;
        std::abort();
    }

    return {min, max};
}
//------------------------------------------------------------------
SFResults* SFPeakFinder::SubtractBackground(TH1D* spectrum, TString path, TString colimator, 
                                            bool verbose, bool tests)
{

    // setting background-subtracted histogram
    TString tmp   = spectrum->GetName();
    TString pname = tmp.Append("_peak");
    TH1D* peak         = (TH1D*)spectrum->Clone(pname);
    peak->Reset();

    SFResults* results = new SFResults(Form("PeakFinderBGS_%s", TString(spectrum->GetName()).Data()));
    
    // getting background + signal function
    TF1* fitted = spectrum->GetFunction(Form("f_%s" ,TString(spectrum->GetName()).Data()));
    
    if (fitted == nullptr)
    {
        SFResults *fitting_results = FindPeakFit(spectrum, path, verbose, tests);
        fitted = spectrum->GetFunction(Form("f_%s" ,TString(spectrum->GetName()).Data()));
        results->AddResult(SFResultTypeNum::kPeakConst, 
                           fitting_results->GetValue(SFResultTypeNum::kPeakConst),
                           fitting_results->GetUncertainty(SFResultTypeNum::kPeakConst));
        results->AddResult(SFResultTypeNum::kPeakPosition, 
                           fitting_results->GetValue(SFResultTypeNum::kPeakPosition),
                           fitting_results->GetUncertainty(SFResultTypeNum::kPeakPosition));
        results->AddResult(SFResultTypeNum::kPeakSigma, 
                           fitting_results->GetValue(SFResultTypeNum::kPeakSigma),
                           fitting_results->GetUncertainty(SFResultTypeNum::kPeakSigma));
        results->AddObject(SFResultTypeObj::kSpectrum, spectrum);
    }
    
    auto tuple = FindPeakRange(spectrum, path, colimator, verbose, tests);
    double peak_min = std::get<0>(tuple);
    double peak_max = std::get<1>(tuple);

    peak_min = peak_min - results->GetValue(SFResultTypeNum::kPeakSigma);
    peak_max = peak_max + results->GetValue(SFResultTypeNum::kPeakSigma);

    TF1* fun_bg = new TF1("fun_bg", "pol0(0)+[1]*TMath::Exp((x-[2])*[3])", 0, 5000);
    fun_bg->SetParameters(fitted->GetParameter(3), fitted->GetParameter(4),
                          fitted->GetParameter(5), fitted->GetParameter(6));

    // background subtraction
    double x, y;

    for (int i = 1; i < spectrum->GetNbinsX() + 1; i++)
    {
        x = spectrum->GetBinCenter(i);
        if (x > peak_min && x < peak_max)
        { 
            y = spectrum->GetBinContent(i) - fun_bg->Eval(x);
        }
        else
            y = 0;
        peak->SetBinContent(i, y);
    }

    results->AddObject(SFResultTypeObj::kPeak, peak);

    return results;
}
//------------------------------------------------------------------
