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

ClassImp(SFPeakFinder);

//------------------------------------------------------------------
/// Default constructor.
SFPeakFinder::SFPeakFinder() : fSpectrum(nullptr),
                               fPeak(nullptr),
                               fFittedFun(nullptr),
                               fVerbose(false),
                               fTests(false),
                               fResults(new SFResults("PeakFinderResults_tmp"))
{
    std::cout << "#### Warning in SFPeakFinder constructor!" << std::endl;
    std::cout << "You are using default constructor!" << std::endl;
    fVerbose = false;
}
//------------------------------------------------------------------
/// Standard constructor.
/// \param spectrum - analyzed spectrum
/// \param verbose - print-outs level. Set to true for verbose and false for quiet mode
/// \param tests - flag for testing mode. Set true for tests and false for runtime.
SFPeakFinder::SFPeakFinder(TH1D* spectrum, bool verbose, bool tests) : fSpectrum(spectrum),
                                                                       fPeak(nullptr),
                                                                       fFittedFun(nullptr),
                                                                       fVerbose(verbose),
                                                                       fTests(tests),
                                                                       fResults(new SFResults("PeakFinderResults_tmp"))
{
}
//------------------------------------------------------------------
/// Standard constructor.
/// \param spectrum - analyzed spectrum
/// \param verbose - print-outs level. Set true for verbose and false for quiet mode.
SFPeakFinder::SFPeakFinder(TH1D* spectrum, bool verbose) : fSpectrum(spectrum),
                                                           fPeak(nullptr),
                                                           fFittedFun(nullptr),
                                                           fVerbose(verbose),
                                                           fTests(false),
                                                           fResults(new SFResults("PeakFinderResults_tmp"))
{
}
//------------------------------------------------------------------
/// Standard constructor. Vebose level set by default to quiet. In order to change verbose
/// level use SetVerbLevel().
/// \param spectrum - analyzed spectrum
SFPeakFinder::SFPeakFinder(TH1D* spectrum) : fSpectrum(spectrum),
                                             fPeak(nullptr),
                                             fFittedFun(nullptr),
                                             fVerbose(false),
                                             fTests(false),
                                             fResults(new SFResults("PeakFinderResults_tmp"))
{
    std::cout << "##### Warning in SFPeakFinder constructor. Quiet mode on, no print outs."
              << std::endl;
    std::cout << "##### To set verbose level use SetVerbLevel(bool)" << std::endl;
}
//------------------------------------------------------------------
/// Default destructor.
SFPeakFinder::~SFPeakFinder()
{
}
//------------------------------------------------------------------
/// If default constructor was used, sets analyzed spectrum histogram.
/// \param spectrum - histogram containing analyzed spectrum.
void SFPeakFinder::SetSpectrum(TH1D* spectrum)
{

    if (spectrum == nullptr)
    {
        std::cerr << "##### Error in SFPeakFinder::SetSpectrum()! " << std::endl;
        std::cerr << "Spectrum cannot be a null pointer!" << std::endl;
        std::abort();
    }
    fSpectrum = spectrum;
    return;
}
//------------------------------------------------------------------
TString SFPeakFinder::Init(void)
{

    TString hname    = fSpectrum->GetName();
    int     ID       = SFTools::GetMeasurementID(hname);
    int     seriesNo = SFTools::GetSeriesNo(hname);

    fResults->SetName(Form("PeakFinderResults_%i", ID));

    SFData* data;
    try
    {
        data = new SFData(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Error in SFPeakFinder::Init()!" << std::endl;
    }

    std::vector<TString> names     = data->GetNames();
    std::vector<int>     measureID = data->GetMeasurementsIDs();
    std::vector<double>  positions = data->GetPositions();
    int                  index     = SFTools::GetIndex(measureID, ID);
    TString              dir_name  = names[index];
    TString              full_path = SFTools::FindData(dir_name);

    TString conf_name = "/fitconfig.txt";

    TString functions = "gaus(0) pol0(3)+[4]*TMath::Exp((x-[5])*[6])";
    TString hnames[5] = {Form("S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[index], ID),
                         Form("S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[index], ID),
                         Form("S%i_pos%.1f_ID%i_PEAverage", seriesNo, positions[index], ID),
                         Form("hEnergyRecoExp_S%i_pos%.1f", seriesNo, positions[index]),
                         Form("hEnergyRecoCorr_S%i_pos%.1f", seriesNo, positions[index])};

    std::fstream test(full_path + conf_name, std::ios::in);

    if (test.fail())
    {
        std::cout << "Fitting config for " << full_path << " doesn't exist..." << std::endl;
        std::cout << "Creating new config file..." << std::endl;

        TSpectrum* spec   = new TSpectrum(10);
        int        npeaks = spec->Search(fSpectrum, 10, "goff", 0.1);
        double*    peaksX = spec->GetPositionX();
        double     peak   = TMath::MaxElement(npeaks, peaksX);

        TString opt;
        if (fVerbose)
            opt = "0R";
        else
            opt = "Q0R";

        TF1* fun_gaus = new TF1("fun_gaus", "gaus", peak - 100, peak + 100);
        fSpectrum->Fit("fun_gaus", opt);

        double par0     = fun_gaus->GetParameter(0);
        double par1     = fun_gaus->GetParameter(1);
        double par2     = fun_gaus->GetParameter(2);
        double par2_min = 0;
        double par2_max = 300;
        double par3     = 30.;

        TF1* fun_expo =
            new TF1("fun_expo", "[0]*TMath::Exp((x-[1])*[2])", par1 - 5 * par2, par1 - 4 * par2);
        fSpectrum->Fit("fun_expo", opt);

        double par4 = fun_expo->GetParameter(0);
        double par5 = fun_expo->GetParameter(1);
        double par6 = fun_expo->GetParameter(2);

        double xmin = par1 - par2 * 2;
        double xmax = par1 + par2 * 3;

        std::fstream config(full_path + conf_name, std::ios::out);
        for (int i = 0; i < 3; i++)
        {
            config << " " << hnames[i] << " " << functions << " " << 0 << " " << xmin << " " << xmax
                   << " " << par0 << " " << par1 << " " << par2 << " : " << par2_min << " "
                   << par2_max << " " << par3 << " " << par4 << " " << par5 << " " << par6 << "\n";
        }
        
        par0 = 200;
        par1 = 511;
        par2 = 30;
        par3 = 10;
        par4 = 500;
        par5 = 100;
        par6 = -0.02;
        xmin = par0 - (par2 * 3.5);
        xmax = par0 + (par2 * 4.5);
        
        for (int i = 3; i < 5; i++)
        {
            config << " " << hnames[i] << " " << functions << " " << 0 << " " << xmin << " " << xmax
                   << " " << par0 << " " << par1 << " " << par2 << " : " << par2_min << " "
                   << par2_max << " " << par3 << " " << par4 << " " << par5 << " " << par6 << "\n";
        }
        
        config.close();
    }
    else
    {
        std::cout << "Fitting config for " << full_path << " exists!" << std::endl;
        test.close();
    }

    return full_path;
}
//------------------------------------------------------------------
/// Finds range of the 511 keV peak. Range is returned as references.
/// For measurements with lead collimator range is defined as position
/// +/- sigma and for measurements with electronic collimator range is
/// defined as position +/- 2 sigma.
bool SFPeakFinder::FindPeakRange(double& min, double& max)
{

    // Getting series attributes
    int     seriesNo = SFTools::GetSeriesNo(fSpectrum->GetName());
    SFData* data;

    try
    {
        data = new SFData(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in SFPeakFinder::FindPeakRange()!" << std::endl;
        std::abort();
    }

    TString type = data->GetCollimator();
    delete data;

    // Calculating peak range
    const double delta = 1E-8;

    double pos = fResults->GetValue(SFResultTypeNum::kPeakPosition);
    double sig = fResults->GetValue(SFResultTypeNum::kPeakSigma);

    if (fabs(sig + 1) < delta || fabs(pos + 1) < delta) 
    {
        std::cout << "##### Warning in SFPeakFinder::FindPeakRange(): peak_pos = "
                  << pos << "\t peak_sigma = " << sig << std::endl; 
        std::cout << "Fitting peak... " << std::endl;
        
        FindPeakFit();
        pos = fResults->GetValue(SFResultTypeNum::kPeakPosition);
        sig = fResults->GetValue(SFResultTypeNum::kPeakSigma);
        
        std::cout << "After fitting: peak_pos = " << pos 
                  << "\t peak_sig = " << sig << std::endl;
    }

    if (type == "Lead")
    {
        min = pos - sig;
        max = pos + 1.5 * sig;
    }
    else if (type == "Electronic")
    {
        min = pos - 2 * sig;
        max = pos + 2 * sig;
    }

    if (max < min || fabs(min + 1) < delta || fabs(max + 1) < delta)
    {
        std::cerr << "##### Error in SFPeakFinder::FindPeakRange(). Incorrect range." << std::endl;
        std::cerr << "min = " << min << "\t max = " << max << std::endl;
        return false;
    }

    return true;
}
//------------------------------------------------------------------
/// Fits function describing spectrum in the area where 511 keV peak
/// is visible. Fitted function: expo + pol0 + gaus. Results can be
/// accessed via SFPeakFinder::GetParameters() function and other
/// defined getters. fPeak histogram is not filled.
bool SFPeakFinder::FindPeakFit(void)
{

    TString data_path = Init();

    FitterFactory fitter;
    fitter.initFactoryFromFile((data_path + "/fitconfig.txt").Data(),
                               (data_path + "/fitparams.out").Data());
    HistFitParams             histFP;
    FitterFactory::FIND_FLAGS fl = fitter.findParams(fSpectrum->GetName(), histFP);
    printf("fl = %d for %s\n", fl, fSpectrum->GetName());
    fitter.fit(histFP, fSpectrum);
    fitter.updateParams(fSpectrum, histFP);
    fitter.exportFactoryToFile();
    fFittedFun = (TF1*)histFP.funSum->Clone();
    //fFittedFun->Print();

    if (fFittedFun == nullptr)
    {
        std::cerr << "##### Error in SFPeakFinder::FindPeak()! Function is null pointer"
                  << std::endl;
        std::abort();
    }

    fResults->AddResult(SFResultTypeNum::kPeakConst, fFittedFun->GetParameter(0),
                        fFittedFun->GetParError(0));
    fResults->AddResult(SFResultTypeNum::kPeakPosition, fFittedFun->GetParameter(1),
                        fFittedFun->GetParError(1));
    fResults->AddResult(SFResultTypeNum::kPeakSigma, fFittedFun->GetParameter(2),
                        fFittedFun->GetParError(2));

    // for tests
    if (fTests)
    {
        double fit_min = 0;
        double fit_max = 0;
        fFittedFun->GetRange(fit_min, fit_max);

        TString opt;
        if (fVerbose)
            opt = "BS+";
        else
            opt = "BSQ+";

        TF1* fun_bg_clone =
            new TF1("fun_bg_clone", "pol0(0)+[1]*TMath::Exp((x-[2])*[3])", fit_min, fit_max);
        fun_bg_clone->SetLineColor(kGreen + 3);
        fun_bg_clone->FixParameter(0, fFittedFun->GetParameter(3));
        fun_bg_clone->FixParameter(1, fFittedFun->GetParameter(4));
        fun_bg_clone->FixParameter(2, fFittedFun->GetParameter(5));
        fun_bg_clone->FixParameter(3, fFittedFun->GetParameter(6));
        fSpectrum->Fit("fun_bg_clone", opt);

        if (fVerbose)
        {
            std::cout << "Fitted functions parameters: \n"
                      << fFittedFun->GetParameter(3) << "\t" << fFittedFun->GetParameter(4) << "\t"
                      << fFittedFun->GetParameter(5) << "\t" << fFittedFun->GetParameter(6)
                      << std::endl;

            std::cout << "Expo clone parameters: \n"
                      << fun_bg_clone->GetParameter(0) << "\t" << fun_bg_clone->GetParameter(1)
                      << "\t" << fun_bg_clone->GetParameter(2) << "\t"
                      << fun_bg_clone->GetParameter(3) << std::endl;
        }

        TF1* fun_gaus_clone = new TF1("fun_gaus_clone", "gaus", fit_min, fit_max);
        fun_gaus_clone->SetLineColor(kMagenta);
        fun_gaus_clone->FixParameter(0, fFittedFun->GetParameter(0));
        fun_gaus_clone->FixParameter(1, fFittedFun->GetParameter(1));
        fun_gaus_clone->FixParameter(2, fFittedFun->GetParameter(2));
        fSpectrum->Fit("fun_gaus_clone", opt);

        if (fVerbose)
        {
            std::cout << "Fitted function parameters: \n"
                      << fFittedFun->GetParameter(0) << "\t" << fFittedFun->GetParameter(1) << "\t"
                      << fFittedFun->GetParameter(2) << std::endl;

            std::cout << "Gaus clone parameters: \n"
                      << fun_gaus_clone->GetParameter(0) << "\t" << fun_gaus_clone->GetParameter(1)
                      << "\t" << fun_gaus_clone->GetParameter(2) << std::endl;
        }
    }

    if (fResults->GetValue(SFResultTypeNum::kPeakPosition) < 0 ||
        fResults->GetValue(SFResultTypeNum::kPeakSigma) < 0)
    {
        std::cerr
            << "##### Error in SFPeakFinder::FitPeak(). Position and Sigma cannot be negative!"
            << std::endl;
        std::cerr << "fPosition = " << fResults->GetValue(SFResultTypeNum::kPeakPosition)
                  << "\t fSigma = " << fResults->GetValue(SFResultTypeNum::kPeakSigma) << std::endl;
        return false;
    }

    fResults->AddObject(SFResultTypeObj::kSpectrum, fSpectrum);

    return true;
}
//------------------------------------------------------------------
/// Finds 511 keV peak via SFPeakFinder::FindPeakFit() method and
/// performs background subtraction. Exponential function is fitted on
/// the left sige of the peak and pol0 function - on the right side. For
/// the background subtraction functions are sewed together in the peak
/// region. fPeak histogram is filled with the background-subtracted peak.
/// Subsequently Gaussian function is fitted in order to describe peak shape.
/// Results are accessible via SFPeakFinder::GetParameters() and other
/// definded getter functions.
/// IMPORTANT - if you want to use this function, all the parameters of
/// the fits need to be adjusted!

bool SFPeakFinder::SubtractBackground(void)
{

    // setting background-subtracted histogram
    TString tmp   = fSpectrum->GetName();
    TString pname = tmp.Append("_peak");
    fPeak         = (TH1D*)fSpectrum->Clone(pname);
    fPeak->Reset();

    // getting background + signal function
    if (fFittedFun == nullptr) FindPeakFit();

    double peak_min, peak_max;
    FindPeakRange(peak_min, peak_max);

    peak_min = peak_min - fResults->GetValue(SFResultTypeNum::kPeakSigma);
    peak_max = peak_max + fResults->GetValue(SFResultTypeNum::kPeakSigma);

    TF1* fun_bg = new TF1("fun_bg", "pol0(0)+[1]*TMath::Exp((x-[2])*[3])", 0, 1000);
    fun_bg->SetParameters(fFittedFun->GetParameter(3), fFittedFun->GetParameter(4),
                          fFittedFun->GetParameter(5), fFittedFun->GetParameter(6));

    // background subtraction
    double x, y;

    for (int i = 1; i < fSpectrum->GetNbinsX() + 1; i++)
    {
        x = fSpectrum->GetBinCenter(i);
        if (x > peak_min && x < peak_max) { y = fSpectrum->GetBinContent(i) - fun_bg->Eval(x); }
        else
            y = 0;
        fPeak->SetBinContent(i, y);
    }

    fResults->AddObject(SFResultTypeObj::kPeak, fPeak);

    return true;
}
//------------------------------------------------------------------
/// Prints details of SFPeakFinder class object.
void SFPeakFinder::Print(void)
{
    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "This is print out of SFPeakFinder class object" << std::endl;
    std::cout << "Analyzed spectrum: " << std::endl;
    if (fSpectrum != nullptr)
        fSpectrum->Print();
    else
        std::cout << "NULL" << std::endl;
}
//------------------------------------------------------------------
