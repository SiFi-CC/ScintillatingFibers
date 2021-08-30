// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            SFTimeConst.cc             *
// *       J. Kasper, K. Rusiecka          *
// *     kasper@physik.rwth-aachen.de      *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "SFTimeConst.hh"

ClassImp(SFTimeConst);

//------------------------------------------------------------------
/// Default constructor.
SFTimeConst::SFTimeConst() : fSeriesNo(-1),
                             fData(nullptr),
                             fPE(-1),
                             fVerb(false),
                             fResults(nullptr)
{

    std::cout << "##### Warning in SFTimeConst constructor! You are using default constructor!" << std::endl;
    std::cout << "Set object attributes via SetDetails()" << std::endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
/// \param seriesNo - number of analyzed series
/// \param PE - PE value for signals selection
/// \param verb - verbose level
SFTimeConst::SFTimeConst(int seriesNo, double PE, bool verb) : fSeriesNo(seriesNo),
                                                               fPE(PE),
                                                               fVerb(verb),
                                                               fResults(nullptr)
{

    bool stat = SetDetails(seriesNo, PE, verb);
    if (stat == false) { throw "##### Exception in SFTimeConst constructor!"; }
}
//------------------------------------------------------------------
/// Default destructor.
SFTimeConst::~SFTimeConst()
{
    if (fData != nullptr) delete fData;
}
//------------------------------------------------------------------
/// Sets values to private members of the class. Loads TProfile histograms
/// of average signals for requested PE vlaue. Creates vectors of SFFitResults
/// objects.
/// \param seriesNo - number of the series
/// \param PE - PE value
/// \param verb - verbose level
bool SFTimeConst::SetDetails(int seriesNo, double PE, bool verb)
{

    if (fSeriesNo == -1)
    {
        fSeriesNo = seriesNo;
        fPE       = PE;
        fVerb     = verb;
    }

    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cout << message << std::endl;
        throw "##### Exception in SFTimeConst::SetDetails()!";
        return false;
    }

    int              npoints         = fData->GetNpoints();
    TString          fiber           = fData->GetFiber();
    std::vector<int> measurementsIDs = fData->GetMeasurementsIDs();
    TString          selection       = Form("ch_0.fPE>%.1f && ch_0.fPE<%.1f", fPE - 0.5, fPE + 0.5);
    TString          results_name;

    int nsig = 0;
    if (fiber.Contains("LuAG"))
        nsig = 50;
    else if (fiber.Contains("LYSO"))
        nsig = 10;
    else if (fiber.Contains("GAGG"))
        nsig = 20;
    else
    {
        std::cerr << "##### Error in SFTimeConst::SetDetails()! Unknown fiber material!"
                  << std::endl;
        return false;
    }

    for (int i = 0; i < npoints; i++)
    {
        fSignalsCh0.push_back(fData->GetSignalAverage(0, measurementsIDs[i], selection, nsig, true));
        results_name = fSignalsCh0[i]->GetName();
        fFitResultsCh0.push_back(new SFFitResults(results_name));

        fSignalsCh1.push_back(fData->GetSignalAverage(1, measurementsIDs[i], selection, nsig, true));
        results_name = fSignalsCh1[i]->GetName();
        fFitResultsCh1.push_back(new SFFitResults(results_name));
    }
    
    fResults = new SFResults(Form("TimeConstResults_S%i_PE%.2f", fSeriesNo, PE));

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);

    return true;
}
//------------------------------------------------------------------
/// Double decay function describing falling slope of the signal.
double funDecayDouble(double* x, double* par)
{
    double fast_dec = par[0] * TMath::Exp(-(x[0] - par[1]) / par[2]);
    double slow_dec = par[3] * TMath::Exp(-(x[0] - par[1]) / par[4]);
    double constant = par[5];
    return fast_dec + slow_dec + constant;
}
//------------------------------------------------------------------
/// Single decay function describing falling slope of the signal.
double funDecaySingle(double* x, double* par)
{
    double dec      = par[0] * TMath::Exp(-(x[0] - par[1]) / par[2]);
    double constant = par[3];
    return dec + constant;
}
//------------------------------------------------------------------
/// This function performs fitting to the given TProfile signal.
/// Single decay function if fitted to the falling slope of the signal.
/// Results of the fit are subsequently written in the SFFitResults
/// class object. Function returns true if fitting was successful and
/// fit results are valid.
/// \param signal - averaged signal (TProfile)
/// \param ID - measurement ID
bool SFTimeConst::FitDecayTimeSingle(TProfile* signal, int ID)
{

    TString opt;
    if (fVerb)
        opt = "R0";
    else
        opt = "QR0";

    int ch = -1;
    TString hname = signal->GetName();
    if (hname.Contains("ch0"))
        ch = 0;
    else if (hname.Contains("ch1"))
        ch = 1;
    else
    {
        std::cerr << "Error in SFTimeConst::FitDecayTimeSingle()!" << std::endl;
        std::cerr << "Could not interpret signal name!" << std::endl;
        return false;
    }

    std::vector<int> measurementsIDs = fData->GetMeasurementsIDs();
    int    index = SFTools::GetIndex(measurementsIDs, ID);
    double xmin  = signal->GetBinCenter(signal->GetMaximumBin()) + 20.;
    double xmax  = signal->GetBinCenter(signal->GetNbinsX());

    TF1* fun_BL = new TF1("fun_BL", "pol0", 0, 50);
    signal->Fit(fun_BL, opt);

    TF1* fun_dec = new TF1("fun_dec", "[0]*exp(-(x-[1])/[2])", xmin, xmin + 100);
    fun_dec->SetParameters(100., 100., 10.);
    fun_dec->SetParNames("A", "t0", "tau");
    signal->Fit(fun_dec, opt);

    TF1* fun_all = new TF1("fall", funDecaySingle, xmin, xmax, 4);

    fun_all->SetParNames("A", "t0", "tau", "const");
    fun_all->SetParameter(0, fun_dec->GetParameter(0));
    fun_all->SetParameter(1, fun_dec->GetParameter(1));
    fun_all->SetParameter(2, fun_dec->GetParameter(2));
    fun_all->FixParameter(3, fun_BL->GetParameter(0));
    int fitStat = signal->Fit(fun_all, opt);

    if (fitStat != 0)
    {
        std::cerr << "##### Warning in SFTimeConst::FitDecayTimeSingle()" << std::endl;
        std::cerr << "\t fit status: " << fitStat << std::endl;
    }

    if (ch == 0)
    {
        fFitResultsCh0[index]->SetFromFunction(fun_all);
        if (fitStat != 0) fFitResultsCh0[index]->SetStat(-1);
        fFitResultsCh0[index]->Print();
    }
    else if (ch == 1)
    {
        fFitResultsCh1[index]->SetFromFunction(fun_all);
        if (fitStat != 0) fFitResultsCh1[index]->SetStat(-1);
        fFitResultsCh1[index]->Print();
    }

    return true;
}
//------------------------------------------------------------------
/// This function performs fitting to the given TProfile signal.
/// Double decay function if fitted to the falling slope of the signal.
/// Results of the fit are subsequently written in the SFFitResults
/// class object. Function returns true if fitting was successful and
/// fit results are valid.
/// \param signal - averaged signal (TProfile)
/// \param ID - measurement ID
bool SFTimeConst::FitDecayTimeDouble(TProfile* signal, int ID)
{

    TString opt;
    if (fVerb)
        opt = "R0";
    else
        opt = "QR0";

    int     ch    = -1;
    TString hname = signal->GetName();
    if (hname.Contains("ch0"))
        ch = 0;
    else if (hname.Contains("ch1"))
        ch = 1;
    else
    {
        std::cerr << "##### Error in SFTimeConst::FitDecayTimeDouble()!" << std::endl;
        std::cerr << "Could not interpret signal name!" << std::endl;
        return false;
    }

    std::vector<int> measurementsIDs = fData->GetMeasurementsIDs();
    int              index           = SFTools::GetIndex(measurementsIDs, ID);
    double           xmin            = signal->GetBinCenter(signal->GetMaximumBin()) + 20;
    double           xmax            = signal->GetBinCenter(signal->GetNbinsX());

    TF1* fun_BL = new TF1("fun_BL", "pol0", 0, 50);
    signal->Fit(fun_BL, opt);

    TF1* fun_fast = new TF1("fun_fast", "[0]*exp(-(x-[1])/[2])", xmin, xmin + 130);
    fun_fast->SetParameters(100., 100., 10.);
    fun_fast->SetParNames("A", "t0", "tau");
    signal->Fit(fun_fast, opt);

    TF1* fun_slow = new TF1("fun_slow", "[0]*exp(-(x-[1])/[2])", 400, xmax);
    fun_slow->SetParameters(100., 100., 400.);
    fun_slow->SetParNames("A", "t0", "tau");
    fun_slow->FixParameter(1, fun_fast->GetParameter(1));
    signal->Fit(fun_slow, opt);

    TF1* fun_all = new TF1("fall", funDecayDouble, xmin, xmax, 6);

    fun_all->SetParNames("A_fast", "t0", "tau_fast", "A_slow", "tau_slow", "const");
    fun_all->SetParameter(0, fun_fast->GetParameter(0));
    fun_all->FixParameter(1, fun_fast->GetParameter(1));
    fun_all->SetParameter(2, fun_fast->GetParameter(2));
    fun_all->SetParameter(3, fun_slow->GetParameter(0));
    fun_all->SetParameter(4, fun_slow->GetParameter(2));
    fun_all->FixParameter(5, fun_BL->GetParameter(0));
    fun_all->FixParameter(1, xmin - 20);
    int fitStat = signal->Fit(fun_all, opt);

    if (fitStat != 0)
    {
        std::cerr << "##### Warning in SFTimeConst::FitDecayTimeDouble()" << std::endl;
        std::cerr << "\t fit status: " << fitStat << std::endl;
    }

    if (ch == 0)
    {
        fFitResultsCh0[index]->SetFromFunction(fun_all);
        if (fitStat != 0) fFitResultsCh0[index]->SetStat(-1);
        fFitResultsCh0[index]->Print();
    }
    else if (ch == 1)
    {
        fFitResultsCh1[index]->SetFromFunction(fun_all);
        if (fitStat != 0) fFitResultsCh1[index]->SetStat(-1);
        fFitResultsCh1[index]->Print();
    }

    return true;
}
//------------------------------------------------------------------
/// Fits all signals of the analyzed series from both channels.
/// This function also calculates average time constants with their
/// uncertainties.
bool SFTimeConst::FitAllSignals(void)
{

    int              n               = fData->GetNpoints();
    std::vector<int> measurementsIDs = fData->GetMeasurementsIDs();
    TString          fiber           = fData->GetFiber();

    if (fiber.Contains("LuAG") || fiber.Contains("GAGG"))
    {
        for (int i = 0; i < n; i++)
        {
            FitDecayTimeDouble(fSignalsCh0[i], measurementsIDs[i]);
            FitDecayTimeDouble(fSignalsCh1[i], measurementsIDs[i]);
        }
    }
    else if (fiber.Contains("LYSO"))
    {
        for (int i = 0; i < n; i++)
        {
            FitDecayTimeSingle(fSignalsCh0[i], measurementsIDs[i]);
            FitDecayTimeSingle(fSignalsCh1[i], measurementsIDs[i]);
        }
    }
    else
    {
        std::cerr << "##### Error in SFTimeConst::FitAllSignals()!" << std::endl;
        std::cerr << "Unknown fiber material!" << std::endl;
        return false;
    }

    //----- Calculating average time constants
    //----- and intensities
    double statCh0, statCh1;
    int    counter = 0;

    //----- double decay
    double fastDec, fastDecErr;
    double slowDec, slowDecErr;
    double fastDecSum    = 0;
    double slowDecSum    = 0;
    double fastDecSumErr = 0;
    double slowDecSumErr = 0;

    double fastAmp, fastAmpErr;
    double slowAmp, slowAmpErr;
    double fastAmpSum    = 0;
    double slowAmpSum    = 0;
    double fastAmpSumErr = 0;
    double slowAmpSumErr = 0;

    if (fiber.Contains("LuAG") || fiber.Contains("GAGG"))
    {

        for (int i = 0; i < n; i++)
        {
            statCh0 = fFitResultsCh0[i]->GetStat();
            statCh1 = fFitResultsCh1[i]->GetStat();

            if (statCh0 == 0)
            {
                fFitResultsCh0[i]->GetFastDecTime(fastDec, fastDecErr);
                fFitResultsCh0[i]->GetSlowDecTime(slowDec, slowDecErr);
                fastDecSum += fastDec * (1. / pow(fastDecErr, 2));
                slowDecSum += slowDec * (1. / pow(slowDecErr, 2));
                fastDecSumErr += 1. / pow(fastDecErr, 2);
                slowDecSumErr += 1. / pow(slowDecErr, 2);

                fFitResultsCh0[i]->GetAmpFast(fastAmp, fastAmpErr);
                fFitResultsCh0[i]->GetAmpSlow(slowAmp, slowAmpErr);
                fastAmpSum += fastAmp * (1. / pow(fastAmpErr, 2));
                slowAmpSum += slowAmp * (1. / pow(slowAmpErr, 2));
                fastAmpSumErr += 1. / pow(fastAmpErr, 2);
                slowAmpSumErr += 1. / pow(slowAmpErr, 2);

                counter++;
            }

            if (statCh1 == 0)
            {
                fFitResultsCh1[i]->GetFastDecTime(fastDec, fastDecErr);
                fFitResultsCh1[i]->GetSlowDecTime(slowDec, slowDecErr);
                fastDecSum += fastDec * (1. / pow(fastDecErr, 2));
                slowDecSum += slowDec * (1. / pow(slowDecErr, 2));
                fastDecSumErr += 1. / pow(fastDecErr, 2);
                slowDecSumErr += 1. / pow(slowDecErr, 2);

                fFitResultsCh1[i]->GetAmpFast(fastAmp, fastAmpErr);
                fFitResultsCh1[i]->GetAmpSlow(slowAmp, slowAmpErr);
                fastAmpSum += fastAmp * (1. / pow(fastAmpErr, 2));
                slowAmpSum += slowAmp * (1. / pow(slowAmpErr, 2));
                fastAmpSumErr += 1. / pow(fastAmpErr, 2);
                slowAmpSumErr += 1. / pow(slowAmpErr, 2);

                counter++;
            }
        }

        double fastDecAv    = fastDecSum / fastDecSumErr;
        double slowDecAv    = slowDecSum / slowDecSumErr;
        double fastDecAvErr = sqrt(1. / fastDecSumErr);
        double slowDecAvErr = sqrt(1. / slowDecSumErr);

        double fastAmpAve    = fastAmpSum / fastAmpSumErr;
        double slowAmpAve    = slowAmpSum / slowAmpSumErr;
        double fastAmpAveErr = sqrt(1. / fastAmpSumErr);
        double slowAmpAveErr = sqrt(1. / slowAmpSumErr);

        double denom   = fastAmpAve * fastDecAv + slowAmpAve * slowDecAv;
        double iFastAv = ((fastAmpAve * fastDecAv) / denom) * 100;
        double iSlowAv = ((slowAmpAve * slowDecAv) / denom) * 100;

        std::cout << "\n\n----------------------------------" << std::endl;
        std::cout << "Average decay constants for whole series:" << std::endl;
        std::cout << "Fast decay: " << fastDecAv << " +/- "
                  << fastDecAvErr << " ns" << std::endl;
        std::cout << "Fast component intensity: " << iFastAv << " %" << std::endl;
        std::cout << "Slow decay: " << slowDecAv << " +/- " 
                  << slowDecAvErr << " ns" << std::endl;
        std::cout << "Slow component intensity: " << iSlowAv << " %" << std::endl;
        std::cout << "Counter: " << counter << std::endl;
        std::cout << "----------------------------------" << std::endl;
        
        fResults->AddResult(SFResultTypeNum::kFastDecay, fastDecAv, fastDecAvErr);
        fResults->AddResult(SFResultTypeNum::kIFast, iFastAv, -1);
        fResults->AddResult(SFResultTypeNum::kSlowDecay, slowDecAv, slowDecAvErr);
        fResults->AddResult(SFResultTypeNum::kISlow, iSlowAv, -1);
    }

    //----- For single decay time
    double dec, decErr;
    double decSum    = 0;
    double decSumErr = 0;

    if (fiber.Contains("LYSO"))
    {

        for (int i = 0; i < n; i++)
        {
            statCh0 = fFitResultsCh0[i]->GetStat();
            statCh1 = fFitResultsCh1[i]->GetStat();

            if (statCh0 == 0)
            {
                fFitResultsCh0[i]->GetDecTime(dec, decErr);
                decSum += dec * (1. / pow(decErr, 2));
                decSumErr += 1. / pow(decErr, 2);
                counter++;
            }

            if (statCh1 == 0)
            {
                fFitResultsCh1[i]->GetDecTime(dec, decErr);
                decSum += dec * (1. / pow(decErr, 2));
                decSumErr += 1. / pow(decErr, 2);
                counter++;
            }
        }

        double fastDecAv    = decSum / decSumErr;
        double fastDecAvErr = sqrt(1. / decSumErr);
        double slowDecAv    = 0;
        double slowDecAvErr = 0;
        double iFastAv      = 100;
        double iSlowAv      = 0;

        std::cout << "\n\n----------------------------------" << std::endl;
        std::cout << "Average decay constant for the whole series:" << std::endl;
        std::cout << "Decay constant: " << fastDecAv << " +/- " 
                  << fastDecAvErr << " ns" << std::endl;
        std::cout << "Counter: " << counter << std::endl;
        std::cout << "----------------------------------" << std::endl;
        
        fResults->AddResult(SFResultTypeNum::kFastDecay, fastDecAv, fastDecAvErr);
        fResults->AddResult(SFResultTypeNum::kIFast, iFastAv, -1);
        fResults->AddResult(SFResultTypeNum::kSlowDecay, slowDecAv, slowDecAvErr);
        fResults->AddResult(SFResultTypeNum::kISlow, iSlowAv, -1);
    }

    return true;
}
//------------------------------------------------------------------
/// Fits all signals of the analyzed series from the chosen channel.
/// \param ch - channel number.
bool SFTimeConst::FitAllSignals(int ch)
{

    int              n               = fData->GetNpoints();
    std::vector<int> measurementsIDs = fData->GetMeasurementsIDs();
    TString          fiber           = fData->GetFiber();

    if (ch != 0 || ch != 1)
    {
        std::cerr << "##### Error in SFTimeConst::FitAllSignals()!" << std::endl;
        std::cerr << "Incorrect channel number. Possible options: 0 or 1" << std::endl;
        return false;
    }

    if (fiber.Contains("LuAG") || fiber.Contains("GAGG"))
    {
        for (int i = 0; i < n; i++)
        {
            if (ch == 0)
                FitDecayTimeDouble(fSignalsCh0[i], measurementsIDs[i]);
            else if (ch == 1)
                FitDecayTimeDouble(fSignalsCh1[i], measurementsIDs[i]);
        }
    }
    else if (fiber.Contains("LYSO"))
    {
        for (int i = 0; i < n; i++)
        {
            if (ch == 0)
                FitDecayTimeSingle(fSignalsCh0[i], measurementsIDs[i]);
            else if (ch == 1)
                FitDecayTimeSingle(fSignalsCh1[i], measurementsIDs[i]);
        }
    }
    else
    {
        std::cerr << "##### Error in SFTimeConst::FitAllSignal()!" << std::endl;
        std::cerr << "Unknown fiber material!" << std::endl;
        return false;
    }

    return true;
}
//------------------------------------------------------------------
/// Returns vector containing signals for whole series and
/// chosen channel
///\param ch - channel number
std::vector<TProfile*> SFTimeConst::GetSignals(int ch)
{

    if ((ch == 0 && fSignalsCh0.empty()) || (ch == 1 && fSignalsCh1.empty()))
    {
        std::cerr << "##### Error in SFTimeConst::GetSignals()!" << std::endl;
        std::cerr << "No signals available!" << std::endl;
        std::abort();
    }

    std::vector<TProfile*> signals;
    
    if (ch == 0)
        signals = fSignalsCh0;
    else if (ch == 1)
        signals = fSignalsCh1;
    else
    {
        std::cerr << "##### Error in SFTimeConst::GetSignals()!" << std::endl;
        std::cerr << "No signals available for channel: " << ch << std::endl;
        std::abort();
    }
    
    return signals;
}
//------------------------------------------------------------------
/// Returns fitting results for a requested channel.
/// \param ch - channel number
std::vector <SFFitResults*> SFTimeConst::GetFitResults(int ch)
{
    if ((ch == 0 && fFitResultsCh0.empty()) ||
        (ch == 1 && fFitResultsCh1.empty()))
    {
        std::cerr << "##### Error in SFTimeConst::GetFitResults()!" << std::endl;
        std::cerr << "No results available!" << std::endl;
        std::abort();
    }
    
    std::vector<SFFitResults*> results;
    
    if (ch == 0)
        results = fFitResultsCh0;
    else if (ch == 1)
        results = fFitResultsCh1;
    else
    {
        std::cerr << "##### Error in SFTimeConst::GetFitResults()!" << std::endl;
        std::cerr << "No results available for channel: " << ch << std::endl;
        std::abort();
    }
    
    return results;
}
//------------------------------------------------------------------
/// Prints details of the SFTimeConst class object.
void SFTimeConst::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFTimeConst class object" << std::endl;
    std::cout << "Exparimental series number: " << fSeriesNo << std::endl;
    std::cout << "Signals PE: " << fPE << std::endl;
    std::cout << "Verbose level: " << fVerb << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
