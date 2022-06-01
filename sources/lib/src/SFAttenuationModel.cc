// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *         SFAttenuationModel.cc         *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFAttenuationModel.hh"

int iparSl[] = {0, 1, 2, 3, 4, 5};
int iparSr[] = {0, 1, 2, 3, 4, 5};

ClassImp(SFAttenuationModel);

//------------------------------------------------------------------
/// Standard constructor.
/// \param seriesNo - number of the experimental series.
SFAttenuationModel::SFAttenuationModel(int seriesNo) : fSeriesNo(seriesNo), 
                                                       fData(nullptr),
                                                       fMAttCh0Graph(nullptr),
                                                       fMAttCh1Graph(nullptr),
                                                       fMAttCh0CorrGraph(nullptr),
                                                       fMAttCh1CorrGraph(nullptr),
                                                       fResults(nullptr),
                                                       fFitterResults(nullptr)
{
    try
    {
        fData = new SFData(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFAttenuationModel constructor!";
    }

    TString desc = fData->GetDescription();
    if (!desc.Contains("Regular series"))
    {
        std::cout << "##### Error in SFAttenuationModel constructor! Non-regular series!"
                  << std::endl;
        throw "##### Exception in SFAttenuationModel constructor!";
    }

    SFAttenuation* att;
    try
    {
        att = new SFAttenuation(fSeriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        throw "##### Exception in SFAttenuationModel constructor!";
    }
    
    att->AttSeparateCh(0);
    att->AttSeparateCh(1);

    std::vector<SFResults*> att_res = att->GetResults();
    fMAttCh0Graph = (TGraphErrors*)att_res[0]->GetObject(SFResultTypeObj::kAttGraph);
    fMAttCh1Graph = (TGraphErrors*)att_res[1]->GetObject(SFResultTypeObj::kAttGraph);

    fResults = new SFResults(Form("ReconstructionResults_S%i_Mod", fSeriesNo));
}
//------------------------------------------------------------------
/// Destructor.
SFAttenuationModel::~SFAttenuationModel()
{
    if (fData != nullptr) delete fData;
};
//------------------------------------------------------------------
/// Calculates partial derivatives of the reconstructed primary components
/// for requested side: L (left) or R (right). Calculated values are returned
/// in an array. Order in the array is the following: [0] - lambda, [1] - 
/// etha right, [2] - etha left, [3] - ksi, [4] - L, [5] - S left, [6] -
/// S right.
double* ConstructDerivativesMatrix(std::vector<double> par, TString side)
{
    /*-----
     * params:
     * 0 - lambda
     * 1 - eta right
     * 2 - eta left
     * 3 - ksi
     * 4 - L
     * 5 - S left
     * 6 - S right
    -----*/

    /*
    std::cout << "lambda: " << par[0] <<std::endl;
    std::cout << "eta right: " << par[1] <<std::endl;
    std::cout << "eta left: " << par[2] <<std::endl;
    std::cout << "ksi: " << par[3] <<std::endl;
    std::cout << "L: " << par[4] <<std::endl;
    std::cout << "S left: " << par[5] <<std::endl;
    std::cout << "S right: " << par[6] <<std::endl;
    */

    double dPdLambda, dPdEtaR, dPdEtaL, dPdKsi;

    double e = exp(par[4] / par[0]);

    if (side == "L")
    {
        dPdLambda = -(e * par[4] * par[1] *
                     (-2 * e * par[3] * par[5] * par[2] + par[6] * (pow(e, 2) + par[1] * par[2]))) /
                     (pow(par[0], 2) * par[3] * pow((pow(e, 2) - par[1] * par[2]), 2));

        dPdEtaR = (pow(e, 2) * (-e * par[6] + par[3] * par[5] * par[2])) /
                  (par[3] * (pow((pow(e, 2) - par[1] * par[2]), 2)));

        dPdEtaL = (e * par[1] * (e * par[3] * par[5] - par[6] * par[1])) /
                  (par[3] * pow((pow(e, 2) - par[1] * par[2]), 2));

        dPdKsi = (e * par[6] * par[1]) / (pow(par[3], 2) * (pow(e, 2) - par[1] * par[2]));
    }
    else if (side == "R")
    {
        dPdLambda = -(e * par[4] * par[2] *
                     (-2 * e * par[6] * par[1] + par[3] * par[5] * (pow(e, 2) + par[1] * par[2]))) /
                     (pow(par[0], 2) * par[3] * (pow((pow(e, 2) - par[1] * par[2]), 2)));

        dPdEtaR = (e * par[2] * (e * par[6] - par[3] * par[5] * par[2])) /
                  (par[3] * (pow((pow(e, 2) - par[1] * par[2]), 2)));

        dPdEtaL = (-pow(e, 3) * par[3] * par[5] + pow(e, 2) * par[6] * par[1]) /
                  (par[3] * (pow((pow(e, 2) - par[1] * par[2]), 2)));

        dPdKsi = -(pow(e, 2) * par[6]) / (pow(par[3], 2) * (pow(e, 2) - par[1] * par[2]));
    }

    double* matrix = new double[6];
    matrix[0]      = 0;
    matrix[1]      = dPdLambda;
    matrix[2]      = dPdEtaR;
    matrix[3]      = dPdEtaL;
    matrix[4]      = dPdKsi;
    matrix[5]      = 0;

    return matrix;
}
//------------------------------------------------------------------
/// Multiplies an array containing partial derivatives of the 
/// reconstructed primary components with the covariance matrix. 
/// Returns final value of the multiplication. 
double Multiply(double* deriv, TMatrixD cov)
{
    // cov - covariance matrix: 6 rows & 6 columns
    // deriv - vector of derivatives: 1 row & 6 columns

    int    dim = 6;
    double var = 0;

    //----- multiply cov x deriv
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            var += cov[i][j] * deriv[i] * deriv[j];
        }
    }

    if (var < 0)
    {
        std::cerr << "##### Error in Multiply: negative variance!" << std::endl;
        std::cerr << "var = " << var << std::endl;
        std::abort();
    }

    return var;
}
//------------------------------------------------------------------
/// Calculates uncertainty of the reconstructed primary component.
double SFAttenuationModel::CalculateUncertainty(std::vector<double> params, TString side)
{
    double uncert = 0;
    
    std::vector<double> parsForDerivatives(7);
    parsForDerivatives[0] = params[0]; //fResults->GetValue(SFResultTypeNum::kLambda);
    parsForDerivatives[1] = params[1]; //fResults->GetValue(SFResultTypeNum::kEtaR);
    parsForDerivatives[2] = params[2]; //fResults->GetValue(SFResultTypeNum::kEtaL);
    parsForDerivatives[3] = params[3]; //fResults->GetValue(SFResultTypeNum::kKsi);
    parsForDerivatives[4] = params[4]; //fResults->GetValue(SFResultTypeNum::kLength);
    parsForDerivatives[5] = params[5]; //SL
    parsForDerivatives[6] = params[6]; //SR
    // params[7] - sigma_SL
    // params[8] - sigma_SR 
    
    double* derivativesL;
    double* derivativesR;
    
    if (side == "L")
    {
        derivativesL = ConstructDerivativesMatrix(parsForDerivatives, "L");

        double sigma_fSL = Multiply(derivativesL, fCovMatrix);

        double dPLdSL = (exp(2 * params[4] / params[0])) /
                        (exp(2 * params[4] / params[0]) -
                        params[1] * params[2]);
        double dPLdSR = (exp(params[4] / params[0]) *
                        params[1]) / (params[3] *
                        (exp(2 * params[4] / params[0]) -
                        params[1] * params[2]));

        double signalLCh0Err = sqrt(sigma_fSL + pow(dPLdSL * params[7], 2) + pow(dPLdSR * params[8], 2));
        uncert = signalLCh0Err;
    }
    else if (side == "R")
    {
        derivativesR = ConstructDerivativesMatrix(parsForDerivatives, "R");
        
        double sigma_fSR = Multiply(derivativesR, fCovMatrix);

        double dPRdSL = (-exp(params[4] / params[0]) *
                        params[2]) / (exp(2 * params[4] /
                        params[0]) - params[1] *
                        params[2]);
        double dPRdSR = (exp(2 * params[4] / params[0])) /
                        (params[3] * (exp(2 * params[4] /
                        params[0]) - params[1] *
                        params[2]));
                        
        double signalRCh1Err = sqrt(sigma_fSR + pow(dPRdSL * params[7], 2) + pow(dPRdSR * params[8], 2));
        uncert = signalRCh1Err;
    }
    else 
    {
        std::cerr << "Error in SFAttenuationModel::CalculateUncertainty()" << std::endl;
        std::cerr << "Incorrect side. Possible optios are: L and R" << std::endl;
        std::abort();
    }
    
    //delete[] derivativesL;
    //delete[] derivativesR;
    
    return uncert;
}
//------------------------------------------------------------------
/// Fits exponential attenuation model with light reflection to the
/// experimental data. Model equations are fitted simultaneously to
/// data sets for left and right side of the fiber. Additionally,
/// primary component is reconstructed and presented in graphs.
bool SFAttenuationModel::FitModel(void)
{
    std::cout << "\n\n----- Inside SFAttenuationModel::FitModel() for series " << fSeriesNo << "\n"
              << std::endl;

    //----- setting formulas start -----
    double fiberLen = fData->GetFiberLength();

    TFormula* form_Pl = new TFormula("Pl", "[0]*exp(-x/[1])");
    TFormula* form_Pr = new TFormula("Pr", "[0]*exp(-([5]-x)/[1])");

    TF1* fun_Pl = new TF1("fun_Pl", "Pl(x)", 0, fiberLen);
    TF1* fun_Pr = new TF1("fun_Pr", "Pr(x)", 0, fiberLen);

    TFormula* form_Rl = new TFormula("Rl", "[2]*Pr(x)*exp(-[5]/[1])");
    TFormula* form_Rr = new TFormula("Rr", "[3]*Pl(x)*exp(-[5]/[1])");

    TF1* fun_Rl = new TF1("fun_Rl", "Rl(x)", 0, fiberLen);
    TF1* fun_Rr = new TF1("fun_Rr", "Rr(x)", 0, fiberLen);
    
    TFormula* form_Sl = new TFormula("Sl", "(Pl(x)+Rl(x))");
    TFormula* form_Sr = new TFormula("Sr", "[4]*(Pr(x)+Rr(x))");

    TF1* fun_Sl = new TF1("fun_Sl", "Sl(x)", 0, fiberLen);
    TF1* fun_Sr = new TF1("fun_Sr", "Sr(x)", 0, fiberLen);

    /*-----
     * [0] - S0
     * [1] - lambda
     * [2] - eta right
     * [3] - eta left
     * [4] - ksi
     * [5] - fiber length
    -----*/

    //----- setting formulas end -----

    //----- fitting start -----
    ROOT::Math::WrappedMultiTF1 wfSl(*fun_Sl, 1);
    ROOT::Math::WrappedMultiTF1 wfSr(*fun_Sr, 1);

    ROOT::Fit::DataOptions opt;

    ROOT::Fit::DataRange rangeSl;
    rangeSl.SetRange(0, fiberLen);
    ROOT::Fit::BinData dataSl(opt, rangeSl);
    ROOT::Fit::FillData(dataSl, fMAttCh0Graph); // ch0 -> L

    ROOT::Fit::DataRange rangeSr;
    rangeSr.SetRange(0, fiberLen);
    ROOT::Fit::BinData dataSr(opt, rangeSr);
    ROOT::Fit::FillData(dataSr, fMAttCh1Graph); // ch1 -> R

    ROOT::Fit::Chi2Function chi2_Sl(dataSl, wfSl);
    ROOT::Fit::Chi2Function chi2_Sr(dataSr, wfSr);

    GlobalChi2Model globalChi2(chi2_Sl, chi2_Sr);

    ROOT::Fit::Fitter fitter;

    const int npar       = 6;
//     double    par0[npar] = {500, 150, 0.5, 0.5, 1, fiberLen};
    double    par0[npar] = {100, 350, 0.5, 0.5, 0.5, fiberLen};

    fitter.Config().SetParamsSettings(npar, par0);

    fitter.Config().ParSettings(0).SetLimits(0, 1E6); //par limits for series 1-262
    fitter.Config().ParSettings(1).SetLimits(0, 1E4);
    fitter.Config().ParSettings(2).SetLimits(-1, 3);
    fitter.Config().ParSettings(3).SetLimits(-1, 3);
    fitter.Config().ParSettings(4).SetLimits(0, 5);
    
//     fitter.Config().ParSettings(0).SetLimits(0, 1E10); //par limits for series 263
//     fitter.Config().ParSettings(1).SetLimits(0, 1E4);
//     fitter.Config().ParSettings(2).SetLimits(-1, 10);
//     fitter.Config().ParSettings(3).SetLimits(-1, 10);
//     fitter.Config().ParSettings(4).SetLimits(0, 5);
    fitter.Config().ParSettings(5).Fix();

    fitter.Config().ParSettings(0).SetName("S0");
    fitter.Config().ParSettings(1).SetName("Latt");
    fitter.Config().ParSettings(2).SetName("EtaR");
    fitter.Config().ParSettings(3).SetName("EtaL");
    fitter.Config().ParSettings(4).SetName("Ksi");
    fitter.Config().ParSettings(5).SetName("L");

    fitter.Config().MinimizerOptions().SetPrintLevel(1);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");
    fitter.Config().MinimizerOptions().SetMaxIterations(1e6);
    fitter.Config().MinimizerOptions().SetMaxFunctionCalls(1e6);

    fitter.FitFCN(npar, globalChi2, 0, dataSl.Size() + dataSr.Size(), true);
    fFitterResults = new ROOT::Fit::FitResult(fitter.Result());
    //----- fitting end -----

    //----- setting numerical results start -----
    const double* params = fFitterResults->GetParams();
    const double* errors = fFitterResults->GetErrors();

    fResults->AddResult(SFResultTypeNum::kS0, params[0], errors[0]);
    fResults->AddResult(SFResultTypeNum::kLambda, params[1], errors[1]);
    fResults->AddResult(SFResultTypeNum::kEtaR, params[2], errors[2]);
    fResults->AddResult(SFResultTypeNum::kEtaL, params[3], errors[3]);
    fResults->AddResult(SFResultTypeNum::kKsi, params[4], errors[4]);
    fResults->AddResult(SFResultTypeNum::kLength, params[5], 0);
    
    fResults->AddResult(SFResultTypeNum::kChi2NDF, 
                        fFitterResults->Chi2() / fFitterResults->Ndf(), -1);

    fun_Pl->SetFitResult(*fFitterResults, iparSl);
    fun_Pl->SetRange(rangeSl().first, rangeSl().second);

    fun_Pr->SetFitResult(*fFitterResults, iparSr);
    fun_Pr->SetRange(rangeSr().first, rangeSr().second);

    fun_Rl->SetFitResult(*fFitterResults, iparSl);
    fun_Rl->SetRange(rangeSl().first, rangeSl().second);

    fun_Rr->SetFitResult(*fFitterResults, iparSr);
    fun_Rr->SetRange(rangeSr().first, rangeSr().second);

    fun_Sl->SetFitResult(*fFitterResults, iparSl);
    fun_Sl->SetRange(rangeSl().first, rangeSl().second);

    fun_Sr->SetFitResult(*fFitterResults, iparSr);
    fun_Sr->SetRange(rangeSr().first, rangeSr().second);
    //----- setting numerical results end -----

    //----- recalculating primary signal start -----
    const int           npoints    = fData->GetNpoints();
    std::vector<double> positions  = fData->GetPositions();
    TString             collimator = fData->GetCollimator();
    TString             testBench  = fData->GetTestBench();

    fMAttCh0CorrGraph = new TGraphErrors(npoints);
    fMAttCh0CorrGraph->SetName("fMAttCh0CorrGraph");
    fMAttCh0CorrGraph->SetTitle("Reconstructed emission Ch0");
    fMAttCh0CorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fMAttCh0CorrGraph->GetYaxis()->SetTitle("reconstructed signal [a.u.]");
    fMAttCh0CorrGraph->SetMarkerStyle(8);

    fMAttCh1CorrGraph = new TGraphErrors(npoints);
    fMAttCh1CorrGraph->SetName("fMAttCh1CorrGraph");
    fMAttCh1CorrGraph->SetTitle("Reconstructed emission Ch1");
    fMAttCh1CorrGraph->GetXaxis()->SetTitle("source position [mm]");
    fMAttCh1CorrGraph->GetYaxis()->SetTitle("reconstructed signal [a.u.]");
    fMAttCh1CorrGraph->SetMarkerStyle(8);

    // x - Sr - 1
    // y - Sl - 0

    TFormula* form_PlReco = new TFormula("PlReco", "(exp([5]/[1])*(exp([5]/[1])*[4]*y - x*[2]))/([4]*(exp(2*[5]/[1]) - [3]*[2])) + [0]*0");

    TF2* fun_PlReco = new TF2("fun_PlReco", "PlReco(x,y)", 0, 1000, 0, 1000);
    fun_PlReco->FixParameter(0, 0);
    fun_PlReco->FixParameter(1, fResults->GetValue(SFResultTypeNum::kLambda));
    fun_PlReco->FixParameter(2, fResults->GetValue(SFResultTypeNum::kEtaR));
    fun_PlReco->FixParameter(3, fResults->GetValue(SFResultTypeNum::kEtaL));
    fun_PlReco->FixParameter(4, fResults->GetValue(SFResultTypeNum::kKsi));
    fun_PlReco->FixParameter(5, fResults->GetValue(SFResultTypeNum::kLength));

    TFormula* form_PrReco = new TFormula("PrReco", "-(exp([5]/[1])*(-exp([5]/[1])*x + [4]*y*[3]))/([4]*(exp(2*[5]/[1]) - [3]*[2])) + [0]*0");

    TF2* fun_PrReco = new TF2("fun_PrReco", "PrReco(x,y)", 0, 1000, 0, 1000);

    fun_PrReco->FixParameter(0, 0);
    fun_PrReco->FixParameter(1, fResults->GetValue(SFResultTypeNum::kLambda));
    fun_PrReco->FixParameter(2, fResults->GetValue(SFResultTypeNum::kEtaR));
    fun_PrReco->FixParameter(3, fResults->GetValue(SFResultTypeNum::kEtaL));
    fun_PrReco->FixParameter(4, fResults->GetValue(SFResultTypeNum::kKsi));
    fun_PrReco->FixParameter(5, fResults->GetValue(SFResultTypeNum::kLength));

    fPlRecoFun = fun_PlReco;
    fPrRecoFun = fun_PrReco;

//  TMatrixD covMatrix(6, 6);
    fCovMatrix.ResizeTo(6, 6);
    fFitterResults->GetCovarianceMatrix(fCovMatrix);

//     double* derivativesL;
//     double* derivativesR;

    std::vector<double> parsForErrors(9);
    parsForErrors[0] = fResults->GetValue(SFResultTypeNum::kLambda);
    parsForErrors[1] = fResults->GetValue(SFResultTypeNum::kEtaR);
    parsForErrors[2] = fResults->GetValue(SFResultTypeNum::kEtaL);
    parsForErrors[3] = fResults->GetValue(SFResultTypeNum::kKsi);
    parsForErrors[4] = fResults->GetValue(SFResultTypeNum::kLength);
    
    //std::vector<double> parsForDerivatives(7);
    //parsForDerivatives[0] = fResults->GetValue(SFResultTypeNum::kLambda);
    //parsForDerivatives[1] = fResults->GetValue(SFResultTypeNum::kEtaR);
    //parsForDerivatives[2] = fResults->GetValue(SFResultTypeNum::kEtaL);
    //parsForDerivatives[3] = fResults->GetValue(SFResultTypeNum::kKsi);
    //parsForDerivatives[4] = fResults->GetValue(SFResultTypeNum::kLength);

    for (int i = 0; i < npoints; i++)
    {
        double SR = fMAttCh1Graph->GetPointY(i);
        double SL = fMAttCh0Graph->GetPointY(i);

        double sigmaSL = fMAttCh0Graph->GetErrorY(i);
        double sigmaSR = fMAttCh1Graph->GetErrorY(i);

        parsForErrors[5] = SL;
        parsForErrors[6] = SR;
        parsForErrors[7] = sigmaSL;
        parsForErrors[8] = sigmaSR;
        
        double signalLCh0 = fun_PlReco->Eval(SR, SL);
        double signalLCh0Err = CalculateUncertainty(parsForErrors, "L");

        //parsForDerivatives[5] = SL;
        //parsForDerivatives[6] = SR;

        //derivativesL = ConstructDerivativesMatrix(parsForDerivatives, "L");
        //derivativesR = ConstructDerivativesMatrix(parsForDerivatives, "R");

        //double sigma_fSL = Multiply(derivativesL, covMatrix);

        //double dPLdSL = (exp(2 * fun_PlReco->GetParameter(5) / fun_PlReco->GetParameter(1))) /
        //                (exp(2 * fun_PlReco->GetParameter(5) / fun_PlReco->GetParameter(1)) -
        //                fun_PlReco->GetParameter(2) * fun_PlReco->GetParameter(3));
        //double dPLdSR = (exp(fun_PlReco->GetParameter(5) / fun_PlReco->GetParameter(1)) *
        //                fun_PlReco->GetParameter(2)) / (fun_PlReco->GetParameter(4) *
        //                (exp(2 * fun_PlReco->GetParameter(5) / fun_PlReco->GetParameter(1)) -
        //                fun_PlReco->GetParameter(2) * fun_PlReco->GetParameter(3)));

        //double signalLCh0Err = sqrt(sigma_fSL + pow(dPLdSL * sigmaSL, 2) + pow(dPLdSR * sigmaSR, 2));

        fMAttCh0CorrGraph->SetPoint(i, positions[i], signalLCh0);
        fMAttCh0CorrGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), signalLCh0Err);

        double signalRCh1 = fun_PrReco->Eval(SR, SL);
        double signalRCh1Err = CalculateUncertainty(parsForErrors, "R");

        //double sigma_fSR = Multiply(derivativesR, covMatrix);

        //double dPRdSL = (-exp(fun_PrReco->GetParameter(5) / fun_PrReco->GetParameter(1)) *
        //                fun_PrReco->GetParameter(3)) / (exp(2 * fun_PrReco->GetParameter(5) /
        //                fun_PrReco->GetParameter(1)) * fun_PrReco->GetParameter(2) *
        //                fun_PrReco->GetParameter(3));

        //double dPRdSR = (exp(2 * fun_PrReco->GetParameter(5) / fun_PrReco->GetParameter(1))) /
        //                (fun_PrReco->GetParameter(4) * (exp(2 * fun_PrReco->GetParameter(5) /
        //                fun_PrReco->GetParameter(1)) - fun_PrReco->GetParameter(2) *
        //                fun_PrReco->GetParameter(3)));

        //double signalRCh1Err = sqrt(sigma_fSR + pow(dPRdSL * sigmaSL, 2) + pow(dPRdSR * sigmaSR, 2));

        fMAttCh1CorrGraph->SetPoint(i, positions[i], signalRCh1);
        fMAttCh1CorrGraph->SetPointError(i, SFTools::GetPosError(collimator, testBench), signalRCh1Err);
    }
    //----- recalculating primary signal end -----

    //----- setting objects start
    fResults->AddObject(SFResultTypeObj::kPlFun, fun_Pl);
    fResults->AddObject(SFResultTypeObj::kPrFun, fun_Pr);
    fResults->AddObject(SFResultTypeObj::kRlFun, fun_Rl);
    fResults->AddObject(SFResultTypeObj::kRrFun, fun_Rr);
    fResults->AddObject(SFResultTypeObj::kSlFun, fun_Sl);
    fResults->AddObject(SFResultTypeObj::kSrFun, fun_Sr);

    fResults->AddObject(SFResultTypeObj::kSlVsPosGraph, fMAttCh0Graph);
    fResults->AddObject(SFResultTypeObj::kSrVsPosGraph, fMAttCh1Graph);
    fResults->AddObject(SFResultTypeObj::kPlVsPosGraph, fMAttCh0CorrGraph);
    fResults->AddObject(SFResultTypeObj::kPrVsPosGraph, fMAttCh1CorrGraph);

    fResults->AddObject(SFResultTypeObj::kPlRecoFun, fun_PlReco);
    fResults->AddObject(SFResultTypeObj::kPrRecoFun, fun_PrReco);
    //----- setting objects end

//     delete[] derivativesL;
//     delete[] derivativesR;

    return true;
}
//------------------------------------------------------------------
/// Prints details of the SFAttenuationModel class object.
void SFAttenuationModel::Print(void)
{
    std::cout << "\n-------------------------------------------" << std::endl;
    std::cout << "This is print out of SFAttenuationModel class object" << std::endl;
    std::cout << "Experimental series number " << fSeriesNo << std::endl;
    std::cout << "-------------------------------------------\n" << std::endl;
}
//------------------------------------------------------------------
