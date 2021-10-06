// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFEnergyRes.cc              *
// *            Jonas Kasper               *
// *      kasper@physik.rwth-aachen.de     *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyRes.hh"
    
//------------------------------------------------------------------    
SFResults* SFEnergyRes::CalculateEnergyResSingle(TH1D* h)
{
 
    auto peakFin = std::make_unique<SFPeakFinder>(h, false);
    peakFin->FindPeakFit();
    auto fitResults = std::unique_ptr<SFResults>(peakFin->GetResults());
    
    double enRes = fitResults->GetValue(SFResultTypeNum::kPeakSigma) /
                   fitResults->GetValue(SFResultTypeNum::kPeakPosition);
    double enResErr = enRes *
                      sqrt(pow(fitResults->GetUncertainty(SFResultTypeNum::kPeakPosition), 2) /
                      pow(fitResults->GetValue(SFResultTypeNum::kPeakPosition), 2) +
                      pow(fitResults->GetUncertainty(SFResultTypeNum::kPeakSigma), 2) /
                      pow(fitResults->GetValue(SFResultTypeNum::kPeakSigma), 2));
    
    enRes    = enRes * 100;
    enResErr = enResErr * 100;
    
    TString r_name = h->GetName() + TString("_results");
    SFResults *results = new SFResults(r_name);
    results->AddResult(kEnergyRes, enRes, enResErr);
    
    return results;
}
//------------------------------------------------------------------
SFResults* SFEnergyRes::CalculateEnergyResSeries(TString suffix, double pos_uncert,
                                                 std::vector<double> positions,
                                                 std::vector<TH1D*> spectra)
{

    int npoints = spectra.size();

    TString       gname = Form("energy_resolution_%s", suffix.Data());
    TGraphErrors* graph = new TGraphErrors(npoints);
    graph->GetXaxis()->SetTitle("source position [mm]");
    graph->GetYaxis()->SetTitle("energy resolution [%]");
    graph->SetTitle(Form("Energy Resolution (%s)", suffix.Data()));
    graph->SetName(gname);
    graph->SetMarkerStyle(4);

    double     enRes       = 0;
    double     enResErr    = 0;
    double     enResAve    = 0;
    double     enResAveErr = 0;

    for (int i = 0; i < npoints; i++)
    { 
        auto singleER = std::unique_ptr<SFResults>(CalculateEnergyResSingle(spectra[i]));
        enRes = singleER->GetValue(kEnergyRes);
        enResErr = singleER->GetUncertainty(kEnergyRes);
        graph->SetPoint(i, positions[i], enRes);
        graph->SetPointError(i, pos_uncert, enResErr);
        enResAve += enRes * (1. / pow(enResErr, 2));
        enResAveErr += (1. / pow(enResErr, 2));
    }

    enResAve    = enResAve / enResAveErr;
    enResAveErr = sqrt(1. / enResAveErr);

    std::cout << "Average energy resolution (" << suffix << "): " 
              << enResAve << " +/- " << enResAveErr << " % \n" << std::endl;

    SFResults *results = new SFResults(Form("AverageEnergyResolution_%s", suffix.Data()));
    
    results->AddResult(SFResultTypeNum::kEnergyRes, enResAve, enResAveErr);
    results->AddObject(SFResultTypeObj::kEnergyResGraph, graph);
    
    return results;
}
