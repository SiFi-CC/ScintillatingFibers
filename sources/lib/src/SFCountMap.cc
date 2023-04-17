// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             SFCountMap.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFCountMap.hh"

//------------------------------------------------------------------
/// Creates graph of counts registered in the fibers as a function of 
/// the fiber number for a chosen fiber side (channel 0/1). Number of 
/// counts is obtained as an integral charge spectra. IntegralAndError()
/// function of ROOT framework is used for that purpose. Average number
/// of counts is calulated as weighted mean with standard error of the 
/// weighted mean for each channel. Additionally, standard deviation of 
/// number of counts is calculated.
/// \param side is fiber side (l or r)
/// \param spectra vector containing charge spectra
SFResults* SFCountMap::DrawCounts(TString side, std::vector<TH1D*> spectra,
                                  std::vector<int> fiberNo)
{
    int npoints = spectra.size();
    
    if (spectra.empty() || spectra[0] == nullptr)
    {
        std::cerr << "##### Error in SFCountMap::" << __func__ << std::endl;
        std::cerr << "No spectra found! Please check!" << std::endl;
        std::abort();
    }

    TGraphErrors *graph = new TGraphErrors(npoints);
    graph->SetName(Form("CountsMap_%s", side.Data()));
    graph->SetTitle(Form("Counts per fiber %s", side.Data()));
    graph->GetXaxis()->SetTitle("fiber number");
    graph->GetYaxis()->SetTitle("counts");
    graph->SetMarkerStyle(4);
    
    double integral    = 0;
    double integralErr = 0;
    int    nbins       = 0;
    
    std::vector<double> integralsAll;
    
    double countsAv    = 0;
    double countsAvErr = 0;
    
    for (int i=0; i<npoints; i++)
    {
        nbins = spectra[i]->GetNbinsX();
        integral = spectra[i]->IntegralAndError(1, nbins, integralErr);
        
        graph->SetPoint(i, fiberNo[i], integral);
        graph->SetPointError(i, 0, integralErr);
        
        countsAv += integral * (1./pow(integralErr, 2));
        countsAvErr += 1./pow(integralErr, 2);
        integralsAll.push_back(integral);
    }
    
    countsAv    = countsAv / countsAvErr;
    countsAvErr = sqrt(1. / countsAvErr);
    
    double stdDev = SFTools::GetStandardDev(integralsAll);
    
    SFResults *results = new SFResults(Form("CountsMap_%s", side.Data()));
    results->AddObject(SFResultTypeObj::kCountsGraph, graph);
    results->AddResult(SFResultTypeNum::kCounts, countsAv, countsAvErr);
    results->AddResult(SFResultTypeNum::kCountsStdDev, stdDev, 0);
    
    std::cout << "\n\tAveraged counts for side " << side << ": " << countsAv 
              << " +/- " << countsAvErr << std::endl;
    std::cout << "\tStandard deviation: " << stdDev << std::endl << std::endl;
    
    return results;
}
