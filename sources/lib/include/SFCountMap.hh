// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             SFCountMap.hh             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************
#ifndef __SFCountMap_H_
#define __SFCounrMap_H_ 1

#include "SFData.hh"
#include "SFResults.hh"

#include <TObject.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>

#include <iostream>

/// This class is meant for analysis of measurements recorded with
/// the entire detector module (test bench: PMI, description: Module 
/// series). It allows to plot number of counts registered in the 
/// fiber as a function of the fiber number. Two fiber sides are plotted 
/// separately. In order to obtain number of counts the charge spectra 
/// are integrated with the ROOT's method IntegralAndError(). 

class SFCountMap : public TObject
{

private:
    int     fSeriesNo; ///< Number of experimental series to be analyzed
    SFData* fData;     ///< SFData object of the experimental series
    
    std::vector<TH1D*> fSpectraCh0; ///< Vector containing spectra (channel 0) 
    std::vector<TH1D*> fSpectraCh1; ///< Vector containing spectra (channel 1)
    
    SFResults* fResultsCh0; ///< Results for channel 0
    SFResults* fResultsCh1; ///< Results for channel 1
    
public:
    SFCountMap(int seriesNo);
    ~SFCountMap();
    
    bool                    DrawCounts(int ch);
    std::vector<TH1D*>      GetSpectra(int ch);
    std::vector<SFResults*> GetResults(void);
    void                    Print(void);
    
    ClassDef(SFCountMap, 1)
};

#endif /* __SFCountMap_H_ */
