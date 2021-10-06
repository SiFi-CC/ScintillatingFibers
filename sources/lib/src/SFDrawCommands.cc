// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           SFDrawCommands.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFDrawCommands.hh"

const double ampMax = 660; // maximum valid amplitude [mV] in measurements 
                           // with the Desktop Digitizer; signals with the amplitude
                           // above this value are saturated and need to be discarded

//------------------------------------------------------------------
/// Checks whether returned selection or selection name is not an
/// empty string.
/// \param string - string to be validated.
void SFDrawCommands::CheckCommand(TString string)
{
    if (string == "" || string == " ")
    {
        std::cerr << "##### Error in SFDrawCommands::GetSelection() or GetSelectionName()!"
                  << std::endl;
        std::cerr << "Empty string!" << std::endl;
        std::abort();
    }

    return;
}
//------------------------------------------------------------------
/// Returns name of chosen selection as TString. Useful when constructing
/// names of histograms.
/// \param selection - selection type.
TString SFDrawCommands::GetSelectionName(SFSelectionType selection)
{
    TString selectionName = "";

    switch (selection)
    {
        case SFSelectionType::kDDPE:
            selectionName = "DDPE";
            break;
        case SFSelectionType::kDDCharge:
            selectionName = "DDCharge";
            break;
        case SFSelectionType::kDDAmplitude:
            selectionName = "DDAmplitude";
            break;
        case SFSelectionType::kDDT0:
            selectionName = "DDT0";
            break;
        case SFSelectionType::kDDTOT:
            selectionName = "DDTOT";
            break;
        case SFSelectionType::kDDLogSqrtPERatio:
            selectionName = "DDLogSqrtPERatio";
            break;
        case SFSelectionType::kDDT0Difference:
            selectionName = "DDT0Difference";
            break;
        case SFSelectionType::kDDPEAverage:
            selectionName = "DDPEAverage";
            break;
        case SFSelectionType::kDDAmplitudeAverage:
            selectionName = "DDAmplitudeAverage";
            break;
        case SFSelectionType::kDDPECorrelation:
            selectionName = "DDPECorrelation";
            break;
        case SFSelectionType::kDDAmplitudeCorrelation:
            selectionName = "DDAmplitudeCorrelation";
            break;
        case SFSelectionType::kDDAmpPECorrelation:
            selectionName = "DDAmpPECorrelation";
            break;
        case SFSelectionType::kDDT0Correlation:
            selectionName = "DDT0Correlation";
            break;
        case SFSelectionType::kDDPEAttCorrected:
            selectionName = "DDPEAttCorrected";
            break;
        case SFSelectionType::kDDPEAttCorrectedSum:
            selectionName = "DDPEAttCorrectedSum";
            break;
        case SFSelectionType::kDDBL:
            selectionName = "DDBL";
            break;
        case SFSelectionType::kDDBLSigma:
            selectionName = "DDBLSigma";
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetSelectionName()!" << std::endl;
            std::cerr << "Unknown selection type! Please check!" << std::endl;
            break;
    }

    CheckCommand(selectionName);
    return selectionName;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object.
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param chAddr - channel address
/// \param customNum - standard vector containing values to be inserted in the selection.
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique, 
                                     SFChAddr chAddr, std::vector<double> customNum)
{
    TString selectionString = "";
    char    side            = chAddr.fSide;

    switch (selection)
    {
        case SFSelectionType::kDDPE:
            selectionString = Form("SDDSamples.data.signal_%c.fPE>>htemp%i(2200, -150, 1500)", side, unique);
            break;
        case SFSelectionType::kDDCharge:
            selectionString = Form("SDDSamples.data.signal_%c.fCharge>>htemp%i(1000, -1E4, 2.5E5)", side, unique);
            break;
        case SFSelectionType::kDDAmplitude:
            selectionString = Form("SDDSamples.data.signal_%c.fAmp>>htemp%i(800, 0, 800)", side, unique);
            break;
        case SFSelectionType::kDDT0:
            selectionString = Form("SDDSamples.data.signal_%c.fT0>>htemp%i(1210, -110, 1100)", side, unique);
            break;
        case SFSelectionType::kDDTOT:
            selectionString = Form("SDDSamples.data.signal_%c.fTOT>>htemp%i(1210, -110, 1100)", side, unique);
            break;
        case SFSelectionType::kDDPEAttCorrected:
            selectionString = Form("SDDSamples.data.signal_%c.fPE/exp(%f/%f)>>htemp%i(1300, -150, 1600)",
                                   side, customNum[0], customNum[1], unique);
            break;
        case SFSelectionType::kDDAmpPECorrelation:
            selectionString = Form("SDDSamples.data.signal_%c.fAmp:SDDSamples.data.signal_%c.fPE>>htemp%i"
                                   "(2200, -150, 1500, 1000, -10, 800)", side, side, unique);
            break;
        case SFSelectionType::kDDBL:
            selectionString = Form("SDDSamples.data.signal_%c.fBL>>htemp%i(2500, 1000, 3500)", side, unique);
            break;
        case SFSelectionType::kDDBLSigma:
            selectionString = Form("SDDSamples.data.signal_%c.fBL_sigma>>htemp%i(200, 0, 50)", side, unique);
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
            std::cerr << "Unknown selection type! Please check!" << std::endl;
            break;
    }

    CheckCommand(selectionString);
    std::cout << "\tUsing selection: " << selectionString << std::endl;

    return selectionString;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object.
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param customNum - standard vector containing values to be inserted in the selection.
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique,
                                     std::vector<double> customNum)
{
    TString selectionString = "";

    switch (selection)
    {
        case SFSelectionType::kDDLogSqrtPERatio:
            selectionString = Form("log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))>>"
                                   "htemp%i(500, -2, 2)", unique);
            break;
        case SFSelectionType::kDDT0Difference:
            selectionString = Form("(SDDSamples.data.signal_l.fT0-SDDSamples.data.signal_r.fT0)>>"
                                   "htemp%i(2500, -50, 50)", unique);
            break;
        case SFSelectionType::kDDPEAverage:
            selectionString = Form("sqrt(SDDSamples.data.signal_l.fPE*SDDSamples.data.signal_r.fPE)>>"
                                   "htemp%i(1350, -150, 1200)", unique);
            break;
        case SFSelectionType::kDDAmplitudeAverage:
            selectionString = Form("sqrt(SDDSamples.data.signal_l.fAmp*SDDSamples.data.signal_r.fAmp)>>"
                                   "htemp%i(800, 0, 800)", unique);
            break;
        case SFSelectionType::kDDPECorrelation:
            selectionString = Form("SDDSamples.data.signal_l.fPE:SDDSamples.data.signal_r.fPE>>"
                                   "htemp%i(3300, -150, 1500, 3300, -150, 1500)", unique);
            break;
        case SFSelectionType::kDDAmplitudeCorrelation:
            selectionString = Form("SDDSamples.data.signal_l.fAmp:SDDSamples.data.signal_r.fAmp>>"
                                   "htemp%i(800, 0, 800, 800, 0, 800)", unique);
            break;
        case SFSelectionType::kDDT0Correlation:
            selectionString = Form("SDDSamples.data.signal_l.fT0:SDDSamples.data.signal_r.fT0>>"
                                   "htemp%i(2420, -110, 1100, 2420, -110, 1100", unique);
            break;
        case SFSelectionType::kDDPEAttCorrectedSum:
            selectionString = Form("SDDSamples.data.signal_l.fPE/exp(%f/%f) + SDDSamples.data."
                                   "signal_r.fPE/exp(%f/%f)>>htemp%i(1500, -150, 4000)",
                                   customNum[0], customNum[1], customNum[2], customNum[3], unique);
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
            std::cerr << "Unknown selection type! Please check!" << std::endl;
            break;
    }

    CheckCommand(selectionString);
    std::cout << "\tUsing selection: " << selectionString << std::endl;

    return selectionString;
}
//------------------------------------------------------------------
/// Returns cut as a TString for ROOT's TTree type object.
/// \param cut - cut type
/// \param chAddr - channel address
/// \param customNum - standard vector containing numbers to be inserted in the cut
TString SFDrawCommands::GetCut(SFCutType cut, SFChAddr chAddr, std::vector<double> customNum)
{
    TString cutString = "";
    
    Int_t  mod  = chAddr.fModule;
    Int_t  lay  = chAddr.fLayer;
    Int_t  fib  = chAddr.fFiber;
    Char_t side = chAddr.fSide;

    switch (cut)
    {
        case SFCutType::kDDSpec:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i && "
                             "SDDSamples.data.signal_%c.fVeto == 0 && "
                             "SDDSamples.data.signal_%c.fBL_sigma<%f && "
                             "SDDSamples.data.signal_%c.fPE>0 && "
                             "SDDSamples.data.signal_%c.fT0>0 && "
                             "SDDSamples.data.signal_%c.fT0<590 && "
                             "SDDSamples.data.signal_%c.fTOT>0 && "
                             "SDDSamples.data.signal_%c.fAmp<%f",
                             mod, lay, fib, side, side, customNum[0], 
                             side, side, side, side, side, ampMax);
            break;
        case SFCutType::kDDSpecA:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i && "
                             "SDDSamples.data.signal_%c.fVeto == 0 && "
                             "SDDSamples.data.signal_%c.fBL_sigma<%f && "
                             "SDDSamples.data.signal_%c.fPE>0 && "
                             "SDDSamples.data.signal_%c.fT0>0 && "
                             "SDDSamples.data.signal_%c.fT0<590 && "
                             "SDDSamples.data.signal_%c.fTOT>0",
                              mod, lay, fib, side, side, customNum[0],
                              side, side, side, side);
            break;
        case SFCutType::kDDCombLR:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i && "
                             "SDDSamples.data.signal_l.fVeto == 0 && "
                             "SDDSamples.data.signal_r.fVeto == 0 && "
                             "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_r.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>0 && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0 && "
                             "SDDSamples.data.signal_r.fPE>0 && "
                             "SDDSamples.data.signal_r.fT0>0 && "
                             "SDDSamples.data.signal_r.fT0<590 && "
                             "SDDSamples.data.signal_r.fTOT>0 && "
                             "SDDSamples.data.signal_l.fAmp<%f && "
                             "SDDSamples.data.signal_r.fAmp<%f",
                             mod, lay, fib, customNum[0], customNum[1], ampMax, ampMax);
            break;
        case SFCutType::kDDT0Diff:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i && "
                             "SDDSamples.data.signal_l.fVeto == 0 && "
                             "SDDSamples.data.signal_r.fVeto == 0 && "
                             "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_r.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>%f && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0 && "
                             "SDDSamples.data.signal_r.fPE>%f && "
                             "SDDSamples.data.signal_r.fT0>0 && "
                             "SDDSamples.data.signal_r.fT0<590 && "
                             "SDDSamples.data.signal_r.fTOT>0 && "
                             "SDDSamples.data.signal_l.fAmp<%f && "
                             "SDDSamples.data.signal_r.fAmp<%f && "
                             "log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))>%f && "
                             "log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))<%f",
                             mod, lay, fib, customNum[0], customNum[1], customNum[2],  
                             customNum[3], ampMax, ampMax, customNum[4], customNum[5]);
            break;
        case SFCutType::kDDT0DiffECut:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i && "
                             "SDDSamples.data.signal_l.fVeto == 0 && "
                             "SDDSamples.data.signal_r.fVeto == 0 && "
                             "SDDSamples.data.signal_l.fAmp<%f && "
                             "SDDSamples.data.signal_r.fAmp<%f && "
                             "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_r.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>%f && "
                             "SDDSamples.data.signal_l.fPE<%f && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0 && "
                             "SDDSamples.data.signal_r.fPE>%f && "
                             "SDDSamples.data.signal_r.fPE<%f && "
                             "SDDSamples.data.signal_r.fT0>0 && "
                             "SDDSamples.data.signal_r.fT0<590 && "
                             "SDDSamples.data.signal_r.fTOT>0 && "
                             "log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))>%f && "
                             "log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))<%f",
                             mod, lay, fib, ampMax, ampMax, customNum[0], customNum[1], customNum[2],
                             customNum[3], customNum[4], customNum[5], customNum[6], customNum[7]);
            break;
        case SFCutType::kDDBL:
            cutString = Form("SDDSamples.data.module == %i && "
                             "SDDSamples.data.layer == %i && "
                             "SDDSamples.data.fiber == %i",
                             mod, lay, fib);
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetCut()!" << std::endl;
            std::cerr << "Unknown cut type! Please check!" << std::endl;
            break;
    }

    std::cout << "\tUsing cut: " << cutString << std::endl;

    return cutString;
}
//------------------------------------------------------------------
