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

ClassImp(SFDrawCommands);

const double ampMax = 660; // maximum valid amplitude [mV] in measurements 
                           // with the Desktop Digitizer; signals with the amplitude
                           // above this value are saturated and need to be discarded

//------------------------------------------------------------------
/// Returns full address of the requested channel, according to the convention
/// from sifi-framework. Address is returned as ChannelAddress object, including
/// address, channel ID, module, layer, fiber and side.
/// \param ch - channel number
ChannelAddress SFDrawCommands::GetChannelAddress(int ch)
{
    ChannelAddress chAddr;

    switch (ch)
    {
        case 0:
            chAddr.fAddress = 0x1000;
            chAddr.fChID    = 0;
            chAddr.fModule  = 0;
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'l';
            break;
        case 1:
            chAddr.fAddress = 0x1001;
            chAddr.fChID    = 1;
            chAddr.fModule  = 0;
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'r';
            break;
        case 2:
            chAddr.fAddress = 0x1002;
            chAddr.fChID    = 2;
            chAddr.fModule  = 1;
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'l';
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetChannelAddress()! " << std::endl;
            std::cerr << "Incorrect channel number: " << ch << std::endl;
            std::cerr << "Correct channel numbers: 0, 1 & 2" << std::endl;
            break;
    }

    return chAddr;
}
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
        case SFSelectionType::kPE:
            selectionName = "PE";
            break;
        case SFSelectionType::kCharge:
            selectionName = "Charge";
            break;
        case SFSelectionType::kAmplitude:
            selectionName = "Amplitude";
            break;
        case SFSelectionType::kT0:
            selectionName = "T0";
            break;
        case SFSelectionType::kTOT:
            selectionName = "TOT";
            break;
        case SFSelectionType::kLogSqrtPERatio:
            selectionName = "LogSqrtPERatio";
            break;
        case SFSelectionType::kT0Difference:
            selectionName = "T0Difference";
            break;
        case SFSelectionType::kPEAverage:
            selectionName = "PEAverage";
            break;
        case SFSelectionType::kAmplitudeAverage:
            selectionName = "AmplitudeAverage";
            break;
        case SFSelectionType::kPECorrelation:
            selectionName = "PECorrelation";
            break;
        case SFSelectionType::kAmplitudeCorrelation:
            selectionName = "AmplitudeCorrelation";
            break;
        case SFSelectionType::kAmpPECorrelation:
            selectionName = "AmpPECorrelation";
            break;
        case SFSelectionType::kT0Correlation:
            selectionName = "T0Correlation";
            break;
        case SFSelectionType::kPEAttCorrected:
            selectionName = "PEAttCorrected";
            break;
        case SFSelectionType::kPEAttCorrectedSum:
            selectionName = "PEAttCorrectedSum";
            break;
        case SFSelectionType::kBL:
            selectionName = "BL";
            break;
        case SFSelectionType::kBLSigma:
            selectionName = "BLSigma";
            break;
        //----- selections for PMI measurements
        case SFSelectionType::kPMICharge:
            selectionName = "kPMICharge";
            break;
        case SFSelectionType::kPMIChargeAverage:
            selectionName = "kPMITChargeAverage";
            break;
        case SFSelectionType::kPMIT0:
            selectionName = "kPMIT0";
            break;
        case SFSelectionType::kPMIChargeCorrelation:
            selectionName = "kPMIChargeCorrelation";
            break;
        case SFSelectionType::kPMIT0Correlation:
            selectionName = "kPMIT0Correlation";
            break;
        case SFSelectionType::kPMIT0Difference:
            selectionName = "kPMIT0Difference";
            break;
        case SFSelectionType::kPMILogSqrtChargeRatio:
            selectionName = "kPMILogSqrtChargeRatio";
            break;
        case SFSelectionType::kPMIChargeAttCorrected:
            selectionName = "kPMIChargeAttCorrected";
            break;
        case SFSelectionType::kPMIChargeAttCorrectedSum:
            selectionName = "kPMIChargeAttCorrectedSum";
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
/// \param ch - channel number
/// \param customNum - standard vector containing values to be inserted in the selection.
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique, int ch,
                                     std::vector<double> customNum)
{
    TString        selectionString = "";
    ChannelAddress chAddr          = GetChannelAddress(ch);
    char           side            = chAddr.fSide;

    switch (selection)
    {
        case SFSelectionType::kPE:
            selectionString = Form("SDDSamples.data.signal_%c.fPE>>htemp%i(2200, -150, 1500)", side, unique);
//             selectionString = Form("SDDSamples.data.signal_%c.fPE>>htemp%i(2500, -150, 6E4)", side, unique);
            break;
        case SFSelectionType::kCharge:
            selectionString = Form("SDDSamples.data.signal_%c.fCharge>>htemp%i(1000, -1E4, 2.5E5)", side, unique);
            break;
        case SFSelectionType::kAmplitude:
            selectionString = Form("SDDSamples.data.signal_%c.fAmp>>htemp%i(800, 0, 800)", side, unique);
            break;
        case SFSelectionType::kT0:
            selectionString = Form("SDDSamples.data.signal_%c.fT0>>htemp%i(1210, -110, 1100)", side, unique);
            break;
        case SFSelectionType::kTOT:
            selectionString = Form("SDDSamples.data.signal_%c.fTOT>>htemp%i(1210, -110, 1100)", side, unique);
            break;
        case SFSelectionType::kPEAttCorrected:
            selectionString = Form("SDDSamples.data.signal_%c.fPE/exp(%f/%f)>>htemp%i(1300, -150, 1600)",
                                   side, customNum[0], customNum[1], unique);
            break;
        case SFSelectionType::kAmpPECorrelation:
            selectionString = Form("SDDSamples.data.signal_%c.fAmp:SDDSamples.data.signal_%c.fPE>>htemp%i"
                                   "(2200, -150, 1500, 1000, -10, 800)", side, side, unique);
            break;
        case SFSelectionType::kBL:
            selectionString = Form("SDDSamples.data.signal_%c.fBL>>htemp%i(2500, 1000, 3500)", side, unique);
            break;
        case SFSelectionType::kBLSigma:
            selectionString = Form("SDDSamples.data.signal_%c.fBL_sigma>>htemp%i(200, 0, 50)", side, unique);
            break;
        //----- selections for PMI measurements
        case SFSelectionType::kPMICharge:
            selectionString = Form("SFibersRaw.data.qdc_%c>>htemp%i(440, 0, 1500)", side, unique);
            break;
        case SFSelectionType::kPMIT0:
            selectionString = Form("SFibersRaw.data.time_%c>>htemp%i(300, -50, 100)", side, unique);
            break;
        case SFSelectionType::kPMIChargeAttCorrected:
            selectionString = Form("SFibersRaw.data.qdc_%c/exp(%f/%f)>>htemp%i(500, 0, 1600)",
                                   side, customNum[0], customNum[1], unique);
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
        case SFSelectionType::kLogSqrtPERatio:
            selectionString = Form("log(sqrt(SDDSamples.data.signal_r.fPE/SDDSamples.data.signal_l.fPE))>>"
                                   "htemp%i(500, -2, 2)", unique);
            break;
        case SFSelectionType::kT0Difference:
            selectionString = Form("(SDDSamples.data.signal_l.fT0-SDDSamples.data.signal_r.fT0)>>"
                                   "htemp%i(2500, -50, 50)", unique);
            break;
        case SFSelectionType::kPEAverage:
            selectionString = Form("sqrt(SDDSamples.data.signal_l.fPE*SDDSamples.data.signal_r.fPE)>>"
                                   "htemp%i(1350, -150, 1200)", unique); 
//             selectionString = Form("sqrt(SDDSamples.data.signal_l.fPE*SDDSamples.data.signal_r.fPE)>>"
//                                    "htemp%i(2500, -150, 6E4)", unique);
            break;
        case SFSelectionType::kAmplitudeAverage:
            selectionString = Form("sqrt(SDDSamples.data.signal_l.fAmp*SDDSamples.data.signal_r.fAmp)>>"
                                   "htemp%i(800, 0, 800)", unique);
            break;
        case SFSelectionType::kPECorrelation:
            selectionString = Form("SDDSamples.data.signal_l.fPE:SDDSamples.data.signal_r.fPE>>"
                                   "htemp%i(3300, -150, 1500, 3300, -150, 1500)", unique);
            break;
        case SFSelectionType::kAmplitudeCorrelation:
            selectionString = Form("SDDSamples.data.signal_l.fAmp:SDDSamples.data.signal_r.fAmp>>"
                                   "htemp%i(800, 0, 800, 800, 0, 800)", unique);
            break;
        case SFSelectionType::kT0Correlation:
            selectionString = Form("SDDSamples.data.signal_l.fT0:SDDSamples.data.signal_r.fT0>>"
                                   "htemp%i(2420, -110, 1100, 2420, -110, 1100", unique);
            break;
        case SFSelectionType::kPEAttCorrectedSum:
            selectionString = Form("SDDSamples.data.signal_l.fPE/exp(%f/%f) + SDDSamples.data."
                                   "signal_r.fPE/exp(%f/%f)>>htemp%i(1500, -150, 4000)",
                                   customNum[0], customNum[1], customNum[2], customNum[3], unique);
            break;
        //----- selections for PMI measurements
        case SFSelectionType::kPMIChargeAverage:
            selectionString = Form("sqrt(SFibersRaw.data.qdc_l*SFibersRaw.data.qdc_r)>>"
                                   "htemp%i(440, 0, 1500)", unique);
            break;
        case SFSelectionType::kPMIChargeCorrelation:
            selectionString = Form("SFibersRaw.data.qdc_l:SFibersRaw.data.qdc_r>>"
                                   "htemp%i(440, 0, 1500, 440, 0, 1500)", unique);
            break;
        case SFSelectionType::kPMIT0Correlation:
            selectionString = Form("SFibersRaw.data.time_l.:SFibersRaw.data.time_r>>"
                                   "htemp%i(300, -50, 50, 600, -50, 50", unique);
            break;
        case SFSelectionType::kPMIT0Difference:
            selectionString = Form("(SFibersRaw.data.time_l-SFibersRaw.data.time_r)>>"
                                   "htemp%i(500, -50, 50)", unique);
            break;
        case SFSelectionType::kPMILogSqrtChargeRatio:
            selectionString = Form("log(sqrt(SFibersRaw.data.qdc_r/SFibersRaw.data.qdc_l))>>"
                                   "htemp%i(200, -2, 2)", unique);
            break;
        case SFSelectionType::kPMIChargeAttCorrectedSum:
            selectionString = Form("SFibersRaw.data.qdc_l/exp(%f/%f) + SFibersRaw.data."
                                   "qdc_r/exp(%f/%f)>>htemp%i(500, 0, 4000)",
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
/// \param customNum - standard vector containing numbers to be inserted in the cut
TString SFDrawCommands::GetCut(SFCutType cut, std::vector<double> customNum)
{
    TString cutString = "";

    switch (cut)
    {
        case SFCutType::kPMISpecCh0:
            cutString = "SFibersRaw.data.module==0";
            break;
        case SFCutType::kPMISpecCh1:
            cutString = "SFibersRaw.data.module==0";
            break;
        case SFCutType::kPMICombCh0Ch1:
            cutString = "SFibersRaw.data.module==0";
            break;
        case SFCutType::kPMIT0Diff:
            cutString = Form("SFibersRaw.data.module==0 && "
                             "log(sqrt(SFibersRaw.data.qdc_r/SFibersRaw.data.qdc_l))>%f && "
                             "log(sqrt(SFibersRaw.data.qdc_r/SFibersRaw.data.qdc_l))<%f", 
                             customNum[0], customNum[1]);
            break;
        case SFCutType::kPMIT0DiffECut:
            cutString = Form("SFibersRaw.data.module==0 && "
                             "SFibersRaw.data.qdc_l>%f && "
                             "SFibersRaw.data.qdc_l<%f && "
                             "SFibersRaw.data.qdc_r>%f && "
                             "SFibersRaw.data.qdc_r<%f && "
                             "log(sqrt(SFibersRaw.data.qdc_r/SFibersRaw.data.qdc_l))>%f && " 
                             "log(sqrt(SFibersRaw.data.qdc_r/SFibersRaw.data.qdc_l))<%f",
                             customNum[0], customNum[1], customNum[2], 
                             customNum[3], customNum[4], customNum[5]);
            break;
        case SFCutType::kSpecCh0:
            cutString = Form("SDDSamples.data.module==0 && "
//                              "SDDSamples.data.signal_l.fVeto == 0 && "
//                              "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>0 && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0 && "
                             "SDDSamples.data.signal_l.fAmp<%f",
                             /*customNum[0],*/ ampMax);
            break;
        case SFCutType::kSpecCh0A:
            cutString = Form("SDDSamples.data.module==0 && "
//                              "SDDSamples.data.signal_l.fVeto == 0 && "
//                              "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>0 && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0"/*,
                             customNum[0]*/);
            break;
        case SFCutType::kSpecCh1:
            cutString = Form("SDDSamples.data.module==0 && "
//                              "SDDSamples.data.signal_r.fVeto == 0 && "
//                              "SDDSamples.data.signal_r.fBL_sigma<%f && "
                             "SDDSamples.data.signal_r.fPE>0 && "
                             "SDDSamples.data.signal_r.fT0>0 && "
                             "SDDSamples.data.signal_r.fT0<590 && "
                             "SDDSamples.data.signal_r.fTOT>0 && "
                             "SDDSamples.data.signal_r.fAmp<%f",
                             /*customNum[0],*/ ampMax);
            break;
        case SFCutType::kSpecCh1A:
            cutString = Form("SDDSamples.data.module==0 && "
//                              "SDDSamples.data.signal_r.fVeto == 0 && "
//                              "SDDSamples.data.signal_r.fBL_sigma<%f && "
                             "SDDSamples.data.signal_r.fPE>0 && "
                             "SDDSamples.data.signal_r.fT0>0 && "
                             "SDDSamples.data.signal_r.fT0<590 && "
                             "SDDSamples.data.signal_r.fTOT>0"/*,
                             customNum[0]*/);
            break;
        case SFCutType::kSpecCh2:
            cutString = Form("SDDSamples.data.module==1 && "
                             "SDDSamples.data.signal_l.fVeto == 0 && "
                             "SDDSamples.data.signal_l.fBL_sigma<%f && "
                             "SDDSamples.data.signal_l.fPE>0 && "
                             "SDDSamples.data.signal_l.fT0>0 && "
                             "SDDSamples.data.signal_l.fT0<590 && "
                             "SDDSamples.data.signal_l.fTOT>0",
                             customNum[0]);
            break;
        case SFCutType::kCombCh0Ch1:
            cutString = Form("SDDSamples.data.module==0 && "
//                              "SDDSamples.data.signal_l.fVeto == 0 && "
//                              "SDDSamples.data.signal_r.fVeto == 0 && "
//                              "SDDSamples.data.signal_l.fBL_sigma<%f && "
//                              "SDDSamples.data.signal_r.fBL_sigma<%f && "
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
                             /*customNum[0], customNum[1],*/ ampMax, ampMax);
            break;
        case SFCutType::kT0Diff:
            cutString =
                Form("SDDSamples.data.module==0 && "
                     //"SDDSamples.data.signal_l.fVeto == 0 && "
                     //"SDDSamples.data.signal_r.fVeto == 0 && "
                     //"SDDSamples.data.signal_l.fBL_sigma<%f && "
                     //"SDDSamples.data.signal_r.fBL_sigma<%f && "
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
                     /*customNum[0], customNum[1],*/ customNum[2], customNum[3], ampMax, ampMax,
                     customNum[4], customNum[5]);
            break;
        case SFCutType::kT0DiffECut:
            cutString =
                Form("SDDSamples.data.module==0 && "
                     //"SDDSamples.data.signal_l.fVeto == 0 && "
                     //"SDDSamples.data.signal_r.fVeto == 0 && "
                     "SDDSamples.data.signal_l.fAmp<%f && "
                     "SDDSamples.data.signal_r.fAmp<%f && "
                     //"SDDSamples.data.signal_l.fBL_sigma<%f && "
                     //"SDDSamples.data.signal_r.fBL_sigma<%f && "
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
                     ampMax, ampMax, /*customNum[0], customNum[1],*/ customNum[2], customNum[3],
                     customNum[4], customNum[5], customNum[6], customNum[7]);
            break;
        case SFCutType::kBLCh0:
            cutString = "SDDSamples.data.module==0";
            break;
        case SFCutType::kBLCh1:
            cutString = "SDDSamples.data.module==0";
            break;
        case SFCutType::kBLCh2:
            cutString = "SDDSamples.data.module==1";
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
/// Prints details of the SFDrawCommands class object.
void SFDrawCommands::Print(void)
{
    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "This is print out of SFSelection class object" << std::endl;
    std::cout << "--------------------------------------------------\n" << std::endl;
    return;
}
//------------------------------------------------------------------
