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

const double ampMax = 660;

//------------------------------------------------------------------
/// Returns full address of the requested channel, according to the convention
/// from sifi-framework. Address is returned as ChannelAddress object, including
/// address, channel ID, module, layer, fiber and side.
/// \param ch - channel number
ChannelAddress SFDrawCommands::GetChannelAddress(int ch)
{
    ChannelAddress chAddr;
    
    switch(ch){
        case 0:
            chAddr.fAddress = 0x1000;
            chAddr.fChID    = 0;
            chAddr.fModule  = 0; 
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'l';
            break;
        case 1:
            chAddr.fAddress = 0x1000;
            chAddr.fChID    = 1;
            chAddr.fModule  = 0; 
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'r';
            break;
        case 2:
            chAddr.fAddress = 0x1000;
            chAddr.fChID    = 2;
            chAddr.fModule  = 1; 
            chAddr.fLayer   = 0;
            chAddr.fFiber   = 0;
            chAddr.fSide    = 'l';
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
  if (string=="" || string==" ")
  {
    std::cerr << "##### Error in SFDrawCommands::GetSelection() or GetSelectionName()!" << std::endl;
    std::cerr << "Empty selection!" << std::endl;
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
        case SFSelectionType::PE:
            selectionName = "PE";
            break;
        case SFSelectionType::Charge:
            selectionName = "Charge";
            break;
        case SFSelectionType::Amplitude:
            selectionName = "Amplitude";
            break;
        case SFSelectionType::T0:
            selectionName = "T0";
            break;
        case SFSelectionType::TOT:
            selectionName = "TOT";
            break;
        case SFSelectionType::LogSqrtPERatio:
            selectionName = "LogSqrtPERatio";
            break;
        case SFSelectionType::T0Difference:
            selectionName = "T0Difference";
            break;
        case SFSelectionType::PEAverage:
            selectionName = "PEAverage";
            break;
        case SFSelectionType::AmplitudeAverage:
            selectionName = "AmplitudeAverage";
            break;
        case SFSelectionType::PECorrelation:
            selectionName = "PECorrelation";
            break;
        case SFSelectionType::AmplitudeCorrelation:
            selectionName = "AmplitudeCorrelation";
            break;
        case SFSelectionType::AmpPECorrelation:
            selectionName = "AmpPECorrelation";
            break;
        case SFSelectionType::T0Correlation:
            selectionName = "T0Correlation";
            break;
        case SFSelectionType::PEAttCorrected:
            selectionName = "PEAttCorrected";
            break;
        case SFSelectionType::PEAttCorrectedSum:
            selectionName = "PEAttCorrectedSum";
            break;
        case SFSelectionType::BL:
            selectionName = "BL";
            break;
        case SFSelectionType::BLSigma:
            selectionName = "BLSigma";
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
                                     std::vector <double> customNum)
{
  TString selectionString = "";
  ChannelAddress chAddr = GetChannelAddress(ch);
  char side = chAddr.fSide;
   
  switch(selection){
      case SFSelectionType::PE:
          selectionString = Form("SFibersStackCal.data.qdc_%c>>htemp%i(2200, -150, 1500)", side, unique);
          break;
      case SFSelectionType::Charge:
          selectionString = Form("SFibersStackRaw.data.qdc_%c>>htemp%i(1000, -1E4, 2.5E5)", side, unique);
          break;
      case SFSelectionType::Amplitude:
          selectionString = Form("SDDSamples.data.signal_%c.fAmp>>htemp%i(800, 0, 800)", side, unique);
          break;
      case SFSelectionType::T0:
          selectionString = Form("SFibersStackCal.data.time_%c.fT0>>htemp%i(1210, -110, 1100)", side, unique);
          break;
      case SFSelectionType::TOT:
          selectionString = Form("SDDSamples.data.signal_%c.fTOT>>htemp%i(1210, -110, 1100)", side, unique);
          break;
      case SFSelectionType::PEAttCorrected:
          selectionString = Form("SFibersStackCal.data.qdc_%c/exp(%f/%f)>>htemp%i(1300, -150, 1600)", 
                            side, customNum[0], customNum[1], unique);
          break;
      case SFSelectionType::AmpPECorrelation:
          selectionString = Form("SDDSamples.data.signal_%c.fAmp:SFibersStackCal.data.qdc_%c>>htemp%i(2200, -150, 1500, 1000, -10, 800)",
                            side, side, unique);
          break;
      case SFSelectionType::BL:
          selectionString = Form("SDDSamples.data.signal_%c.fBL>>htemp%i(2500, 1000, 3500)", side, unique);
          break;
      case SFSelectionType::BLSigma:
          selectionString = Form("SDDSamples.data.signal_%c.fBL_sigma>>htemp%i(200, 0, 50)", side, unique);
          break;
      default:
          std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
          std::cerr << "Unknown selection type! Please check!" << std::endl;
          break;
  }
    
  CheckCommand(selectionString);
  std::cout << "\n\tUsing selection: " << selectionString << std::endl;
  
  return selectionString;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object. 
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param customNum - standard vector containing values to be inserted in the selection.
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique,
                                     std::vector <double> customNum)
{
  TString selectionString = "";

  switch(selection){
      case SFSelectionType::LogSqrtPERatio:
          selectionString = Form("log(sqrt(SFibersStackCal.data.qdc_r/SFibersStackCal.data.qdc_l))>>htemp%i(500, -2, 2)", unique);
          break;
      case SFSelectionType::T0Difference:
          selectionString = Form("(SFibersStackCal.data.time_l-SFibersStackCal.data.time_r)>>htemp%i(2500, -50, 50)", unique);
          break;
      case SFSelectionType::PEAverage:
          selectionString = Form("sqrt(SFibersStackCal.data.qdc_l*SFibersStackCal.data.qdc_r)>>htemp%i(1350, -150, 1200)", unique);
          break;
      case SFSelectionType::AmplitudeAverage:
          selectionString = Form("sqrt(SDDSamples.data.signal_l.fAmp*SDDSamples.data.signal_r.fAmp)>>htemp%i(800, 0, 800)",unique);
          break;
      case SFSelectionType::PECorrelation:
          selectionString = Form("SFibersStackCal.data.qdc_l:SFibersStackCal.data.qdc_r>>htemp%i(3300, -150, 1500, 3300, -150, 1500)", unique);
          break;
      case SFSelectionType::AmplitudeCorrelation:
          selectionString = Form("SDDSamples.data.signal_l.fAmp:SDDSamples.data.signal_r.fAmp>>htemp%i(800, 0, 800, 800, 0, 800)", unique);
          break;
      case SFSelectionType::T0Correlation:
          selectionString = Form("SFibersStackCal.data.time_l:SFibersStackCal.data.time_r>>htemp%i(2420, -110, 1100, 2420, -110, 1100", unique);
          break;
      case SFSelectionType::PEAttCorrectedSum:
          selectionString = Form("SFibersStackCal.data.qdc_l/exp(%f/%f) + SFibersStackCal.data.qdc_r/exp(%f/%f)>>htemp%i(1500, -150, 4000)",
                            customNum[0], customNum[1], customNum[2], customNum[3], unique);
          break;
      default:
          std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
          std::cerr << "Unknown selection type! Please check!" << std::endl;
          break;
  }
    
  CheckCommand(selectionString);
  std::cout << "\n\tUsing selection: " << selectionString << std::endl;
  
  return selectionString;
}
//------------------------------------------------------------------
/// Returns cut as a TString for ROOT's TTree type object.
/// \param cut - cut type
/// \param customNum - standard vector containing numbers to be inserted in the cut 
TString SFDrawCommands::GetCut(SFCutType cut, std::vector <double> customNum)
{
    TString cutString = "";
    
    switch(cut){
        case SFCutType::SpecCh0:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_l.fBL_sigma<%f &&"
                             "SFibersStackCal.data.qdc_l>0 &&"
                             "SFibersStackCal.data.time_l>0 &&"
                             "SFibersStackCal.data.time_l<590 &&" "SDDSamples.data.signal_l.fTOT>0 &&" "SDDSamples.data.signal_l.fAmp<%f",
                             customNum[0], ampMax);
            break;
        case SFCutType::SpecCh0A:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_l.fBL_sigma<%f &&"
                             "SFibersStackCal.data.qdc_l>0 &&"
                             "SFibersStackCal.data.time_l>0 &&"
                             "SFibersStackCal.data.time_l<590 &&"
                             "SDDSamples.data.signal_l.fTOT>0", 
                             customNum[0]);
            break;
        case SFCutType::SpecCh1:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_r.fBL_sigma<%f &&"
                             "SFibersStackCal.data.qdc_r>0 &&"
                             "SFibersStackCal.data.time_r>0 &&"
                             "SFibersStackCal.data.time_r<590 &&"
                             "SDDSamples.data.signal_r.fTOT>0 &&"
                             "SDDSamples.data.signal_r.fAmp<%f", 
                             customNum[0], ampMax);
            break;
        case SFCutType::SpecCh1A:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_r.fBL_sigma<%f &&"
                             "SFibersStackCal.data.qdc_r>0 &&"
                             "SFibersStackCal.data.time_r>0 &&"
                             "SFibersStackCal.data.time_r<590 &&"
                             "SDDSamples.data.signal_r.fTOT>0", 
                             customNum[0]);
            break;
        case SFCutType::SpecCh2:
            cutString = Form("SDDSamples.data.module==1 && SDDSamples.data.signal_l.fBL_sigma<%f && SFibersStackCal.data.qdc_l>0 && SFibersStackCal.data.time_l>0 && SFibersStackCal.data.time_l<590 && SDDSamples.data.signal_l.fTOT>0", customNum[0]);
            break;
        case SFCutType::CombCh0Ch1:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_l.fBL_sigma<%f &&"
                             "SDDSamples.data.signal_r.fBL_sigma<%f &&" "SFibersStackCal.data.qdc_l>0 &&"
                             "SFibersStackCal.data.time_l>0 &&"
                             "SFibersStackCal.data.time_l<590 &&"
                             "SDDSamples.data.signal_l.fTOT>0 &&"
                             "SFibersStackCal.data.qdc_r>0 &&"
                             "SFibersStackCal.data.time_r>0 &&"
                             "SFibersStackCal.data.time_r<590 &&"
                             "SDDSamples.data.signal_r.fTOT>0 &&"
                             "SDDSamples.data.signal_l.fAmp<%f &&"
                             "SDDSamples.data.signal_r.fAmp<%f", 
                             customNum[0], customNum[1], ampMax, ampMax);
            break;
        case SFCutType::T0Diff:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_l.fBL_sigma<%f &&"
                             "SDDSamples.data.signal_r.fBL_sigma<%f &&"
                             "SFibersStackCal.data.qdc_l>%f &&"
                             "SFibersStackCal.data.time_l>0 &&"
                             "SFibersStackCal.data.time_l<590 &&"
                             "SDDSamples.data.signal_l.fTOT>0 &&"
                             "SFibersStackCal.data.qdc_r>%f &&"
                             "SFibersStackCal.data.time_r>0 &&"
                             "SFibersStackCal.data.time_r<590 &&"
                             "SDDSamples.data.signal_r.fTOT>0 &&"
                             "SDDSamples.data.signal_l.fAmp<%f &&"
                             "SDDSamples.data.signal_r.fAmp<%f &&"
                             "log(sqrt(SFibersStackCal.data.qdc_r/SFibersStackCal.data.qdc_l))>%f &&"
                             "log(sqrt(SFibersStackCal.data.qdc_r/SFibersStackCal.data.qdc_l))<%f",
                             customNum[0], customNum[1], customNum[2], customNum[3],
                             ampMax, ampMax, customNum[4], customNum[5]);
            break;
        case SFCutType::T0DiffECut:
            cutString = Form("SDDSamples.data.module==0 &&"
                             "SDDSamples.data.signal_l.fAmp<%f &&"
                             "SDDSamples.data.signal_r.fAmp<%f &&"
                             "SDDSamples.data.signal_l.fBL_sigma<%f &&"
                             "SDDSamples.data.signal_r.fBL_sigma<%f &&" "SFibersStackCal.data.qdc_l>%f &&"
                             "SFibersStackCal.data.qdc_l.<%f &&"
                             "SFibersStackCal.data.time_l>0 &&"
                             "SFibersStackCal.data.time_l<590 &&"
                             "SDDSamples.data.signal_l.fTOT>0 &&"
                             "SFibersStackCal.data.qdc_r>%f &&"
                             "SFibersStackCal.data.qdc_r<%f &&"
                             "SFibersStackCal.data.time_r>0 &&"
                             "SFibersStackCal.data.time_r<590 &&"
                             "SDDSamples.data.signal_r.fTOT>0 &&" "log(sqrt(SFibersStackCal.data.qdc_r/ SFibersStackCal.data.qdc_l))>%f &&"
                             "log(sqrt(SFibersStackCal.data.qdc_r/SFibersStackCal.data.qdc_l))<%f",
                             ampMax, ampMax, customNum[0], customNum[1], 
                             customNum[2], customNum[3], customNum[4], 
                             customNum[5], customNum[6], customNum[7]);
            break;
        case SFCutType::BLCh0:
            cutString = "SDDSamples.data.module==0";
            break;
        case SFCutType::BLCh1:
            cutString = "SDDSamples.data.module==0";
            break;
        case SFCutType::BLCh2:
            cutString = "SDDSamples.data.module==1";
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetCut()!" << std::endl;
            std::cerr << "Unknown cut type! Please check!" << std::endl;
            break;
    }
    
    std::cout << "\n\tUsing cut: " << cutString << std::endl;
    
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
