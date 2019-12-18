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

//------------------------------------------------------------------
/// Checks whether returned selection or selection name is not an 
/// empty string. 
void SFDrawCommands::CheckSelection(TString string){
    
  if(string=="" || string==" "){
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
TString SFDrawCommands::GetSelectionName(SFSelectionType selection){
 
    TString selectionName = "";
    
    switch(selection){
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
        case SFSelectionType::LogSqrtRatioCut:
            selectionName = "LogSqrtRatioCut";
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
        case SFSelectionType::T0Correlation:
            selectionName = "T0Correlation";
            break;
        case SFSelectionType::PEAttCorrected:
            selectionName = "PEAttCorrected";
            break;
        case SFSelectionType::PEAttCorrectedSum:
            selectionName = "PEAttCorrectedSum";
            break;
        default:
            std::cerr << "##### Error in SFDrawCommands::GetSelectionName()!" << std::endl;
            std::cerr << "Unknown selection type! Please check!" << std::endl;
            break;
    }

    CheckSelection(selectionName);
    return selectionName;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object. 
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param ch - channel number
/// \param customNum - standard vector containing values to be inserted in the selection
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique, int ch,
                                     std::vector <double> customNum){
    
  TString selectionString = "";
   
  switch(selection){
      case SFSelectionType::PE:
          selectionString = Form("ch_%i.fPE>>htemp%i(1000,-150,1500)", ch, unique);
          break;
      case SFSelectionType::Charge:
          selectionString = Form("ch_%i.fCharge>>htemp%i(1000,-1E4,2.5E5)", ch, unique);
          break;
      case SFSelectionType::Amplitude:
          selectionString = Form("ch_%i.fAmp>>htemp%i(1000,0,700)", ch, unique);
          break;
      case SFSelectionType::T0:
          selectionString = Form("ch_%i.fT0>>htemp%.i(1000,-110,1100)", ch, unique);
          break;
      case SFSelectionType::TOT:
          selectionString = Form("ch_%i.fTOT>>htemp%.i(1000,-110,1100)", ch, unique);
          break;
      case SFSelectionType::PEAttCorrected:
          selectionString = Form("ch_%i.fPE/exp(%f/%f)>>htemp%i(1300,-150,1600)", ch,
                            customNum[0], customNum[1], unique);
          break;
      default:
          std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
          std::cerr << "Unknown selection type! Please check!" << std::endl;
          break;
  }
    
  CheckSelection(selectionString);
  return selectionString;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object. 
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param customNum - standard vector containing values to be inserted in the selection
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique,
                                     std::vector <double> customNum){
    
  TString selectionString = "";
  TString LogSqrtQ1Q0 = "log(sqrt(ch_1.fPE/ch_0.fPE))";

  switch(selection){
      case SFSelectionType::LogSqrtPERatio:
          selectionString = Form("log(sqrt(ch_1.fPE/ch_0.fPE))>>htemp%i(500, -2, 2)", unique);
          break;
      case SFSelectionType::LogSqrtRatioCut:
          selectionString = Form("log(sqrt(ch_1.fPE/ch_0.fPE))>>htemp%i(500, -0.5, 0.5)", unique);
          break;
      case SFSelectionType::T0Difference:
          selectionString = Form("(ch_0.fT0-ch_1.fT0)>>htemp%i(500,-50,50)", unique);
          break;
      case SFSelectionType::PEAverage:
          selectionString = Form("sqrt(ch_0.fPE*ch_1.fPE)>>htemp%i(1000,-150,1200)", unique);
          break;
      case SFSelectionType::AmplitudeAverage:
          selectionString = Form("sqrt(ch_0.fAmp*ch_1.fAmp)>>htemp%i(1000,0,800)",unique);
          break;
      case SFSelectionType::PECorrelation:
          selectionString = Form("ch_0.fPE:ch_1.fPE>>htemp%i(1000,-150,1200,1000,-150,1200)", unique);
          break;
      case SFSelectionType::AmplitudeCorrelation:
          selectionString = Form("ch_0.fAmp:ch_1.fAmp>>htemp%i(1000,0,800,1000,0,800)", unique);
          break;
      case SFSelectionType::T0Correlation:
          selectionString = Form("ch_0.fT0:ch_1.fT0>>htemp%i(1000,-110,1100,1000,-110,1100", unique);
          break;
      case SFSelectionType::PEAttCorrectedSum:
          selectionString = Form("ch_0.fPE/exp(%f/%f) + ch_1.fPE/exp(%f/%f)>>htemp%i(1500,-150,4000)",
                            customNum[0], customNum[1], customNum[2], customNum[3],
                            unique);
          break;
      default:
          std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
          std::cerr << "Unknown selection type! Please check!" << std::endl;
          break;
  }
    
  CheckSelection(selectionString);
  return selectionString;
}
//------------------------------------------------------------------
/// Prints details of the SFDrawCommands class object.
void SFDrawCommands::Print(void){
  std::cout << "\n------------------------------------------------" << std::endl;
  std::cout << "This is print out of SFSelection class object" << std::endl; 
  std::cout << "--------------------------------------------------\n" << std::endl; 
  return;
}
//------------------------------------------------------------------
