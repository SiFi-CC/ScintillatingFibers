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
        default:
            std::cerr << "##### Error in SFDrawCommands::GetSelectionName()!" << std::endl;
            std::cerr << "Unknown selection type! Please check!" << std::endl;
            break;
    }

    CheckSelection(selectionName);
    return selectionName;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object. This
/// method provides selections for simple histograms: Charge, PE, Amplitude
/// T0 and TOT.
/// \param selection - selection type
/// \param unique - unique histogram ID
/// \param ch - channel number
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique, int ch){
    
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
      default:
          std::cerr << "##### Error in SFDrawCommands::GetSelection()!" << std::endl;
          std::cerr << "Unknown selection type! Please check!" << std::endl;
          break;
  }
    
  CheckSelection(selectionString);
  return selectionString;
}
//------------------------------------------------------------------
/// Returns selection as a TString for ROOT's TTree type object. This 
/// method provides selections for the following histograms: LogSqrtPERatio,
/// T0Difference, PEAverage, AmplitudeAverage, PECorrelation, AmplitudeCorrelation,
/// T0Correlation.
/// \param selection - selection type
/// \param unique - unique histogram ID
TString SFDrawCommands::GetSelection(SFSelectionType selection, int unique){
 
  TString selectionString = "";
  
  switch(selection){
      case SFSelectionType::LogSqrtPERatio:
          selectionString = Form("log(sqrt(ch_1.fPE/ch_0.fPE))>>htemp%i(500,-5,5)", unique);
          break;
      case SFSelectionType::T0Difference:
          selectionString = Form("(ch_0.fT0-ch_1.fT0)>>htemp%i(500,-50,50)", unique);
          break;
      case SFSelectionType::PEAverage:
          selectionString = Form("sqrt(ch_0.fPE*ch_1.fPE)>>htemp%i(1000,-150,1200)", unique);
          break;
      case SFSelectionType::AmplitudeAverage:
          selectionString = Form("sqrt(ch_0.fAmp*ch_1.fAmp)>>htemp%i(1000,0,700)",unique);
          break;
      case SFSelectionType::PECorrelation:
          selectionString = Form("ch_0.fPE:ch_1.fPE>>htemp%i(1000,-150,1200,1000,-150,1200)", unique);
          break;
      case SFSelectionType::AmplitudeCorrelation:
          selectionString = Form("ch_0.fAmp:ch_1.fAmp>>htemp%i(1000,0,700,1000,0,700)", unique);
          break;
      case SFSelectionType::T0Correlation:
          selectionString = Form("ch_0.fT0:ch_1.fT0>>htemp%i(1000,-110,1100,1000,-110,1100", unique);
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
/// Returns cut as a TString for ROOT's TTree type object. This method provides 
/// fixed cuts, where no additional parameters are needed.
/// \param cut - cut type.
TString SFDrawCommands::GetCut(SFCutType cut){
    
  TString cutString = "";
   
  switch(cut){
      case SFCutType::PEAbove0:
          cutString = "ch_0.fPE>0 && ch_1.fPE>0";
          break;
      case SFCutType::T0Above0:
          cutString = "ch_0.fT0>0 && ch_1.fT0>0";
          break;
      case SFCutType::T0Valid:
          cutString = "ch_0.fT0>-100 && ch_1.fT0>-100";
          break;
      case SFCutType::TOTValid:
          cutString = "ch_0.fTOT>-100 && ch_1.fTOT>-100";
          break;
      case SFCutType::Empty:
          cutString = "";
          break;
      default:
          std::cerr << "##### Error in SFDrawCommands::GetCut()!" << std::endl;
          std::cerr << "Unknown cut type! Please check!" << std::endl;
          break;
  }
  
  return cutString;
}
//------------------------------------------------------------------
/// Returns cut as a TString for ROOT's TTree type object. This method provides 
/// custom cuts for chosen channel.
/// \param cut - cut type
/// \param customNumbers - standard vector containing values to be inserted in the cut
/// \param ch - channel number.
TString SFDrawCommands::GetCut(SFCutType cut, std::vector <double> customNumbers, int ch){
    
  TString cutString = "";
    
  switch(cut){
      case SFCutType::T0ChBelowMax:
        cutString = Form("ch_0.fT0<%f", ch, customNumbers[0]);
        break;
      case SFCutType::Only511:
        cutString = Form("ch_%i.fPE>%f && ch_%i.fPE<%f", ch, customNumbers[0], ch, customNumbers[1]);
        break;
      default:
        std::cerr << "##### Error in SFDrawCommands::GetCut()!" << std::endl;
        std::cerr << "Unknown cut type! Please check!" << std::endl;
        break;
  }
    
  return cutString;
}
//------------------------------------------------------------------
/// Returns cut as a TString for ROOT's TTree type object. This method provides
/// custom cut for combination of channels 0 and 1.
/// \param cut - cut type
/// \param customNumbers - standard vector containing values to be inserted in the cut.
TString SFDrawCommands::GetCut(SFCutType cut, std::vector <double> customNumbers){
 
  TString cutString = "";
    
  switch(cut){
      case SFCutType::T0BelowMax:
        cutString = Form("ch_0.fT0<%f && ch_1.fT0<%f", customNumbers[0], customNumbers[1]);
        break;
      case SFCutType::ScatteredEvents:
          cutString = Form("log(sqrt(ch_1.fPE/ch_0.fPE))>%f && log(sqrt(ch_1.fPE/ch_0.fPE))<%f",
                           customNumbers[0], customNumbers[1]);
          break;
      default:
        std::cerr << "##### Error in SFDrawCommands::GetCut!" << std::endl;
        std::cerr << "Unknown cut type! Please check!" << std::endl;
        break;
  }
  
  return cutString;
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
