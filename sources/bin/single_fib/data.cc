// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               data.cc                 *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFDrawCommands.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>

#include <sys/stat.h>
#include <sys/types.h>

const double ampMax = 660;

int main(int argc, char** argv)
{

    gROOT->SetBatch(true);

    TString outdir;
    TString dbase;
    int     seriesNo = -1;

    int ret = parse_common_options(argc, argv, outdir, dbase, seriesNo);
    if (ret != 0) exit(ret);

    if (argc < 2)
    {
        std::cout << "to run type: ./data seriesNo ";
        std::cout << "-out path/to/output -db database" << std::endl;
        return 1;
    }

    SFData* data;

    try
    {
        data = new SFData(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in data.cc!" << std::endl;
        return 1;
    }

    int                 npoints         = data->GetNpoints();
    std::vector<double> positions       = data->GetPositions();
    std::vector<int>    measurementsIDs = data->GetMeasurementsIDs();
    TString             description     = data->GetDescription();
    TString             collimator      = data->GetCollimator();
    TString             fiber           = data->GetFiber();
    TString             testBench       = data->GetTestBench();
    data->Print();
    
    if (description.Contains("Module series"))
    {
        std::cerr << "##### Error in data.cc! This is not regular/monitoring series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << description << std::endl;
        return 1;
    }
    
    int colCh0 = kPink - 8;
    int colCh1 = kAzure - 6;

    TString dirName;
    TString stringCh0;
    TString stringCh1;
    TString string;
    TLatex  textCh0;
    TLatex  textCh1;
    TLatex  text;
    textCh0.SetTextSize(0.033);
    textCh0.SetTextColor(colCh0);
    textCh0.SetTextFont(42);
    textCh0.SetNDC(true);
    textCh1.SetTextSize(0.033);
    textCh1.SetTextColor(colCh1);
    textCh1.SetTextFont(42);
    textCh1.SetNDC(true);
    text.SetTextSize(0.031);
    text.SetTextColor(kGray + 2);
    text.SetTextFont(42);
    text.SetNDC(true);

    TLine line;
    line.SetLineColor(kRed);
    line.SetLineStyle(9);

    std::vector<TPaveStats*> paves;

    TF1* fdiag;

    //----- accessing spectra
    double              s     = SFTools::GetSigmaBL(data->GetSiPM());
    std::vector<double> sigma = {s, s};

    TString cutCh0    = " ";
    TString cutCh1    = " ";
    TString cutCh0Ch1 = " ";
    
    if (testBench == "PMI")
    {
        cutCh0    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh0);
        cutCh1    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh1);
        cutCh0Ch1 = SFDrawCommands::GetCut(SFCutType::kPMICombCh0Ch1);
    }
    else 
    {
        cutCh0    = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigma);
        cutCh1    = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigma);
        cutCh0Ch1 = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigma);
    }
    
    TString cutCh2    = SFDrawCommands::GetCut(SFCutType::kSpecCh2, sigma);
    TString cutCh0A   = SFDrawCommands::GetCut(SFCutType::kSpecCh0A, sigma);
    TString cutCh1A   = SFDrawCommands::GetCut(SFCutType::kSpecCh1A, sigma);
    TString cutBLCh0  = SFDrawCommands::GetCut(SFCutType::kBLCh0);
    TString cutBLCh1  = SFDrawCommands::GetCut(SFCutType::kBLCh1);
    TString cutBLCh2  = SFDrawCommands::GetCut(SFCutType::kBLCh2);

    //----- saving
    TString fname       = Form("data_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in data.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    /*********/ //----- Amplitude spectra -----//
    
    if(testBench != "PMI")
    {
        std::vector<TH1D*> hAmpCh0 = data->GetSpectra(0, SFSelectionType::kAmplitude, cutCh0A);
        std::vector<TH1D*> hAmpCh1 = data->GetSpectra(1, SFSelectionType::kAmplitude, cutCh1A);

        TCanvas* can_ampl = new TCanvas("data_ampl", "data_ampl", 2000, 1200);
        can_ampl->DivideSquare(npoints);

        double amp_min_yaxis = 0.;
        double amp_max_yaxis =  SFTools::FindMaxYaxis(hAmpCh1[0]);
    
        for (int i = 0; i < npoints; i++)
        {
            can_ampl->cd(i + 1);
            gPad->SetGrid(1, 1);
            stringCh0 = hAmpCh0[i]->GetTitle();
            stringCh1 = hAmpCh1[i]->GetTitle();
            //ctx.configureFromJson("hAmpCh0");
            //hAmpCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
            //maxCh0 = hAmpCh0[i]->GetBinContent(hAmpCh0[i]->GetMaximumBin());
            //maxCh1 = hAmpCh1[i]->GetBinContent(hAmpCh1[i]->GetMaximumBin());
            //max_tmp = std::max(maxCh0, maxCh1);
            //hAmpCh0[i]->GetYaxis()->SetRangeUser(0, max_tmp + 0.1 * max_tmp);
            hAmpCh0[i]->GetYaxis()->SetRangeUser(amp_min_yaxis, amp_max_yaxis);
            hAmpCh0[i]->SetTitle(Form("Amplitude spectrum, source position %.2f mm", positions[i]));
            hAmpCh0[i]->GetXaxis()->SetTitle("signal amplitude [mV]");
            hAmpCh0[i]->GetYaxis()->SetTitle("counts");
            hAmpCh0[i]->GetYaxis()->SetMaxDigits(2);
            hAmpCh0[i]->SetStats(false);
            hAmpCh0[i]->SetLineColor(colCh0);
            hAmpCh1[i]->SetLineColor(colCh1);
            hAmpCh1[i]->SetStats(false);
            hAmpCh0[i]->Draw();
            hAmpCh1[i]->Draw("same");
            //line.DrawLine(ampMax, ctx.y.min, ampMax, ctx.y.max);
            line.DrawLine(ampMax, amp_min_yaxis, ampMax, amp_max_yaxis);
            textCh0.DrawLatex(0.3, 0.8, stringCh0);
            textCh1.DrawLatex(0.3, 0.75, stringCh1);
        }
        file->cd();
        can_ampl->Write();
        
        dirName = "AmpSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_1 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hAmpCh0)
            h->Write();
       
        for (auto h : hAmpCh1)
            h->Write();
        
        for (auto h : hAmpCh0)
            delete h;
        
        for (auto h : hAmpCh1)
            delete h;
        
        delete can_ampl;
    }
    
    /*********/ //----- Charge spectra -----//

    std::vector<TH1D*> hChargeCh0;
    std::vector<TH1D*> hChargeCh1; 
    
    if(testBench == "PMI")
    {
        hChargeCh0 = data->GetSpectra(0, SFSelectionType::kPMICharge, cutCh0);
        hChargeCh1 = data->GetSpectra(1, SFSelectionType::kPMICharge, cutCh1);
    }
    else 
    {
        hChargeCh0 = data->GetSpectra(0, SFSelectionType::kPE, cutCh0);
        hChargeCh1 = data->GetSpectra(1, SFSelectionType::kPE, cutCh1);
    }

    TCanvas* can_charge = new TCanvas("data_charge", "data_charge", 2000, 1200);
    can_charge->DivideSquare(npoints);

    double q_min_yaxis = 0.;
    double q_max_yaxis =  SFTools::FindMaxYaxis(hChargeCh1[0]);
    
    double q_min_xaxis = 10.;
    double q_max_xaxis = SFTools::FindMaxXaxis(hChargeCh0[0]);

    for (int i = 0; i < npoints; i++)
    {
        can_charge->cd(i + 1);
        gPad->SetGrid(1, 1);
        stringCh0 = hChargeCh0[i]->GetTitle();
        stringCh1 = hChargeCh1[i]->GetTitle();
        //ctx.configureFromJson("hChargeCh0");
        //hChargeCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //hChargeCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //maxCh0 = hChargeCh0[i]->GetBinContent(hChargeCh0[i]->GetMaximumBin());
        //maxCh1 = hChargeCh1[i]->GetBinContent(hChargeCh1[i]->GetMaximumBin());
        //max_q = std::max(maxCh0, maxCh1);
        //hChargeCh0[i]->GetYaxis()->SetRangeUser(0, max_q + max_q * 0.1);
        hChargeCh0[i]->GetYaxis()->SetRangeUser(q_min_yaxis, q_max_yaxis);
        hChargeCh0[i]->GetXaxis()->SetRangeUser(q_min_xaxis, q_max_xaxis);
        hChargeCh0[i]->SetTitle(Form("Charge spectrum, source position %.2f mm", positions[i]));
        hChargeCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        hChargeCh0[i]->GetYaxis()->SetTitle("counts");
        hChargeCh0[i]->GetYaxis()->SetMaxDigits(2);
        hChargeCh0[i]->SetStats(false);
        hChargeCh0[i]->SetLineColor(colCh0);
        hChargeCh1[i]->SetStats(false);
        hChargeCh1[i]->SetLineColor(colCh1);
        hChargeCh0[i]->Draw();
        hChargeCh1[i]->Draw("same");
        textCh0.DrawLatex(0.3, 0.8, stringCh0);
        textCh1.DrawLatex(0.3, 0.75, stringCh1);
    }
    file->cd();
    can_charge->Write();
    
    dirName = "PESpectra";
    if (!file->cd(dirName)) 
        TDirectory *dir_2 = file->mkdir(dirName);
        
    file->cd(dirName);
        
    for (auto h : hChargeCh0)
        h->Write();
       
    for (auto h : hChargeCh1)
        h->Write();
            
    for (auto h : hChargeCh0)
        delete h;
    
    for (auto h : hChargeCh1)
        delete h;
    
    delete can_charge;

    /*********/ //----- Charge correlation spectra -----//

    std::vector<TH2D*> hCorrPE;
    
    if (testBench == "PMI")
        hCorrPE = data->GetCorrHistograms(SFSelectionType::kPMIChargeCorrelation, cutCh0Ch1);
    else
        hCorrPE = data->GetCorrHistograms(SFSelectionType::kPECorrelation, cutCh0Ch1);
 
    TCanvas* can_charge_corr = new TCanvas("data_charge_corr", "data_charge_corr", 2000, 1200);
    can_charge_corr->DivideSquare(npoints);

    for (int i = 0; i < npoints; i++)
    {
        can_charge_corr->cd(i + 1);
        gPad->SetGrid(1, 1);
        string = hCorrPE[i]->GetTitle();
        hCorrPE[i]->SetTitle(Form("Charge correlation spectrum, source "
                                  "position %.2f mm", positions[i]));
        hCorrPE[i]->GetXaxis()->SetTitle("Ch1 charge [P.E.]");
        hCorrPE[i]->GetYaxis()->SetTitle("Ch0 charge [P.E.]");
        //ctx.configureFromJson("hCorrPE");
        //hCorrPE[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //hCorrPE[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //hCorrPE[i]->GetXaxis()->SetRangeUser(0, max_q + 0.1 * max_q);
        //hCorrPE[i]->GetYaxis()->SetRangeUser(0, max_q + 0.1 * max_q);
        hCorrPE[i]->GetXaxis()->SetRangeUser(0, q_max_xaxis);
        hCorrPE[i]->GetYaxis()->SetRangeUser(0, q_max_xaxis);
        hCorrPE[i]->SetStats(false);
        hCorrPE[i]->Draw("colz");
        fdiag = new TF1("fdiag", "x[0]", q_min_xaxis, q_max_xaxis);
        fdiag->Draw("same");
        text.DrawLatex(0.15, 0.85, string);
     }
     file->cd();
     can_charge_corr->Write();
     
     dirName = "PECorrSpectra";
     if (!file->cd(dirName)) 
         TDirectory *dir_3 = file->mkdir(dirName);
        
     file->cd(dirName);
        
     for (auto h : hCorrPE)
         h->Write();
    
     for (auto h : hCorrPE)
         delete h;
     
     delete can_charge_corr;

    /*********/ //----- T0 spectra -----//

    std::vector<TH1D*> hT0Ch0; 
    std::vector<TH1D*> hT0Ch1; 

    if (testBench == "PMI")
    {
        hT0Ch0 = data->GetSpectra(0, SFSelectionType::kPMIT0, cutCh0);
        hT0Ch1 = data->GetSpectra(1, SFSelectionType::kPMIT0, cutCh1);
    }
    else
    {
        hT0Ch0 = data->GetSpectra(0, SFSelectionType::kT0, cutCh0);
        hT0Ch1 = data->GetSpectra(1, SFSelectionType::kT0, cutCh1);
    }
    
    TCanvas* can_t0 = new TCanvas("data_t0", "data_t0", 2000, 1200);
    can_t0->DivideSquare(npoints);

    double t0_min_yaxis = 0.;
    double t0_max_yaxis =  SFTools::FindMaxYaxis(hT0Ch1[4]);

    for (int i = 0; i < npoints; i++)
    {
        can_t0->cd(i + 1);
        gPad->SetGrid(1, 1);
        stringCh0       = hT0Ch0[i]->GetTitle();
        stringCh1       = hT0Ch0[i]->GetTitle();
        //double maxCh0   = hT0Ch0[i]->GetBinContent(hT0Ch0[i]->GetMaximumBin());
        //double maxCh1   = hT0Ch1[i]->GetBinContent(hT0Ch1[i]->GetMaximumBin());
        //double maxYaxis = std::max(maxCh0, maxCh1);
        //maxYaxis += maxYaxis * 0.1;
        //hT0Ch0[i]->GetYaxis()->SetRangeUser(0, maxYaxis);
        if (testBench != "PMI") 
            hT0Ch0[i]->GetYaxis()->SetRangeUser(t0_min_yaxis, t0_max_yaxis);
        else
            hT0Ch0[i]->GetYaxis()->SetRangeUser(0, 
                     TMath::Max(hT0Ch0[i]->GetBinContent(hT0Ch0[i]->GetMaximumBin()), 
                     hT0Ch1[i]->GetBinContent(hT0Ch1[i]->GetMaximumBin())) + 100);
        hT0Ch0[i]->SetTitle(Form("T_{0} spectrum, source position %.2f mm", positions[i]));
        hT0Ch0[i]->GetXaxis()->SetTitle("time [ns]");
        hT0Ch0[i]->GetYaxis()->SetTitle("counts");
        if (testBench != "PMI") 
            hT0Ch0[i]->GetXaxis()->SetRangeUser(150, 350);
        hT0Ch0[i]->SetLineColor(colCh0);
        hT0Ch1[i]->SetLineColor(colCh1);
        hT0Ch0[i]->Draw();
        gPad->Update();
        paves.push_back((TPaveStats*)hT0Ch0[i]->FindObject("stats"));
        if (paves[i] == nullptr) std::cout << "Warning " << i << std::endl;
        paves[i]->SetY1NDC(0.55);
        paves[i]->SetY2NDC(0.71);
        hT0Ch1[i]->Draw("sames");
        gPad->Update();
        textCh0.DrawLatex(0.2, 0.3, stringCh0);
        textCh1.DrawLatex(0.2, 0.25, stringCh1);
    }
    file->cd();
    can_t0->Write();
     
    dirName = "T0Spectra";
    if (!file->cd(dirName)) 
        TDirectory *dir_4 = file->mkdir(dirName);
        
    file->cd(dirName);
        
    for (auto h : hT0Ch0)
        h->Write();
        
    for (auto h : hT0Ch1)
        h->Write();
    
    for (auto h : hT0Ch0)
         delete h;
    
    for (auto h : hT0Ch1)
         delete h;
    
    delete can_t0;

    /*********/ //----- TOT spectra -----//

    if(testBench != "PMI")
    {
        std::vector<TH1D*> hTOTCh0 = data->GetSpectra(0, SFSelectionType::kTOT, cutCh0);
        std::vector<TH1D*> hTOTCh1 = data->GetSpectra(1, SFSelectionType::kTOT, cutCh1);

        TCanvas* can_tot = new TCanvas("data_tot", "data_tot", 2000, 1200);
        can_tot->DivideSquare(npoints);

        double tot_min_yaxis = 0.;
        double tot_max_yaxis =  SFTools::FindMaxYaxis(hTOTCh0[4]);
        
        double tot_min_xaxis = 10.;
        double tot_max_xaxis = SFTools::FindMaxXaxis(hTOTCh0[4]);

        for (int i = 0; i < npoints; i++)
        {
            can_tot->cd(i + 1);
            gPad->SetGrid(1, 1);
            stringCh0 = hTOTCh0[i]->GetTitle();
            stringCh1 = hTOTCh1[i]->GetTitle();
            //ctx.configureFromJson("hTOTCh0");
            //hTOTCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
            //hTOTCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
            //maxCh0 = hTOTCh0[i]->GetBinContent(hTOTCh0[i]->GetMaximumBin());
            //maxCh1 = hTOTCh1[i]->GetBinContent(hTOTCh1[i]->GetMaximumBin());
            //max_tmp = std::max(maxCh0, maxCh1);
            //hTOTCh0[i]->GetYaxis()->SetRangeUser(0, max_tmp + 0.1 * max_tmp);
            hTOTCh0[i]->GetYaxis()->SetRangeUser(tot_min_yaxis, tot_max_yaxis);
            hTOTCh0[i]->GetXaxis()->SetRangeUser(tot_min_xaxis, tot_max_xaxis);
            hTOTCh0[i]->SetTitle(Form("TOT spectrum, source position %.2f mm", positions[i]));
            hTOTCh0[i]->GetXaxis()->SetTitle("time [ns]");
            hTOTCh0[i]->GetYaxis()->SetTitle("counts");
            hTOTCh0[i]->SetStats(false);
            hTOTCh0[i]->SetLineColor(colCh0);
            hTOTCh1[i]->SetStats(false);
            hTOTCh1[i]->SetLineColor(colCh1);
            hTOTCh0[i]->Draw();
            hTOTCh1[i]->Draw("same");
            textCh0.DrawLatex(0.4, 0.8, stringCh0);
            textCh1.DrawLatex(0.4, 0.75, stringCh1);
        }
        file->cd();
        can_tot->Write();
        
        dirName = "TOTSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_5 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hTOTCh0)
            h->Write();
        
        for (auto h : hTOTCh1)
            h->Write();
        
        for (auto h : hTOTCh0)
            delete h;
        
        for (auto h : hTOTCh1)
            delete h;
        
        delete can_tot;
    }
    
    /*********/ //----- Base line spectrum -----//

    if(testBench != "PMI")
    {
        std::vector<TH1D*> hBLCh0 = data->GetSpectra(0, SFSelectionType::kBL, cutBLCh0);
        std::vector<TH1D*> hBLCh1 = data->GetSpectra(1, SFSelectionType::kBL, cutBLCh1);
        std::vector<TH1D*> hBLCh2 = data->GetSpectra(2, SFSelectionType::kBL, cutBLCh2);

        TCanvas* can_bl_ch0 = new TCanvas("data_bl_ch0", "can_bl_ch0", 2000, 1200);
        can_bl_ch0->Divide(npoints);
        
        TCanvas* can_bl_ch1 = new TCanvas("data_bl_ch1", "can_bl_ch1", 2000, 1200);
        can_bl_ch1->Divide(npoints);
        
        TCanvas* can_bl_ch2 = new TCanvas("data_bl_ch2", "can_bl_ch2", 2000, 1200);
        can_bl_ch2->Divide(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_bl_ch0->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLCh0[i]->Draw();
            hBLCh0[i]->GetXaxis()->SetName("baseline [ADC channels]");
            hBLCh0[i]->GetYaxis()->SetName("counts");
            can_bl_ch1->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLCh1[i]->Draw();
            hBLCh1[i]->GetXaxis()->SetName("baseline [ADC channels]");
            hBLCh1[i]->GetYaxis()->SetName("counts");
            can_bl_ch2->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLCh2[i]->Draw();
            hBLCh2[i]->GetXaxis()->SetName("baseline [ADC channels]");
            hBLCh2[i]->GetYaxis()->SetName("counts");
        }
        file->cd();
        can_bl_ch0->Write();
        can_bl_ch1->Write();
        can_bl_ch2->Write();
        
        dirName = "BLSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_6 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hBLCh0)
            delete h;
        
        for (auto h : hBLCh1)
            h->Write();
        
        for (auto h : hBLCh2)
            h->Write();
        
        for (auto h : hBLCh0)
            delete h;
        
        for (auto h : hBLCh1)
            delete h;
        
        for (auto h : hBLCh2)
            delete h;
        
        delete can_bl_ch0;
        delete can_bl_ch1;
        delete can_bl_ch2;
    }

    /*********/ //----- Base line sigma spectrum -----//

    if(testBench != "PMI")
    {
        std::vector<TH1D*> hBLSigmaCh0 = data->GetSpectra(0, SFSelectionType::kBLSigma, cutBLCh0);
        std::vector<TH1D*> hBLSigmaCh1 = data->GetSpectra(1, SFSelectionType::kBLSigma, cutBLCh1);
        std::vector<TH1D*> hBLSigmaCh2 = data->GetSpectra(2, SFSelectionType::kBLSigma, cutBLCh2);

        TCanvas* can_bls_ch0 = new TCanvas("data_bls_ch0", "data_bls_ch0", 2000, 1200);
        can_bls_ch0->Divide(npoints);

        TCanvas* can_bls_ch1 = new TCanvas("data_bls_ch1", "data_bls_ch1", 2000, 1200);
        can_bls_ch1->Divide(npoints);

        TCanvas* can_bls_ch2 = new TCanvas("data_bls_ch2", "data_bls_ch2", 2000, 1200);
        can_bls_ch2->Divide(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_bls_ch0->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLSigmaCh0[i]->Draw();
            hBLSigmaCh0[i]->GetXaxis()->SetTitle("baseline sigma [ADC channels]");
            hBLSigmaCh0[i]->GetYaxis()->SetTitle("counts");
            can_bls_ch1->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLSigmaCh1[i]->Draw();
            hBLSigmaCh1[i]->GetXaxis()->SetTitle("baseline sigma [ADC channels]");
            hBLSigmaCh1[i]->GetYaxis()->SetTitle("counts");
            can_bls_ch2->cd(i + 1);
            gPad->SetGrid(1, 1);
            hBLSigmaCh2[i]->Draw();
            hBLSigmaCh2[i]->GetXaxis()->SetTitle("baseline sigma [ADC channels]");
            hBLSigmaCh2[i]->GetYaxis()->SetTitle("counts");
        }
        file->cd();
        can_bls_ch0->Write();
        can_bls_ch1->Write();
        can_bls_ch2->Write();
        
        dirName = "BLSigSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_7 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hBLSigmaCh0)
            h->Write();
        
        for (auto h : hBLSigmaCh1)
            h->Write();
        
        for (auto h : hBLSigmaCh2)
            h->Write();
        
        for (auto h : hBLSigmaCh0)
            delete h;
        
        for (auto h : hBLSigmaCh1)
            delete h;
        
        for (auto h : hBLSigmaCh2)
            delete h;
        
        delete can_bls_ch0;
        delete can_bls_ch1;
        delete can_bls_ch2;
    }

    /*********/ //----- Amplitude correlation spectra -----//

    if(testBench != "PMI")
    {
        std::vector<TH2D*> hCorrAmp = data->GetCorrHistograms(SFSelectionType::kAmplitudeCorrelation, cutCh0Ch1);

        TCanvas* can_ampl_corr = new TCanvas("data_ampl_corr", "data_ampl_corr", 2000, 1200);
        can_ampl_corr->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_ampl_corr->cd(i + 1);
            gPad->SetGrid(1, 1);
            string = hCorrAmp[i]->GetTitle();
            hCorrAmp[i]->SetTitle(Form("Amplitude correlation spectrum, source "
                                "position %.2f mm", positions[i]));
            hCorrAmp[i]->GetXaxis()->SetTitle("Ch1 amplitude [mV]");
            hCorrAmp[i]->GetYaxis()->SetTitle("Ch0 amplitude [mV]");
            hCorrAmp[i]->SetStats(false);
            hCorrAmp[i]->Draw("colz");
            text.DrawLatex(0.15, 0.85, string);
        }
        file->cd();
        can_ampl_corr->Write();
        
        dirName = "AmpCorrSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_8 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hCorrAmp)
            h->Write();
        
        for (auto h : hCorrAmp)
            delete h;
        
        delete can_ampl_corr;
    }
    
    /*********/ //----- T0 correlation spectra -----//

    std::vector<TH2D*> hCorrT0;
    
    if (testBench == "PMI")
        hCorrT0 = data->GetCorrHistograms(SFSelectionType::kPMIT0Correlation, cutCh0Ch1);
    else
        hCorrT0 = data->GetCorrHistograms(SFSelectionType::kT0Correlation, cutCh0Ch1);

    TCanvas* can_t0_corr = new TCanvas("data_t0_corr", "data_t0_corr", 2000, 1200);
    can_t0_corr->DivideSquare(npoints);

    for (int i = 0; i < npoints; i++)
    {
        can_t0_corr->cd(i + 1);
        gPad->SetGrid(1, 1);
        string = hCorrT0[i]->GetTitle();
        hCorrT0[i]->SetTitle(Form("T0 correlation spectrum, source "
                             "position %.2f mm", positions[i]));
        hCorrT0[i]->GetXaxis()->SetTitle("Ch1 T0 [ns]");
        hCorrT0[i]->GetYaxis()->SetTitle("Ch0 T0 [ns]");
        if (testBench != "PMI")
        {
            hCorrT0[i]->GetXaxis()->SetRangeUser(0, 400);
            hCorrT0[i]->GetYaxis()->SetRangeUser(0, 400);
        }
        else
        {
            hCorrT0[i]->GetXaxis()->SetRangeUser(0, 20);
            hCorrT0[i]->GetYaxis()->SetRangeUser(0, 20);
        }
        hCorrT0[i]->SetStats(false);
        hCorrT0[i]->Draw("colz");
        text.DrawLatex(0.15, 0.85, string);
    }
    file->cd();
    can_t0_corr->Write();
     
    dirName = "T0CorrSpectra";
    if (!file->cd(dirName)) 
        TDirectory *dir_9 = file->mkdir(dirName);
        
    file->cd(dirName);
        
    for (auto h : hCorrT0)
        h->Write();
       
    for (auto h : hCorrT0)
        delete h;
    
    delete can_t0_corr;
    
    /*********/ //----- Amplitude vs charge correlation spectra -----//

    if(testBench != "PMI")
    {
        std::vector<TH2D*> hAmpPECh0 = data->GetCorrHistograms(SFSelectionType::kAmpPECorrelation, cutCh0, 0);

        TCanvas* can_amp_pe_ch0 = new TCanvas("data_amp_pe_ch0", "data_amp_pe_ch0", 2000, 1200);
        can_amp_pe_ch0->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_amp_pe_ch0->cd(i + 1);
            gPad->SetGrid(1, 1);
            string = hAmpPECh0[i]->GetTitle();
            hAmpPECh0[i]->SetTitle(Form("Amplitude vs. Charge correlation spectrum Ch0, source "
                                        "position %.2f mm", positions[i]));
            hAmpPECh0[i]->GetXaxis()->SetTitle("Charge [PE]");
            hAmpPECh0[i]->GetYaxis()->SetTitle("Amplitude [mV]");
            //hAmpPECh0[i]->GetXaxis()->SetRangeUser(-10, ctx.x.max);
            //hAmpPECh0[i]->GetXaxis()->SetRangeUser(-10, max_q + 0.1 * max_q);
            hAmpPECh0[i]->GetXaxis()->SetRangeUser(0, q_max_xaxis);
            hAmpPECh0[i]->SetStats(false);
            hAmpPECh0[i]->Draw("colz");
            text.DrawLatex(0.15, 0.85, string);
        }
        file->cd();
        can_amp_pe_ch0->Write();
        
        dirName = "AmpPECh0CorrSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_10 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hAmpPECh0)
            h->Write();
        
        for (auto h : hAmpPECh0)
            delete h;
        
        delete can_amp_pe_ch0;
        
        /*********/

        std::vector<TH2D*> hAmpPECh1 =
            data->GetCorrHistograms(SFSelectionType::kAmpPECorrelation, cutCh1, 1);

        TCanvas* can_amp_pe_ch1 = new TCanvas("data_amp_pe_ch1", "data_amp_pe_ch1", 2000, 1200);
        can_amp_pe_ch1->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_amp_pe_ch1->cd(i + 1);
            gPad->SetGrid(1, 1);
            string = hAmpPECh1[i]->GetTitle();
            hAmpPECh1[i]->SetTitle(Form("Amplitude vs. Charge correlation spectrum Ch1, source "
                                        "position %.2f mm", positions[i]));
            hAmpPECh1[i]->GetXaxis()->SetTitle("Charge [PE]");
            hAmpPECh1[i]->GetYaxis()->SetTitle("Amplitude [mV]");
            //hAmpPECh1[i]->GetXaxis()->SetRangeUser(-10, ctx.x.max);
            //hAmpPECh1[i]->GetXaxis()->SetRangeUser(-10, max_q + 0.1 * max_q);
            hAmpPECh1[i]->GetXaxis()->SetRangeUser(0, q_max_xaxis);
            hAmpPECh1[i]->GetYaxis()->SetRangeUser(-10, 800);
            hAmpPECh1[i]->SetStats(false);
            hAmpPECh1[i]->Draw("colz");
            text.DrawLatex(0.15, 0.85, string);
        }
        file->cd();
        can_amp_pe_ch1->Write();
        
        dirName = "AmpPECh1CorrSpectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_11 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hAmpPECh1)
            h->Write();
        
        for (auto h : hAmpPECh1)
            delete h;
        
        delete can_amp_pe_ch1;
    }

    /*********/ //----- Reference channel spectra -----//
    
    if (collimator.Contains("Electronic"))
    {
        /*********/

        std::vector<TH1D*> hChargeCh2 = data->GetSpectra(2, SFSelectionType::kCharge, cutCh2);
        TCanvas*           can_ref    = new TCanvas("data_ref", "data_ref", 2000, 1200);
        can_ref->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {

            can_ref->cd(i + 1);
            gPad->SetGrid(1, 1);
            hChargeCh2[i]->SetStats(false);
            string = hChargeCh2[i]->GetTitle();
            hChargeCh2[i]->GetXaxis()->SetTitle("charge [a.u.]");
            hChargeCh2[i]->GetYaxis()->SetTitle("counts");
            hChargeCh2[i]->GetXaxis()->SetRangeUser(0, 120E3);
            hChargeCh2[i]->SetTitle(Form("Charge spectrum, reference detector, "
                                         "position %.2f mm", positions[i]));
            hChargeCh2[i]->Draw();
            text.DrawLatex(0.3, 0.8, string);
        }
        file->cd();
        can_ref->Write();
        
        dirName = "PECh2Spectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_12 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hChargeCh2)
            h->Write();
        
        for (auto h : hChargeCh2)
            delete h;
        
        delete can_ref;

        /*********/

        std::vector<TH2D*> hChargeCh0Ch2 = data->GetRefCorrHistograms(0);
        TCanvas*           can_ref_ch0   = new TCanvas("data_ref_ch0", "data_ref_ch0", 2000, 1200);
        can_ref_ch0->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {

            can_ref_ch0->cd(i + 1);
            gPad->SetGrid(1, 1);
            string = hChargeCh0Ch2[i]->GetTitle();
            hChargeCh0Ch2[i]->SetTitle(Form("Charge correlation spectrum Ch2 vs. Ch0, source "
                                            "position %.2f mm", positions[i]));
            hChargeCh0Ch2[i]->GetXaxis()->SetTitle("Ch2 charge [PE]");
            hChargeCh0Ch2[i]->GetYaxis()->SetTitle("Ch0 charge [a.u.]");
            //ctx.configureFromJson("hChargeChXCh2");
            hChargeCh0Ch2[i]->GetXaxis()->SetRangeUser(0, 120E3);
            //hChargeCh0Ch2[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
            hChargeCh0Ch2[i]->GetYaxis()->SetRangeUser(0, q_max_xaxis);
            hChargeCh0Ch2[i]->SetStats(false);
            hChargeCh0Ch2[i]->Draw("colz");
            text.DrawLatex(0.15, 0.85, string);
        }
        file->cd();
        can_ref_ch0->Write();
        
        dirName = "PECh0PECh2Spectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_13 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hChargeCh0Ch2)
            h->Write();
        
        for (auto h : hChargeCh0Ch2)
            delete h;
        
        delete can_ref_ch0;

        /*********/

        std::vector<TH2D*> hChargeCh1Ch2 = data->GetRefCorrHistograms(1);
        TCanvas*           can_ref_ch1   = new TCanvas("data_ref_ch1", "data_ref_ch1", 2000, 1200);
        can_ref_ch1->DivideSquare(npoints);

        for (int i = 0; i < npoints; i++)
        {
            can_ref_ch1->cd(i + 1);
            gPad->SetGrid(1, 1);
            string = hChargeCh1Ch2[i]->GetTitle();
            hChargeCh1Ch2[i]->SetTitle(Form("Charge correlation spectrum Ch2 vs. Ch1, source "
                                            "position %.2f mm", positions[i]));
            hChargeCh1Ch2[i]->GetXaxis()->SetTitle("Ch2 charge [PE]");
            hChargeCh1Ch2[i]->GetYaxis()->SetTitle("Ch1 charge [a.u.]");
            hChargeCh1Ch2[i]->GetXaxis()->SetRangeUser(0, 120E3);
            //hChargeCh1Ch2[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
            hChargeCh1Ch2[i]->GetYaxis()->SetRangeUser(0, q_max_xaxis);
            hChargeCh1Ch2[i]->SetStats(false);
            hChargeCh1Ch2[i]->Draw("colz");
            text.DrawLatex(0.15, 0.85, string);
        }
        file->cd();
        can_ref_ch1->Write();
        
        dirName = "PECh1PECh2Spectra";
        if (!file->cd(dirName)) 
            TDirectory *dir_14 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hChargeCh1Ch2)
            h->Write();
        
        for (auto h : hChargeCh1Ch2)
            delete h;
        
        delete can_ref_ch1;
    }
    
    /*********/ //----- Signals -----//
    
    if(testBench != "PMI")
    {
        const int          nsig   = 6;
        int                number = 0;
        std::vector<TH1D*> hSigCh0(nsig);
        std::vector<TH1D*> hSigCh1(nsig);

        for (int i = 0; i < nsig / 2; i++)
        {
            number     = 100 * (i + 1);
            hSigCh0[i] = data->GetSignal(0, measurementsIDs[2], "", number, true);
            hSigCh0[i + (nsig / 2)] =
                data->GetSignal(0, measurementsIDs[npoints - 1], "", number, true);
            hSigCh1[i] = data->GetSignal(1, measurementsIDs[2], "", number, true);
            hSigCh1[i + (nsig / 2)] =
                data->GetSignal(1, measurementsIDs[npoints - 1], "", number, true);
        }

        TCanvas* can_sig = new TCanvas("data_sig", "data_sig", 1800, 800);
        can_sig->Divide(3, 2);

        int    polarity = 0;
        double min      = hSigCh0[0]->GetBinContent(hSigCh0[0]->GetMinimumBin());
        double max      = hSigCh0[0]->GetBinContent(hSigCh0[0]->GetMaximumBin());

        if (fabs(min) > max)
        {
            polarity = -1;
            std::cout << "min = " << min << "\t max = " << max << std::endl;
            std::cout << "Signals are negative, am I right?" << std::endl;
        }
        else
        {
            polarity = 1;
            std::cout << "min = " << min << "\t max = " << max << std::endl;
            std::cout << "Signals are positive, am I right?" << std::endl;
        }

        for (int i = 0; i < nsig; i++)
        {
            can_sig->cd(i + 1);
            gPad->SetGrid(1, 1);

            stringCh0 = hSigCh0[i]->GetTitle();
            stringCh1 = hSigCh1[i]->GetTitle();
            hSigCh0[i]->SetLineColor(colCh0);
            hSigCh0[i]->SetTitle(" ");
            hSigCh0[i]->GetXaxis()->SetTitle("time [ns]");
            hSigCh0[i]->GetYaxis()->SetTitle("amplitude [mV]");
            hSigCh0[i]->SetStats(false);
            hSigCh1[i]->SetLineColor(colCh1);
            hSigCh1[i]->SetStats(false);
            hSigCh0[i]->Draw();
            hSigCh1[i]->Draw("same");
            if (polarity == 1)
            {
                double maxCh0   = hSigCh0[i]->GetBinContent(hSigCh0[i]->GetMaximumBin());
                double maxCh1   = hSigCh1[i]->GetBinContent(hSigCh1[i]->GetMaximumBin());
                double maxYaxis = std::max(maxCh0, maxCh1) + 10.;
                hSigCh0[i]->GetYaxis()->SetRangeUser(-2, maxYaxis);
            }
            else if (polarity == -1)
            {
                double minCh0   = hSigCh0[i]->GetBinContent(hSigCh0[i]->GetMinimumBin());
                double minCh1   = hSigCh1[i]->GetBinContent(hSigCh1[i]->GetMinimumBin());
                double minYaxis = std::min(minCh0, minCh1) - 10.;
                hSigCh0[i]->GetYaxis()->SetRangeUser(minYaxis, 20);
            }
            textCh0.DrawLatex(0.5, 0.6, stringCh0);
            textCh1.DrawLatex(0.5, 0.55, stringCh1);
        }
        
        file->cd();
        can_sig->Write();
        
        dirName = "Signals";
        if (!file->cd(dirName)) 
            TDirectory *dir_15 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hSigCh0)
            h->Write();
        
        for (auto h : hSigCh1)
            h->Write();
        
        for (auto h : hSigCh0)
            delete h;
        
        for (auto h : hSigCh1)
            delete h;
        
        delete can_sig;
        
        /*********/
        
        const int              nsigav = 3;
        std::vector<TProfile*> hSigAvCh0(nsigav);
        std::vector<TProfile*> hSigAvCh1(nsigav);

        double PE[3];
        if (collimator.Contains("Electronic"))
        {
            if (fiber.Contains("LYSO") || fiber.Contains("GAGG"))
            {
                PE[0] = 150.;
                PE[1] = 200.;
                PE[2] = 280.;
            }
            else if (fiber.Contains("LuAG"))
            {
                PE[0] = 50.;
                PE[1] = 100.;
                PE[2] = 150.;
            }
            else
            {
                std::cerr << "##### Error in data.cc! Unknown fiber material!" << std::endl;
                std::abort();
            }
        }
        else if (collimator.Contains("Lead"))
        {
            PE[0] = 60.;
            PE[1] = 20.;
            PE[2] = 400.;
        }

        //int ID = SFTools::GetMeasurementID(seriesNo, 45.0);
        //int ID = SFTools::GetMeasurementID(seriesNo, 12.0);
        int ID = 0;
        
        if(seriesNo == 194 || seriesNo == 196)
            ID = SFTools::GetMeasurementID(seriesNo, 45.0);
        else 
            ID = SFTools::GetMeasurementID(seriesNo, 50.0);

        hSigAvCh0[0] = data->GetSignalAverage(
            0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", PE[0] - 0.5, PE[0] + 0.5), 20, true);
        hSigAvCh0[1] = data->GetSignalAverage(
            0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", PE[1] - 0.5, PE[1] + 0.5), 20, true);
        hSigAvCh0[2] = data->GetSignalAverage(
            0, ID, Form("ch_0.fPE>%f && ch_0.fPE<%f", PE[2] - 0.5, PE[2] + 0.5), 20, true);

        hSigAvCh1[0] = data->GetSignalAverage(
            1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", PE[0] - 0.5, PE[0] + 0.5), 20, true);
        hSigAvCh1[1] = data->GetSignalAverage(
            1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", PE[1] - 0.5, PE[1] + 0.5), 20, true);
        hSigAvCh1[2] = data->GetSignalAverage(
            1, ID, Form("ch_1.fPE>%f && ch_1.fPE<%f", PE[2] - 0.5, PE[2] + 0.5), 20, true);
        
        TCanvas* can_sigav = new TCanvas("data_sigav", "data_sigav", 1800, 800);
        can_sigav->Divide(3, 2);

        textCh0.SetTextColor(kGray + 2);
        textCh0.SetTextSize(0.025);
        textCh1.SetTextColor(kGray + 2);
        textCh1.SetTextSize(0.025);

        for (int i = 0; i < nsigav; i++)
        {
            can_sigav->cd(i + 1);
            gPad->SetGrid(1, 1);
            stringCh0 = hSigAvCh0[i]->GetTitle();
            hSigAvCh0[i]->SetTitle(" ");
            hSigAvCh0[i]->GetXaxis()->SetTitle("time [ns]");
            hSigAvCh0[i]->GetYaxis()->SetTitle("amplitude [mV]");
            hSigAvCh0[i]->SetStats(false);
            hSigAvCh0[i]->Draw();
            if (polarity == 1)
            {
                double maxYaxis = hSigAvCh0[i]->GetBinContent(hSigAvCh0[i]->GetMaximumBin());
                maxYaxis        = maxYaxis + 0.2 * maxYaxis;
                hSigAvCh0[i]->GetYaxis()->SetRangeUser(-10, maxYaxis);
            }
            else if (polarity == -1)
            {
                double minYaxis = hSigAvCh0[i]->GetBinContent(hSigAvCh0[i]->GetMinimumBin());
                minYaxis        = minYaxis + 0.2 * minYaxis;
                hSigAvCh0[i]->GetYaxis()->SetRangeUser(minYaxis, 100);
            }
            textCh0.DrawLatex(0.15, 0.20, stringCh0);

            can_sigav->cd(i + 1 + nsigav);
            gPad->SetGrid(1, 1);
            stringCh1 = hSigAvCh1[i]->GetTitle();
            hSigAvCh1[i]->SetTitle(" ");
            hSigAvCh1[i]->GetXaxis()->SetTitle("time [ns]");
            hSigAvCh1[i]->GetYaxis()->SetTitle("amplitude [mV]");
            if (polarity == 1)
            {
                double maxYaxis = hSigAvCh1[i]->GetBinContent(hSigAvCh1[i]->GetMaximumBin());
                maxYaxis        = maxYaxis + 0.2 * maxYaxis;
                hSigAvCh1[i]->GetYaxis()->SetRangeUser(-10, maxYaxis);
            }
            else if (polarity == -1)
            {
                double minYaxis = hSigAvCh1[i]->GetBinContent(hSigAvCh1[i]->GetMinimumBin());
                minYaxis        = minYaxis + 0.2 * minYaxis;
                hSigAvCh1[i]->GetYaxis()->SetRangeUser(minYaxis, 100);
            }
            hSigAvCh1[i]->SetStats(false);
            hSigAvCh1[i]->Draw();
            textCh1.DrawLatex(0.15, 0.20, stringCh1);
        }

        file->cd();
        can_sigav->Write();
        
        dirName = "SignalsAve";
        if (!file->cd(dirName)) 
            TDirectory *dir_16 = file->mkdir(dirName);
        
        file->cd(dirName);
        
        for (auto h : hSigAvCh0)
            h->Write();
        
        for (auto h : hSigAvCh1)
            h->Write();
        
        for (auto h : hSigAvCh0)
            delete h;
        
        for (auto h : hSigAvCh1)
            delete h;
        
        delete can_sigav;
    }

    file->Close();

    //----- writing results to the data base
    TString table = "DATA";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES (%i, '%s')",
                         table.Data(), seriesNo, fname_full.Data());

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- data writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- data writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;

    return 0;
}
