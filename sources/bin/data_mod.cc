// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             data_mod.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFDrawCommands.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TROOT.h>

#include <sys/stat.h>
#include <sys/types.h>

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
        std::cout << "to run type: ./data_mod seriesNo ";
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
        std::cerr << "##### Exception in data_mod.cc!" << std::endl;
        return 1;
    }
    
    int                 npoints         = data->GetNpoints();
    std::vector<double> fibers          = data->GetPositions();
    TString             description     = data->GetDescription();
    data->Print();
    
    if (! description.Contains("Module series"))
    {
        std::cerr << "##### Error in data.cc! This is not module series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << description << std::endl;
        return 1;
    }
    
    const int ncans   = 4;
    int       counter = 0;
    
    int colCh0 = kPink - 8;
    int colCh1 = kAzure - 6;
    
    TString cutCh0    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh0);
    TString cutCh1    = SFDrawCommands::GetCut(SFCutType::kPMISpecCh1);
    TString cutCh0Ch1 = SFDrawCommands::GetCut(SFCutType::kPMICombCh0Ch1);

    TString dirName   = "";
    TString stringCh0 = "";
    TString stringCh1 = "";
    TString string    = "";
    
    TLatex textCh0;
    textCh0.SetTextSize(0.033);
    textCh0.SetTextColor(colCh0);
    textCh0.SetTextFont(42);
    textCh0.SetNDC(true);
    
    TLatex textCh1;
    textCh1.SetTextSize(0.033);
    textCh1.SetTextColor(colCh1);
    textCh1.SetTextFont(42);
    textCh1.SetNDC(true);
    
    TLatex text;
    textCh1.SetTextSize(0.033);
    textCh1.SetTextColor(kGray + 2);
    textCh1.SetTextFont(42);
    textCh1.SetNDC(true);
    
    //----- creating/opening results file 
    TString fname       = Form("data_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in data_mod.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }
    
    /*********/ //----- Photon count spectra -----//

    std::vector<TH1D*> hPhCountCh0 = data->GetSpectra(0, SFSelectionType::kPMICharge, cutCh0);
    std::vector<TH1D*> hPhCountCh1 = data->GetSpectra(1, SFSelectionType::kPMICharge, cutCh1);
    
    std::vector<TCanvas*> can_phcnt(ncans);
    
    double phcnt_min_yaxis = 0;
    double phcnt_max_yaxis = SFTools::FindMaxYaxis(hPhCountCh1[0]);
    
    double phcnt_min_xaxis = 0;
    double phcnt_max_xaxis = SFTools::FindMaxXaxis(hPhCountCh0[0]);
    
    for (int ican = 0; ican < ncans; ican++)
    {
        TString cname = Form("data_phcnt_pt%i", ican+1);
        can_phcnt[ican] = new TCanvas(cname, cname, 2000, 1200);
        can_phcnt[ican]->DivideSquare(npoints/ncans);
        
        for (int i = 0; i < (npoints/ncans); i++)
        {
            can_phcnt[ican]->cd(i+1);
            gPad->SetGrid(1, 1);
            
            stringCh0 = hPhCountCh0[counter]->GetTitle();
            stringCh1 = hPhCountCh1[counter]->GetTitle();
            
            hPhCountCh0[counter]->SetTitle(Form("Photon count spectrum, fiber %.0f", fibers[counter]));
            hPhCountCh0[counter]->GetXaxis()->SetRangeUser(phcnt_min_xaxis, phcnt_max_xaxis);
            hPhCountCh0[counter]->GetYaxis()->SetRangeUser(phcnt_min_yaxis, phcnt_max_yaxis);
            hPhCountCh0[counter]->GetXaxis()->SetTitle("photon count");
            hPhCountCh0[counter]->GetYaxis()->SetTitle("counts");
            hPhCountCh0[counter]->GetYaxis()->SetMaxDigits(2);
            hPhCountCh0[counter]->SetLineColor(colCh0);
            hPhCountCh0[counter]->SetStats(false);
            hPhCountCh0[counter]->Draw();
            
            hPhCountCh1[counter]->SetStats(false);
            hPhCountCh1[counter]->SetLineColor(colCh1);
            hPhCountCh1[counter]->Draw("same");
            
            textCh0.DrawLatex(0.3, 0.8, stringCh0);
            textCh1.DrawLatex(0.3, 0.7, stringCh1);
            
            counter++;
        }
    }
    
    //----- saving 
    file->cd();
    
    for (auto c : can_phcnt)
        c->Write();
    
    dirName = "PhCountSpectra";
    if (!file->cd(dirName))
        TDirectory *dir_1 = file->mkdir(dirName);
    
    file->cd(dirName);
    
    for (auto h : hPhCountCh0)
        h->Write();
    
    for (auto h : hPhCountCh1)
        h->Write();
    
    //----- deleting
    for (auto c : can_phcnt)
        delete c;
    
    for (auto h : hPhCountCh0)
        delete h;
    
    for (auto h : hPhCountCh1)
        delete h;
    
    /*********/ //----- Photon count correlation spectra -----//
    
    TF1* fdiag;
    counter = 0;
    
    std::vector<TH2D*> hCorrPhCount = data->GetCorrHistograms(SFSelectionType::kPMIChargeCorrelation, cutCh0Ch1);
    
    std::vector<TCanvas*> can_phcnt_corr(ncans);
    
    for (int ican = 0; ican < ncans; ican++)
    {
        TString cname = Form("data_phcnt_corr_pt%i", ican+1);
        can_phcnt_corr[ican] = new TCanvas(cname, cname, 2000, 1200);
        can_phcnt_corr[ican]->DivideSquare(npoints/ncans);
        
        for (int i = 0; i < (npoints/ncans); i++)
        {
            can_phcnt_corr[ican]->cd(counter+1);
            gPad->SetGrid(1, 1);
            
            string = hCorrPhCount[counter]->GetTitle();
            
            hCorrPhCount[counter]->SetTitle(Form("Photon count correlation spectrum, fiber %.0f", fibers[counter]));
            hCorrPhCount[counter]->GetXaxis()->SetTitle("Ch1 photon count");
            hCorrPhCount[counter]->GetYaxis()->SetTitle("Ch0 photon count");
            hCorrPhCount[counter]->GetXaxis()->SetRangeUser(phcnt_min_xaxis, phcnt_max_xaxis);
            hCorrPhCount[counter]->GetYaxis()->SetRangeUser(phcnt_min_xaxis, phcnt_max_xaxis);
            hCorrPhCount[counter]->SetStats(false);
            hCorrPhCount[counter]->Draw("colz");
            
            fdiag = new TF1("fdiag", "x[0]", phcnt_min_xaxis, phcnt_max_xaxis);
            fdiag->Draw("same");
            
            text.DrawLatex(0.15, 0.85, string);
            
            counter++;
        }
    }
    
    //----- saving
    file->cd();
    
    for (auto c : can_phcnt_corr)
        c->Write();
    
    dirName = "PhCountCorrSpecta";
    if (!file->cd(dirName))
        TDirectory *dir_2 = file->mkdir(dirName);
    
    file->cd(dirName);
    
    for (auto h : hCorrPhCount)
        h->Write();
    
    //----- deleting
    for (auto c : can_phcnt_corr)
        delete c;
    
    for (auto h : hCorrPhCount)
        delete h;
    
    /*********/ //----- T0 spectra -----//
    
    counter = 0;
    
    std::vector<TH1D*> hT0Ch0 = data->GetSpectra(0, SFSelectionType::kPMIT0, cutCh0);
    std::vector<TH1D*> hT0Ch1 = data->GetSpectra(1, SFSelectionType::kPMIT0, cutCh1);
    
    std::vector<TCanvas*> can_t0(ncans);
    
    double t0_min_yaxis = 0.;
    double t0_max_yaxis = SFTools::FindMaxYaxis(hT0Ch0[0]);
    double t0_min_xaxis = 0.;
    double t0_max_xaxis = SFTools::FindMaxXaxis(hT0Ch0[0]);
    
    for (int ican = 0; ican < ncans; ican++)
    {
        TString cname = Form("data_t0_pt%i", ican+1);
        can_t0[ican] = new TCanvas(cname, cname, 2000, 1200);
        can_t0[ican]->DivideSquare(npoints/ncans);
        
        for (int i = 0; i< (npoints/ncans); i++)
        {
            can_t0[ican]->cd(i+1);
            gPad->SetGrid(1, 1);
            
            stringCh0 = hT0Ch0[counter]->GetTitle();
            stringCh1 = hT0Ch1[counter]->GetTitle();
            
            hT0Ch0[counter]->SetTitle(Form("T0 spectrum, fiber %.0f", fibers[counter]));
            hT0Ch0[counter]->GetXaxis()->SetTitle("time [ns]");
            hT0Ch0[counter]->GetYaxis()->SetTitle("counts");
            hT0Ch0[counter]->GetYaxis()->SetRangeUser(t0_min_yaxis, t0_max_yaxis);
            hT0Ch0[counter]->GetYaxis()->SetRangeUser(t0_min_yaxis, t0_max_yaxis);
            hT0Ch0[counter]->SetStats(false);
            hT0Ch0[counter]->SetLineColor(colCh0);
            hT0Ch0[counter]->Draw();
            
            hT0Ch1[counter]->SetStats(false);
            hT0Ch1[counter]->SetLineColor(colCh1);
            hT0Ch1[counter]->Draw("same");
            
            textCh0.DrawLatex(0.2, 0.3, stringCh0);
            textCh1.DrawLatex(0.2, 0.3, stringCh1);
        }
    }
    
    //----- saving 
    
    file->cd();
    
    for (auto c : can_t0)
        c->Write();
    
    dirName = "T0Spectra";
    if (!file->cd(dirName))
        TDirectory *dir_3 = file->mkdir(dirName);
    
    file->cd(dirName);
    
    for (auto h : hT0Ch0)
        h->Write();
    
    for (auto h : hT0Ch1)
        h->Write();
    
    //----- deleting
    
    for (auto c : can_t0)
        delete c;
    
    for (auto h : hT0Ch0)
        delete h;
    
    for (auto h : hT0Ch1)
        delete h;
    
    /*********/ //----- T0 correlation spectra -----//
    
    counter = 0;
    
    std::vector<TH2D*> hCorrT0 = data->GetCorrHistograms(SFSelectionType::kPMIT0Correlation, cutCh0Ch1);
    
    std::vector<TCanvas*> can_corr_t0(ncans);
    
    for (int ican = 0; ican < ncans; ican++)
    {
        TString cname = Form("data_t0_corr_pt%i", ican+1);
        can_corr_t0[ican] = new TCanvas(cname, cname, 2000, 1200);
        can_corr_t0[ican]->DivideSquare(npoints/ncans);
        
        for (int i = 0; i < (npoints/ncans); i++)
        {
            can_corr_t0[ican]->cd(i+1);
            gPad->SetGrid(1, 1);
            
            string = hCorrT0[counter]->GetTitle();
            
            hCorrT0[counter]->SetTitle(Form("T0 correlation spectrum, fiber %.0f", fibers[counter]));
            hCorrT0[counter]->GetXaxis()->SetTitle("Ch1 T0 [ns]");
            hCorrT0[counter]->GetYaxis()->SetTitle("Ch0 T0 [ns]");
            hCorrT0[counter]->GetXaxis()->SetRangeUser(t0_min_xaxis, t0_max_xaxis);
            hCorrT0[counter]->GetYaxis()->SetRangeUser(t0_min_xaxis, t0_max_xaxis);
            hCorrT0[counter]->SetStats(false);
            hCorrT0[counter]->Draw("colz");
            
            text.DrawLatex(0.15, 0.85, string);
            
            counter++;
        }
    }
    
    //----- saving
    
    file->cd();
    
    for (auto c : can_corr_t0)
        c->Write();
    
    dirName = "T0CorrSpectra";
    if (!file->cd(dirName))
        TDirectory *dir_4 = file->mkdir(dirName);
    
    file->cd(dirName);
    
    for (auto h : hCorrT0)
        h->Write();
    
    //----- deleting 
    
    for (auto c : can_corr_t0)
        delete c;
    
    for (auto h : hCorrT0)
        delete h;
    
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
        std::cout << "----- data_mod writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- data_mod writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    
    return 0;
}
