// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             counts_mod.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFCountMap.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMathBase.h>
#include <TMath.h>

#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{
    TString outdir;
    TString dbase;
    int seriesNo = -1;
    
    int ret = parse_common_options(argc, argv, outdir, dbase, seriesNo);
    if (ret != 0) exit(ret);

    if (argc < 2)
    {
        std::cout << "to run type: ./counts_mod seriesNo";
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
        std::cerr << "##### Exception in counts_mod.cc!" << std::endl;
        return 1;
    }

    data->Print();
    
    int     npoints             = data->GetNpoints();
    TString testBench           = data->GetTestBench();
    TString desc                = data->GetDescription();
    std::vector<double> fiberNo = data->GetPositions();
    
    if (! desc.Contains("Module series"))
    {
        std::cerr << "##### Error in counts_mod.cc! This is not module series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << desc << std::endl;
        return 1;
    }

    SFCountMap* counts = nullptr;
    try
    {
        counts = new SFCountMap(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in counts_mod.cc!" << std::endl;
        return 1;
    }
    
    counts->Print();
    
    counts->DrawCounts(0);
    counts->DrawCounts(1);
    
    std::vector<SFResults*> results = counts->GetResults();
    results[0]->Print(); //ch0
    results[1]->Print(); //ch1
    
    //----- getting results
    TGraphErrors *gCountsCh0 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kCountsGraph);
    TGraphErrors *gCountsCh1 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kCountsGraph);
    
    std::vector<TH1D*> hSpecCh0 = counts->GetSpectra(0);
    std::vector<TH1D*> hSpecCh1 = counts->GetSpectra(1);
    
    int col_ch0 = kPink - 8;
    int col_ch1 = kAzure - 6;
    
    //----- drawing counts
    TCanvas *can_counts = new TCanvas("counts", "counts", 700, 500);
    gPad->SetGrid(1, 1);
    
    gCountsCh0->SetMarkerColor(col_ch0);
    gCountsCh0->SetLineColor(col_ch0);
    gCountsCh0->Draw("AP");
    
    gCountsCh1->SetMarkerColor(col_ch1);
    gCountsCh1->SetLineColor(col_ch1);
    gCountsCh1->Draw("P");
    
    double ymin = TMath::Min(TMath::MinElement(npoints, gCountsCh0->GetY()),
                             TMath::MinElement(npoints, gCountsCh1->GetY())); 
    double ymax = TMath::Max(TMath::MaxElement(npoints, gCountsCh0->GetY()),
                             TMath::MaxElement(npoints, gCountsCh1->GetY())); 
    
    gCountsCh0->GetYaxis()->SetRangeUser(ymin-0.1*ymin, ymax+0.1*ymax);
    
    TLegend *leg = new TLegend(0.68, 0.77, 0.9, 0.9);
    leg->AddEntry(gCountsCh0, "channel 0", "PE");
    leg->AddEntry(gCountsCh1, "channel 1", "PE");
    leg->Draw();
    
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.04);
    
    text.SetTextColor(col_ch0);
    text.DrawLatex(0.4, 0.85, Form("#bar{c} = %.0f #sigma_{c} = %.0f",
                   results[0]->GetValue(SFResultTypeNum::kCounts),
                   results[0]->GetValue(SFResultTypeNum::kCountsStdDev)));
    text.SetTextColor(col_ch1);
    text.DrawLatex(0.4, 0.8, Form("#bar{c} = %.0f #sigma_{c} = %.0f",
                   results[1]->GetValue(SFResultTypeNum::kCounts),
                   results[1]->GetValue(SFResultTypeNum::kCountsStdDev)));
    
    //----- drawing spectra
    const int ncans   = 4;
    int       counter = 0; 
    
    std::vector<TCanvas*> can_ch0(ncans);
    std::vector<TCanvas*> can_ch1(ncans);
    
    text.SetTextColor(kBlack);
    
    for(int c=0; c<ncans; c++)
    {
        TString cname = Form("counts_ch0_pt%i", c+1);
        can_ch0[c] = new TCanvas(cname, cname, 2000, 1200);
        can_ch0[c]->DivideSquare(npoints/ncans);
        
        cname = Form("counts_ch1_pt%i", c+1);
        can_ch1[c] = new TCanvas(cname, cname, 2000, 1200);
        can_ch1[c]->DivideSquare(npoints/ncans);
        
            for(int i=0; i<(npoints/ncans); i++)
            {
                can_ch0[c]->cd(i+1);
                gPad->SetGrid(1, 1);
                hSpecCh0[counter]->SetStats(0);
                hSpecCh0[counter]->GetXaxis()->SetTitle("photon count");
                hSpecCh0[counter]->GetYaxis()->SetTitle("counts");
                hSpecCh0[counter]->SetTitle(Form("Spectrum fib%.0f ch0", fiberNo[counter]));
                hSpecCh0[counter]->Draw();
                text.DrawLatex(0.5, 0.8, Form("c = %.0f  err = %.0f",
                               gCountsCh0->GetPointY(counter),
                               gCountsCh0->GetErrorY(counter)));
                
                can_ch1[c]->cd(i+1);
                gPad->SetGrid(1, 1);
                hSpecCh1[counter]->SetStats(0);
                hSpecCh1[counter]->GetXaxis()->SetTitle("photon count");
                hSpecCh1[counter]->GetYaxis()->SetTitle("counts");
                hSpecCh1[counter]->SetTitle(Form("Spectrum fib%.0f ch1", fiberNo[counter]));
                hSpecCh1[counter]->Draw();
                text.DrawLatex(0.5, 0.8, Form("c = %.0f  err = %.0f",
                               gCountsCh1->GetPointY(counter),
                               gCountsCh1->GetErrorY(counter)));
                
                counter++;
            }
    }
    
    //----- saving
    TString fname       = Form("counts_mod_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in counts_mod.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }
    
    can_counts->Write();
    for (auto h : can_ch0)
        h->Write();
    for (auto h : can_ch1)
        h->Write();

    file->Close();
    
    //----- writing results to the data base
    TString table = "COUNTS";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, "
                         "COUNTS_CH0, COUNTS_ERR_CH0, COUNTS_STDDEV_CH0, COUNTS_CH1, " "COUNTS_ERR_CH1, COUNTS_STDDEV_CH1) VALUES (%i, '%s', %f, %f, " "%f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), 
                         results[0]->GetValue(SFResultTypeNum::kCounts),
                         results[0]->GetUncertainty(SFResultTypeNum::kCounts),
                         results[0]->GetValue(SFResultTypeNum::kCountsStdDev),
                         results[1]->GetValue(SFResultTypeNum::kCounts),
                         results[1]->GetUncertainty(SFResultTypeNum::kCounts),
                         results[1]->GetValue(SFResultTypeNum::kCountsStdDev));
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;
    
    do
    {
        std::cout << "----- counts_mod writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
        
        if (stat) break;
        
        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- counts_mod writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    for (auto h : hSpecCh0)
        delete h;
    
    for (auto h : hSpecCh1)
        delete h;
    
    delete data;
    delete counts;
    
    return 0;
}
