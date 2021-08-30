// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              peakfin.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFAttenuation.hh"
#include "SFData.hh"
#include "SFPeakFinder.hh"
#include "common_options.h"

// #include <DistributionContext.h>

#include <TCanvas.h>
#include <TLatex.h>

#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{

    TString outdir;
    TString dbase;
    int     seriesNo = -1;

    int ret = parse_common_options(argc, argv, outdir, dbase, seriesNo);
    if (ret != 0) exit(ret);

    if (argc < 2)
    {
        std::cout << "to run type: ./posres seriesNo";
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
        std::cerr << "##### Exception in peakfin.cc!" << std::endl;
        return 1;
    }

    data->Print();

    double              s         = SFTools::GetSigmaBL(data->GetSiPM());
    std::vector<double> sigmas    = {s, s};
    TString             cutCh0    = SFDrawCommands::GetCut(SFCutType::kSpecCh0, sigmas);
    TString             cutCh1    = SFDrawCommands::GetCut(SFCutType::kSpecCh1, sigmas);
    TString             cutCh0Ch1 = SFDrawCommands::GetCut(SFCutType::kCombCh0Ch1, sigmas);

    int                 npoints         = data->GetNpoints();
    //int               anaGroup        = data->GetAnalysisGroup();
    TString             desc            = data->GetDescription();
    std::vector<double> positions       = data->GetPositions();
    std::vector <int>   ID              = data->GetMeasurementsIDs();
    std::vector<int>    measurementsIDs = data->GetMeasurementsIDs();
    std::vector<TH1D*>  hSpecCh0        = data->GetSpectra(0, SFSelectionType::kPE, cutCh0);
    std::vector<TH1D*>  hSpecCh1        = data->GetSpectra(1, SFSelectionType::kPE, cutCh1);
    std::vector<TH1D*>  hSpecAv = data->GetCustomHistograms(SFSelectionType::kPEAverage, cutCh0Ch1);
    
    std::vector<TH1D*> hPeakCh0(npoints);
    std::vector<TH1D*> hPeakCh1(npoints);
    std::vector<TH1D*> hPeakAv(npoints);

    std::vector<SFPeakFinder*> pfCh0;
    std::vector<SFPeakFinder*> pfCh1;
    std::vector<SFPeakFinder*> pfAve;

//     DistributionContext ctx;
//     ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));
// 
//     ctx.dim   = DIM1;
//     ctx.x.min = 0;
//     ctx.x.max = 100;
//     ctx.y.min = 0;
//     ctx.y.max = 1000;

    for (int i = 0; i < npoints; i++)
    {
        pfAve.push_back(new SFPeakFinder(hSpecAv[i], 0, 1)); // verbose = 0, tests = 1
        pfAve[i]->FindPeakFit();
        pfAve[i]->SubtractBackground();
        pfAve[i]->GetResults()->Print();
        hPeakAv[i] = (TH1D*)pfAve[i]->GetResults()->GetObject(SFResultTypeObj::kPeak);
        
        pfCh0.push_back(new SFPeakFinder(hSpecCh0[i], 0, 1));
        pfCh0[i]->FindPeakFit();
        pfCh0[i]->SubtractBackground();
        pfCh0[i]->GetResults()->Print();
        hPeakCh0[i] = (TH1D*)pfCh0[i]->GetResults()->GetObject(SFResultTypeObj::kPeak);
        
        pfCh1.push_back(new SFPeakFinder(hSpecCh1[i], 0, 1));
        pfCh1[i]->FindPeakFit();
        pfCh1[i]->SubtractBackground();
        pfCh1[i]->GetResults()->Print();
        hPeakCh1[i] = (TH1D*)pfCh1[i]->GetResults()->GetObject(SFResultTypeObj::kPeak);
    }

    //----- results of fitting
    TCanvas* canCh0 = new TCanvas("pf_ch0", "pf_ch0", 2000, 1200);
    canCh0->DivideSquare(npoints);

    TCanvas* canCh1 = new TCanvas("pf_ch1", "pf_ch1", 2000, 1200);
    canCh1->DivideSquare(npoints);

    TCanvas* canAve = new TCanvas("pf_ave", "pf_ave", 2000, 1200);
    canAve->DivideSquare(npoints);

    //----- results of background subtracting
    TCanvas* canCh0_bgs = new TCanvas("pf_ch0_bgs", "pf_ch0_bgs", 2000, 1200);
    canCh0_bgs->DivideSquare(npoints);

    TCanvas* canCh1_bgs = new TCanvas("pf_ch1_bgs", "pf_ch1_bgs", 2000, 1200);
    canCh1_bgs->DivideSquare(npoints);

    TCanvas* canAve_bgs = new TCanvas("pf_ave_bgs", "pf_ave_bgs", 2000, 1200);
    canAve_bgs->DivideSquare(npoints);

    //----- drawing

    //double maxval = 0;
    double min_yaxis    = 0;
    double max_yaxis_0  = SFTools::FindMaxYaxis(hSpecCh0[npoints-1]);
    double max_yaxis_1  = SFTools::FindMaxYaxis(hSpecCh1[0]);
    double max_yaxis_av = SFTools::FindMaxYaxis(hSpecAv[4]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(hSpecCh0[0]);
    
    TLatex text;
    text.SetNDC(true);
    
    TString fun_name; 
   // TF1* fun = nullptr;
    
    for (int i = 0; i < npoints; i++)
    {
        canCh0->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecCh0[i]->SetStats(0);
        hSpecCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        hSpecCh0[i]->GetYaxis()->SetTitle("counts");
        hSpecCh0[i]->SetTitle(Form("PE spectrum: S%i ch0 %.1f mm", seriesNo, positions[i]));
        //maxval = hSpecCh0[i]->GetBinContent(hSpecCh0[i]->GetMaximumBin());
        //hSpecCh0[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
        //ctx.configureFromJson("hSpec");
        //ctx.print();
        //hSpecCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //hSpecCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        hSpecCh0[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        hSpecCh0[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        hSpecCh0[i]->DrawClone();
        
        SFResults* resultCh0 = pfCh0[i]->GetResults();

        text.DrawLatex(0.60, 0.75, Form("c = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakConst), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.60, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.60, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.60, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kChi2NDF)));

        canCh0_bgs->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecCh0[i]->Draw();

        hSpecCh0[i]->GetFunction("fun_bg_clone")->Delete();
        hSpecCh0[i]->GetFunction("fun_gaus_clone")->Delete();
        hPeakCh0[i]->SetLineColor(kMagenta);
        hPeakCh0[i]->Draw("same");


        text.DrawLatex(0.65, 0.75, Form("c = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakConst), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.65, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.65, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultCh0->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.65, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultCh0->GetValue(SFResultTypeNum::kChi2NDF)));

        canCh1->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecCh1[i]->SetStats(0);
        hSpecCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
        hSpecCh1[i]->GetYaxis()->SetTitle("counts");
        hSpecCh1[i]->SetTitle(Form("PE spectrum: S%i ch1 %.1f mm", seriesNo, positions[i]));
        //maxval = hSpecCh1[i]->GetBinContent(hSpecCh1[i]->GetMaximumBin());
        //hSpecCh1[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
        //hSpecCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //hSpecCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        hSpecCh1[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        hSpecCh1[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        hSpecCh1[i]->DrawClone();
        
        SFResults* resultCh1 = pfCh1[i]->GetResults();

        text.DrawLatex(0.60, 0.75, Form("c = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakConst), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.60, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.60, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.60, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kChi2NDF)));

        canCh1_bgs->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecCh1[i]->Draw();

        hSpecCh1[i]->GetFunction("fun_bg_clone")->Delete();
        hSpecCh1[i]->GetFunction("fun_gaus_clone")->Delete();
        hPeakCh1[i]->SetLineColor(kMagenta);
        hPeakCh1[i]->Draw("same");
        
        text.DrawLatex(0.60, 0.75, Form("c = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakConst), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.60, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.60, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultCh1->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.60, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultCh1->GetValue(SFResultTypeNum::kChi2NDF)));

        canAve->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecAv[i]->SetStats(0);
        hSpecAv[i]->GetXaxis()->SetTitle("charge [P.E.]");
        hSpecAv[i]->GetYaxis()->SetTitle("counts");
        hSpecAv[i]->SetTitle(Form("Average charge spectrum: S%i %.1f mm", seriesNo, positions[i]));
        //maxval = hSpecAv[i]->GetBinContent(hSpecAv[i]->GetMaximumBin());
        //hSpecAv[i]->GetYaxis()->SetRangeUser(-10, maxval+0.1*maxval);
        //ctx.configureFromJson("hSpecAv");
        //ctx.print();
        //hSpecAv[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //hSpecAv[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        hSpecAv[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        hSpecAv[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_av);
        hSpecAv[i]->DrawClone();
        
        SFResults* resultAv = pfAve[i]->GetResults();

        text.DrawLatex(0.60, 0.75, Form("c = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakConst), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.60, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.60, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.60, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kChi2NDF)));

        canAve_bgs->cd(i + 1);
        gPad->SetGrid(1, 1);
        hSpecAv[i]->Draw();

        hSpecAv[i]->GetFunction("fun_bg_clone")->Delete();
        hSpecAv[i]->GetFunction("fun_gaus_clone")->Delete();
        hPeakAv[i]->SetLineColor(kMagenta);
        hPeakAv[i]->Draw("same");
        
        text.DrawLatex(0.60, 0.75, Form("c = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakConst), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.60, 0.70, Form("#mu = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakPosition), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.60, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kPeakSigma), 
                       resultAv->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.60, 0.60, Form("#chi^{2}/NDF = %.3f", 
                       resultAv->GetValue(SFResultTypeNum::kChi2NDF)));
    }

    //----- saving
    TString fname       = Form("peakfin_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in peakfin.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    canCh0->Write();
    canCh1->Write();
    canAve->Write();
    canCh0_bgs->Write();
    canCh1_bgs->Write();
    canAve_bgs->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "PEAK_FINDER";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES(%i, '%s')",
                         table.Data(), seriesNo, fname_full.Data());

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- peakfin writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- peakfin writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    // delete att;

    return 0;
}
