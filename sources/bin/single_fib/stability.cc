// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             stability.cc              *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFStabilityMon.hh"
#include "SFTools.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>

// #include <DistributionContext.h>

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
        std::cout << "to run type: ./stability seriesNo";
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
        std::cerr << "##### Exception in stability.cc!" << std::endl;
        return 1;
    }

    data->Print();
    int                 npoints   = data->GetNpoints();
    //int                 anaGroup  = data->GetAnalysisGroup();
    std::vector<double> positions = data->GetPositions();
    std::vector<int>    ID        = data->GetMeasurementsIDs();

    TString desc = data->GetDescription();
    if (!desc.Contains("Stability monitoring"))
    {
        std::cerr << "##### Error in stability.cc! This is not stability monitoring series!"
                  << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << desc << std::endl;
        return 1;
    }

    SFStabilityMon* stab;
    try
    {
        stab = new SFStabilityMon(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Error in stability.cc" << std::endl;
        return 1;
    }
/*
    DistributionContext ctx;
    ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

    ctx.dim   = DIM1;
    ctx.x.min = 0;
    ctx.x.max = 100;
    ctx.y.min = 0;
    ctx.y.max = 1000;
*/
    //----- stability analysis
    stab->AnalyzeStability(0);
    stab->AnalyzeStability(1);

    //----- accessing results
    std::vector<SFResults*> results = stab->GetResults();
    results[0]->Print(); //results for channel 0
    results[1]->Print(); //results for channel 1

    TGraphErrors* gCh0PeakPos = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPeakPosGraph);
    TGraphErrors* gCh1PeakPos = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPeakPosGraph);

    TGraphErrors* gCh0Residual = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kSMResidualGraph);
    TGraphErrors* gCh1Residual = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kSMResidualGraph);

    std::vector<TH1D*> specCh1 = stab->GetSpectra(1);
    std::vector<TH1D*> specCh0 = stab->GetSpectra(0);

    TCanvas* can         = new TCanvas("stab", "stab", 1000, 700);
    TPad*    pad_peakPos = new TPad("pad_peakPos", "pad_peakPos", 0, 0.3, 1, 1, 10, 0);
    TPad*    pad_res     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);

    double minCh0, minCh1;
    double maxCh0, maxCh1;

    can->cd(0);
    pad_peakPos->Draw();
    pad_peakPos->cd();
    pad_peakPos->SetGrid(1, 1);
    gCh0PeakPos->Draw("AP");
    gCh0PeakPos->SetTitle("511 keV peak position stability");
    gCh0PeakPos->SetMarkerColor(kPink - 8);
    gCh0PeakPos->SetLineColor(kPink - 8);
    gCh0PeakPos->GetFunction("funPol0")->SetLineColor(kPink - 8);
    minCh0 = TMath::MinElement(npoints, gCh0PeakPos->GetY());
    maxCh0 = TMath::MaxElement(npoints, gCh0PeakPos->GetY());

    gCh1PeakPos->Draw("P");
    gCh1PeakPos->SetMarkerColor(kAzure - 6);
    gCh1PeakPos->SetLineColor(kAzure - 6);
    gCh1PeakPos->GetFunction("funPol0")->SetLineColor(kAzure - 6);
    minCh1 = TMath::MinElement(npoints, gCh1PeakPos->GetY());
    maxCh1 = TMath::MaxElement(npoints, gCh1PeakPos->GetY());

    double min = TMath::Min(minCh0, minCh1);
    double max = TMath::Max(maxCh0, maxCh1);
    gCh0PeakPos->GetYaxis()->SetRangeUser(min - 2, max + 2);

    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.04);
    text.SetTextColor(kPink - 8);
    text.DrawLatex(0.2, 0.5,
                   Form("#bar{PP}_{ch0} = %.2f +/- %.2f",
                        results[0]->GetValue(SFResultTypeNum::kAveragePeakPos),
                        results[0]->GetUncertainty(SFResultTypeNum::kAveragePeakPos)));
    text.SetTextColor(kAzure - 6);
    text.DrawLatex(0.2, 0.4,
                   Form("#bar{PP}_{ch1} = %.2f +/- %.2f",
                        results[1]->GetValue(SFResultTypeNum::kAveragePeakPos),
                        results[1]->GetUncertainty(SFResultTypeNum::kAveragePeakPos)));

    can->cd(0);
    pad_res->Draw();
    pad_res->cd();
    pad_res->SetGrid(1, 1);
    gCh0Residual->Draw("AP");
    gCh0Residual->SetTitle("Residuals graph");
    gCh0Residual->GetYaxis()->SetRangeUser(-2, 2);
    gCh0Residual->SetMarkerColor(kPink - 8);
    gCh0Residual->SetLineColor(kPink - 8);
    gCh1Residual->Draw("P");
    gCh1Residual->SetMarkerColor(kAzure - 6);
    gCh1Residual->SetLineColor(kAzure - 6);

    TF1* funPol0 = new TF1("funPol0", "pol0", 0, npoints);
    funPol0->FixParameter(0, 0);
    funPol0->SetLineColor(kGray + 2);
    funPol0->SetLineStyle(9);
    gCh0Residual->Fit(funPol0, "Q");

    TCanvas* can_ch0 = new TCanvas("stab_ch0", "stab_ch0", 2000, 1200);
    can_ch0->DivideSquare(npoints);

    TCanvas* can_ch1 = new TCanvas("stab_ch1", "stab_ch1", 2000, 1200);
    can_ch1->DivideSquare(npoints);

    double min_yaxis = 0.;
    double max_yaxis_0 =  SFTools::FindMaxYaxis(specCh0[npoints-1]);
    double max_yaxis_1 =  SFTools::FindMaxYaxis(specCh1[0]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(specCh0[0]);

    TString fun_name;
    
    for (int i = 0; i < npoints; i++)
    {
        can_ch0->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh0[i]->SetTitle(Form("PE spectrum, ch0, position %.2f", positions[i]));
        specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specCh0[i]->GetYaxis()->SetTitle("counts");
        specCh0[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh0[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        specCh0[i]->SetStats(false);
        //ctx.configureFromJson("hSpec");
        //ctx.print();
        //specCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        specCh0[i]->Draw();
        
        fun_name = Form("f_S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.5, 0.75, Form("#mu = %.2f +/- %.2f", 
                       specCh0[i]->GetFunction(fun_name)->GetParameter(1),
                       specCh0[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.5, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       specCh0[i]->GetFunction(fun_name)->GetParameter(2),
                       specCh0[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.5, 0.60, Form("#chi^{2}/NDF = %.3f",
                       specCh0[i]->GetFunction(fun_name)->GetChisquare() / 
                       specCh0[i]->GetFunction(fun_name)->GetNDF()));

        can_ch1->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh1[i]->SetTitle(Form("PE spectrum, ch1, position %.2f", positions[i]));
        specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specCh1[i]->GetYaxis()->SetTitle("counts");
        specCh1[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh1[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        //specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        specCh1[i]->SetStats(false);
        specCh1[i]->Draw();
        
        fun_name = Form("f_S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.5, 0.75, Form("#mu = %.2f +/- %.2f", 
                       specCh1[i]->GetFunction(fun_name)->GetParameter(1),
                       specCh1[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.5, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       specCh1[i]->GetFunction(fun_name)->GetParameter(2),
                       specCh1[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.5, 0.60, Form("#chi^{2}/NDF = %.3f",
                       specCh1[i]->GetFunction(fun_name)->GetChisquare() / 
                       specCh1[i]->GetFunction(fun_name)->GetNDF()));
    }

    //----- saving
    TString fname       = Form("stability_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in stability.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can->Write();
    can_ch0->Write();
    can_ch1->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "STABILITY_MON";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, CH0_MEAN, "
                         "CH0_STDDEV, CH1_MEAN, CH1_STDDEV) VALUES (%i, '%s', %f, %f, %f, "
                         "%f)",
                         table.Data(), seriesNo, fname_full.Data(),
                         results[0]->GetValue(SFResultTypeNum::kAveragePeakPos),
                         results[0]->GetUncertainty(SFResultTypeNum::kAveragePeakPos),
                         results[1]->GetValue(SFResultTypeNum::kAveragePeakPos),
                         results[1]->GetUncertainty(SFResultTypeNum::kAveragePeakPos));

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- stability writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- stability writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    delete stab;

    return 0;
}
