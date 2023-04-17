// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              timeres.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFTimingRes.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

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
        std::cout << "to run type: ./timeres seriesNo";
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
        std::cout << "##### Exception in timeres.cc!" << std::endl;
        return 1;
    }

    TString desc = data->GetDescription();
    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in timeres.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << desc << std::endl;
        return 1;
    }

    int                 npoints    = data->GetNpoints();
    //int                 anaGroup   = data->GetAnalysisGroup();
    TString             collimator = data->GetCollimator();
    TString             sipm       = data->GetSiPM();
    std::vector<double> positions  = data->GetPositions();
    std::vector<int>    ID         = data->GetMeasurementsIDs();
    data->Print();

    SFTimingRes* timeres;

    try
    {
        timeres = new SFTimingRes(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in timeres.cc!" << std::endl;
        return 1;
    }

    timeres->Print();
    timeres->AnalyzeNoECut();
    timeres->AnalyzeWithECut();
/*
    DistributionContext ctx;
    ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

    ctx.dim   = DIM1;
    ctx.x.min = 0;
    ctx.x.max = 100;
    ctx.y.min = 0;
    ctx.y.max = 1000;
*/
    //----- timing resolution, no enery cut
    std::vector<TH1D*> T0diff   = timeres->GetT0Diff(0);
    std::vector<TH1D*> ratio    = timeres->GetRatios();

    //----- timing resolution, with energy cut
    std::vector<TH1D*> T0diffECut   = timeres->GetT0Diff(1);
    std::vector<TH1D*> specCh0      = timeres->GetSpectra(0);
    std::vector<TH1D*> specCh1      = timeres->GetSpectra(1);

    //----- results
    std::vector<SFResults*> results = timeres->GetResults();
    results[0]->Print(); // no cuts
    results[1]->Print(); // with energy cut
    
    TGraphErrors* gTimeRes = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kTimeResGraph);
    TGraphErrors* gTimeSig = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kTimeSigGraph);
    
    TGraphErrors* gTimeResECut = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kTimeResGraph);
    TGraphErrors* gTimeSigECut = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kTimeSigGraph);


    //----- drawing
    TCanvas* canTDiff = new TCanvas("tr_TDiff", "tr_TDiff", 2000, 1200);
    canTDiff->DivideSquare(npoints);

    TCanvas* canTDiffECut = new TCanvas("tr_TDiffECut", "tr_TDiffECut", 2000, 1200);
    canTDiffECut->DivideSquare(npoints);

    TCanvas* canRatio = new TCanvas("tr_ratio", "tr_ratio", 2000, 1200);
    canRatio->DivideSquare(npoints);

    TCanvas* canSpecCh0 = new TCanvas("tr_spec_ch0", "tr_spec_ch0", 2000, 1200);
    canSpecCh0->DivideSquare(npoints);

    TCanvas* canSpecCh1 = new TCanvas("tr_spec_ch1", "tr_spec_ch1", 2000, 1200);
    canSpecCh1->DivideSquare(npoints);

    TLatex  text;
    TString string;
    TString title;
    text.SetNDC(true);
    text.SetTextSize(0.045);

    TLine line;
    line.SetLineColor(kRed);
    line.SetLineStyle(9);

    std::vector<SFPeakFinder*> peakfinCh0;
    std::vector<SFPeakFinder*> peakfinCh1;

    double mean, sigma;
    double max;
    double xmin, xmax;
    double center, delta; // changed here for smaller cut

    TF1* fT0thin  = new TF1("fT0thin", "gaus", -50, 50);
    TF1* fT0thick = new TF1("fT0thick", "gaus", -50, 50);
    fT0thin->SetLineColor(kMagenta);
    fT0thick->SetLineColor(kMagenta - 10);

    TF1* fRthin  = new TF1("fRthin", "gaus", -50, 50);
    TF1* fRthick = new TF1("fRthick", "gaus", -50, 50);
    TF1* fClone  = new TF1("fColone", "gaus", -50, 50);
    fRthin->SetLineColor(kMagenta);
    fRthick->SetLineColor(kMagenta - 10);
    fClone->SetLineColor(kRed);

    double min_yaxis = 0.;
    double max_yaxis_0 =  SFTools::FindMaxYaxis(specCh0[npoints-1]);
    double max_yaxis_1 =  SFTools::FindMaxYaxis(specCh1[0]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(specCh0[0]);
    
    TString fun_name;

    for (int i = 0; i < npoints; i++)
    {
        title = Form("ch_0.fT0 - ch_1.fT0, S%i Source Position %.2f mm", 
                     seriesNo, positions[i]);

        canTDiff->cd(i + 1);
        gPad->SetGrid(1, 1);
        string = Form("TR = (%.3f +/- %.3f) ns",
                      gTimeSig->GetPointY(i),
                      gTimeSig->GetErrorY(i));
        T0diff[i]->SetTitle(title);
        T0diff[i]->GetXaxis()->SetTitle("ch_0.fT0-ch_1.fT0 [ns]");
        T0diff[i]->Draw();
        text.DrawLatex(0.15, 0.8, string);
        // if((collimator=="Lead" && sipm=="Hamamatsu") ||
        //   (collimator=="Electronic" && sipm=="SensL")){
        fT0thin->FixParameter(0, T0diff[i]->GetFunction("fun")->GetParameter(0));

        fT0thin->FixParameter(1, T0diff[i]->GetFunction("fun")->GetParameter(1));
        fT0thin->FixParameter(2, T0diff[i]->GetFunction("fun")->GetParameter(2));
        fT0thick->FixParameter(0, T0diff[i]->GetFunction("fun")->GetParameter(3));
        fT0thick->FixParameter(1, T0diff[i]->GetFunction("fun")->GetParameter(4));
        fT0thick->FixParameter(2, T0diff[i]->GetFunction("fun")->GetParameter(5));
        fT0thin->DrawClone("same");
        fT0thick->DrawClone("same");
        //}
        if (sipm == "Hamamatsu")
            T0diff[i]->GetXaxis()->SetRangeUser(-10, 10);
        else if (sipm == "SensL")
            T0diff[i]->GetXaxis()->SetRangeUser(-20, 20);

        canTDiffECut->cd(i + 1);
        gPad->SetGrid(1, 1);
        string = Form("TR = (%.3f +/- %.3f) ns",
                      gTimeSigECut->GetPointY(i),
                      gTimeSigECut->GetErrorY(i));
        T0diffECut[i]->SetTitle(title);
        T0diffECut[i]->GetXaxis()->SetTitle("ch_0.fT0-ch_1.fT0 [ns]");
        T0diffECut[i]->GetXaxis()->SetRangeUser(-20, 20);
        T0diffECut[i]->Draw();
        text.DrawLatex(0.15, 0.8, string);

        canRatio->cd(i + 1);
        gPad->SetGrid(1, 1);
        ratio[i]->GetXaxis()->SetTitle("ln(#sqrt{ch1/ch0})");
        ratio[i]->SetTitle(Form("ln(#sqrt{ch1/ch0}), Source Position %.2f mm", positions[i]));
        max = ratio[i]->GetBinContent(ratio[i]->GetMaximumBin());

        int parNo = -1;

        ratio[i]->Draw();

        if (collimator == "Lead" && sipm == "Hamamatsu")
        {
            if (ratio[i]->GetFunction("fDGauss")->GetParameter(0) >
                ratio[i]->GetFunction("fDGauss")->GetParameter(3))
                parNo = 1;
            else
                parNo = 4;
            mean  = ratio[i]->GetFunction("fDGauss")->GetParameter(parNo);
            sigma = ratio[i]->GetFunction("fDGauss")->GetParameter(parNo + 1);
            line.DrawLine(mean - 0.5 * sigma, 0, mean - 0.5 * sigma, max);
            line.DrawLine(mean + 0.5 * sigma, 0, mean + 0.5 * sigma, max);
            fRthin->FixParameter(0, ratio[i]->GetFunction("fDGauss")->GetParameter(0));
            fRthin->FixParameter(1, ratio[i]->GetFunction("fDGauss")->GetParameter(1));
            fRthin->FixParameter(2, ratio[i]->GetFunction("fDGauss")->GetParameter(2));
            fRthick->FixParameter(0, ratio[i]->GetFunction("fDGauss")->GetParameter(3));
            fRthick->FixParameter(1, ratio[i]->GetFunction("fDGauss")->GetParameter(4));
            fRthick->FixParameter(2, ratio[i]->GetFunction("fDGauss")->GetParameter(5));
            fRthin->DrawClone("same");
            fRthick->DrawClone("same");
        }
        else if (collimator == "Electronic" && sipm == "SensL")
        {
            if (ratio[i]->GetFunction("fDGauss")->GetParameter(0) >
                ratio[i]->GetFunction("fDGauss")->GetParameter(3))
                parNo = 1;
            else
                parNo = 4;
            mean  = ratio[i]->GetFunction("fDGauss")->GetParameter(parNo);
            sigma = ratio[i]->GetFunction("fDGauss")->GetParameter(parNo + 1);
            fRthin->FixParameter(0, ratio[i]->GetFunction("fDGauss")->GetParameter(0));
            fRthin->FixParameter(1, ratio[i]->GetFunction("fDGauss")->GetParameter(1));
            fRthin->FixParameter(2, ratio[i]->GetFunction("fDGauss")->GetParameter(2));
            fRthick->FixParameter(0, ratio[i]->GetFunction("fDGauss")->GetParameter(3));
            fRthick->FixParameter(1, ratio[i]->GetFunction("fDGauss")->GetParameter(4));
            fRthick->FixParameter(2, ratio[i]->GetFunction("fDGauss")->GetParameter(5));
            fRthin->DrawClone("same");
            fRthick->DrawClone("same");
            line.DrawLine(mean - 2 * sigma, 0, mean - 2 * sigma, max);
            line.DrawLine(mean + 2 * sigma, 0, mean + 2 * sigma, max);
        }
        else if (collimator == "Electronic" && sipm == "Hamamatsu")
        {
            mean  = ratio[i]->GetFunction("fGauss")->GetParameter(1);
            sigma = ratio[i]->GetFunction("fGauss")->GetParameter(2);
            fClone->FixParameter(0, ratio[i]->GetFunction("fGauss")->GetParameter(0));
            fClone->FixParameter(1, ratio[i]->GetFunction("fGauss")->GetParameter(1));
            fClone->FixParameter(2, ratio[i]->GetFunction("fGauss")->GetParameter(2));
            fClone->DrawClone("same");
            line.DrawLine(mean - 3 * sigma, 0, mean - 3 * sigma, max);
            line.DrawLine(mean + 3 * sigma, 0, mean + 3 * sigma, max);
        }

        ratio[i]->GetXaxis()->SetRangeUser(-1, 1);

        canSpecCh0->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh0[i]->GetXaxis()->SetTitle("charge P.E.");
        specCh0[i]->GetYaxis()->SetTitle("counts");
        specCh0[i]->SetTitle(
            Form("PE spectrum S%i Ch0, Source Position %.2f mm", seriesNo, positions[i]));
        specCh0[i]->SetStats(false);
        specCh0[i]->GetYaxis()->SetMaxDigits(2);
        peakfinCh0.push_back(new SFPeakFinder(specCh0[i], 0, 0));
        peakfinCh0[i]->FindPeakRange(xmin, xmax);
        center = xmin + (xmax - xmin) / 2.; // changed here for smaller cut
        delta = (xmax - xmin) / 6;          //
        //max   = specCh0[i]->GetBinContent(specCh0[i]->GetMaximumBin());
        specCh0[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh0[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        specCh0[i]->Draw();
        
        fun_name = Form("f_S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       specCh0[i]->GetFunction(fun_name)->GetParameter(1),
                       specCh0[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       specCh0[i]->GetFunction(fun_name)->GetParameter(2),
                       specCh0[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       specCh0[i]->GetFunction(fun_name)->GetChisquare() / 
                       specCh0[i]->GetFunction(fun_name)->GetNDF()));
        //ctx.configureFromJson("hSpecAv");
        //ctx.print();
        //specCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        if (collimator == "Lead")
        {
            line.DrawLine(center - delta, min_yaxis, center - delta, max_yaxis_0); // changed here for smaller cut
            line.DrawLine(center + delta, min_yaxis, center + delta, max_yaxis_0); //
        }
        else if (collimator == "Electronic")
        {
            // line.DrawLine(center-3*delta, 0, center-3*delta, max);  //changed here for smaller
            // cut line.DrawLine(center+3*delta, 0, center+3*delta, max);  //
            
            //line.DrawLine(center - 3 * delta, ctx.y.min, center - 3 * delta,
            //              ctx.y.max); // changed here for smaller cut
            //line.DrawLine(center + 3 * delta, ctx.y.min, center + 3 * delta, ctx.y.max); //
            
            line.DrawLine(center - 3 * delta, min_yaxis, center - 3 * delta, max_yaxis_0); // changed here for smaller cut
            line.DrawLine(center + 3 * delta, min_yaxis, center + 3 * delta, max_yaxis_0); //
        }

        canSpecCh1->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh1[i]->GetXaxis()->SetTitle("charge P.E.");
        specCh1[i]->GetYaxis()->SetTitle("counts");
        specCh1[i]->SetTitle(
            Form("PE spectrum S%i Ch1, Source Position %.2f mm", seriesNo, positions[i]));
        specCh1[i]->SetStats(false);
        specCh1[i]->GetYaxis()->SetMaxDigits(2);
        peakfinCh1.push_back(new SFPeakFinder(specCh1[i], 0, 0));
        peakfinCh1[i]->FindPeakRange(xmin, xmax);
        center = xmin + (xmax - xmin) / 2.; // changed here for smaller cut
        delta = (xmax - xmin) / 6;          //
        //max   = specCh1[i]->GetBinContent(specCh1[i]->GetMaximumBin());
        specCh1[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh1[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        specCh1[i]->Draw();
        
        fun_name = Form("f_S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       specCh1[i]->GetFunction(fun_name)->GetParameter(1),
                       specCh1[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       specCh1[i]->GetFunction(fun_name)->GetParameter(2),
                       specCh1[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       specCh1[i]->GetFunction(fun_name)->GetChisquare() / 
                       specCh1[i]->GetFunction(fun_name)->GetNDF()));
        //specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        if (collimator == "Lead")
        {
            line.DrawLine(center - delta, min_yaxis, center - delta, max_yaxis_1); // changed here for smaller cut
            line.DrawLine(center + delta, min_yaxis, center + delta, max_yaxis_1); //
        }
        else if (collimator == "Electronic")
        {
            // line.DrawLine(center-3*delta, 0, center-3*delta, max);  //changed here for smaller
            // cut line.DrawLine(center+3*delta, 0, center+3*delta, max);  //
            
            //line.DrawLine(center - 3 * delta, ctx.y.min, center - 3 * delta,
            //              ctx.y.max); // changed here for smaller cut
            //line.DrawLine(center + 3 * delta, ctx.y.min, center + 3 * delta, ctx.y.max); //
            
            line.DrawLine(center - 3 * delta, min_yaxis, center - 3 * delta, max_yaxis_1); // changed here for smaller cut
            line.DrawLine(center + 3 * delta, min_yaxis, center + 3 * delta, max_yaxis_1); //
        }
    }

    TCanvas* canTimingRes = new TCanvas("tr_TimingRes", "tr_TimingRes", 1200, 600);
    canTimingRes->Divide(2, 1);

    canTimingRes->cd(1);
    gPad->SetGrid(1, 1);
    gTimeRes->Draw("AP");
    gTimeRes->GetXaxis()->SetLabelSize(0.035);
    gTimeRes->GetYaxis()->SetLabelSize(0.035);
    gTimeRes->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.15, 0.8, Form("TR = (%.4f +/- %.4f) ns", 
                   results[0]->GetValue(SFResultTypeNum::kTimeRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kTimeRes)));

    canTimingRes->cd(2);
    gPad->SetGrid(1, 1);
    gTimeResECut->Draw("AP");
    gTimeResECut->GetXaxis()->SetLabelSize(0.035);
    gTimeResECut->GetYaxis()->SetLabelSize(0.035);
    gTimeResECut->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.15, 0.8, Form("TR = (%.4f +/- %.4f) ns",
                   results[1]->GetValue(SFResultTypeNum::kTimeRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kTimeRes)));

    TCanvas *canTimeSig = new TCanvas("tr_TimeSigma", "tr_TimeSigma", 1200, 600);
    canTimeSig->Divide(2, 1);
    
    canTimeSig->cd(1);
    gPad->SetGrid(1, 1);
    gTimeSig->Draw("AP");
    gTimeSig->GetXaxis()->SetLabelSize(0.035);
    gTimeSig->GetYaxis()->SetLabelSize(0.035);
    gTimeSig->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.15, 0.8, Form("TR = (%.4f +/- %.4f) ns", 
                   results[0]->GetValue(SFResultTypeNum::kTimeRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kTimeRes)));
    
    canTimeSig->cd(2);
    gPad->SetGrid(1, 1);
    gTimeSigECut->Draw("AP");
    gTimeSigECut->GetXaxis()->SetLabelSize(0.035);
    gTimeSigECut->GetYaxis()->SetLabelSize(0.035);
    gTimeSigECut->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.15, 0.8, Form("TR = (%.4f +/- %.4f) ns",
                   results[1]->GetValue(SFResultTypeNum::kTimeRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kTimeRes)));
    
    //----- saving
    TString fname       = Form("timeres_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in lightout.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    canTDiff->Write();
    canTDiffECut->Write();
    canRatio->Write();
    canSpecCh0->Write();
    canSpecCh1->Write();
    canTimingRes->Write();
    canTimeSig->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "TIMING_RESOLUTION";
    TString query =
        Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, TIMERES, TIMERES_ERR, "
             "TIMERES_ECUT, TIMERES_ECUT_ERR) VALUES(%i, '%s', %f, %f, %f, %f)",
             table.Data(), seriesNo, fname_full.Data(), 
             results[0]->GetValue(SFResultTypeNum::kTimeRes),
             results[0]->GetUncertainty(SFResultTypeNum::kTimeRes),
             results[1]->GetValue(SFResultTypeNum::kTimeRes),
             results[1]->GetUncertainty(SFResultTypeNum::kTimeRes));

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- timeres writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- timeres writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete timeres;
    delete data;

    return 0;
}
