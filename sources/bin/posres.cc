// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               posres.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2019              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFPositionRes.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>

#include <DistributionContext.h>

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
        std::cerr << "##### Exception in attenuation.cc!" << std::endl;
        return 1;
    }

    data->Print();

    TString desc = data->GetDescription();
    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in attenuation.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cout << "Description: " << desc << std::endl;
        return 1;
    }

    int npoints  = data->GetNpoints();
    int anaGroup = data->GetAnalysisGroup();

    SFPositionRes* posres;

    try
    {
        posres = new SFPositionRes(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in posres.cc!" << std::endl;
        return 1;
    }

    DistributionContext ctx;
    ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

    ctx.dim   = DIM1;
    ctx.x.min = 0;
    ctx.x.max = 100;
    ctx.y.min = 0;
    ctx.y.max = 1000;

    posres->AnalyzePositionRes();

    SFResults* results = posres->GetResults();
    posres->Print();

    TGraphErrors* gPosRecoVsPos = (TGraphErrors*)results->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPosResVsPos = (TGraphErrors*)results->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gAttenuation = (TGraphErrors*)results->GetObject(SFResultTypeObj::kMLRvsPosGraph);
    TGraphErrors* gResiduals   = (TGraphErrors*)results->GetObject(SFResultTypeObj::kResidualGraph);

    double* posResAll    = gPosResVsPos->GetY();
    double* posResAllErr = gPosResVsPos->GetEY();
    double* posReco      = gPosRecoVsPos->GetY();
    double* posRecoErr   = gPosRecoVsPos->GetEY();

    std::vector<TH1D*> spec     = posres->GetSpectra();
    std::vector<TH1D*> hPosReco = posres->GetPositionRecoDist();

    std::vector<SFPeakFinder*> peakFinAv;

    std::vector<double> xmin(npoints);
    std::vector<double> xmax(npoints);

    for (int i = 0; i < npoints; i++)
    {
        peakFinAv.push_back(new SFPeakFinder(spec[i], 0));
        peakFinAv[i]->FindPeakRange(xmin[i], xmax[i]);
    }

    TLatex text;
    text.SetNDC(true);

    TLine line;
    line.SetLineColor(kRed);
    line.SetLineWidth(1);

    TCanvas* can_posreco_dist = new TCanvas("pr_posreco_dist", "pr_posreco_dist", 2000, 1200);
    can_posreco_dist->DivideSquare(npoints);

    TCanvas* can_spec = new TCanvas("pr_spec", "pr_spec", 2000, 1200);
    can_spec->DivideSquare(npoints);

    for (int i = 0; i < npoints; i++)
    {
        can_posreco_dist->cd(i + 1);
        gPad->SetGrid(1, 1);
        hPosReco[i]->Draw();
        text.DrawLatex(0.3, 0.8, Form("#mu = (%.2f +/- %.2f) mm", posReco[i], posRecoErr[i]));
        text.DrawLatex(0.3, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", posResAll[i], posResAllErr[i]));

        can_spec->cd(i + 1);
        gPad->SetGrid(1, 1);
        spec[i]->Draw();
        ctx.configureFromJson("hSpecAv");
        //     ctx.print();
        spec[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        spec[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        line.DrawLine(xmin[i], ctx.y.min, xmin[i], ctx.y.max);
        line.DrawLine(xmax[i], ctx.y.min, xmax[i], ctx.y.max);
    }

    TCanvas* can_posreco = new TCanvas("pr_posreco", "pr_posreco", 1000, 700);
    TPad*    pad_posreco = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0);
    TPad*    pad_res     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);

    can_posreco->cd(1);
    pad_posreco->Draw();
    pad_posreco->cd();
    pad_posreco->SetGrid(1, 1);
    gPosRecoVsPos->Draw("AP");
    // gPosRecoVsPos->GetXaxis()->SetLabelSize(0.035);
    // gPosRecoVsPos->GetYaxis()->SetLabelSize(0.035);
    // gPosRecoVsPos->GetYaxis()->SetTitleOffset(1.4);

    TF1* funpol1 = gPosRecoVsPos->GetFunction("funpol1");

    text.SetTextSize(0.03);
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.2f +/- %.2f", 
                   funpol1->GetParameter(0), funpol1->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.2f +/- %.2f", 
                   funpol1->GetParameter(1), funpol1->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   funpol1->GetChisquare() / funpol1->GetNDF()));

    can_posreco->cd(0);
    pad_res->Draw();
    pad_res->cd();
    pad_res->SetGrid(1, 1);
    gResiduals->Draw("AP");

    TCanvas* can_posres = new TCanvas("pr_posres", "pr_pores", 700, 500);
    can_posres->cd();

    text.SetTextSize(0.05);
    can_posres->cd(2);
    gPad->SetGrid(1, 1);
    gPosResVsPos->Draw("AP");
    // gPosResVsPos->GetXaxis()->SetLabelSize(0.035);
    // gPosResVsPos->GetYaxis()->SetLabelSize(0.035);
    // gPosResVsPos->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.3, 0.8, Form("PR = (%.2f +/- %.2f) mm",
                   results->GetValue(SFResultTypeNum::kPositionRes),
                   results->GetUncertainty(SFResultTypeNum::kPositionRes)));

    TCanvas* can_att = new TCanvas("pr_att", "pr_att", 700, 500); 
    gPad->SetGrid(1, 1);
    gAttenuation->Draw("AP");
    TF1* fPol3 = (TF1*)gAttenuation->GetFunction("funpol3");
    text.SetTextSize(0.025);
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x + p_{2}x^{2} + p_{3}x^{3}");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.4e +/- %.4e", 
                   fPol3->GetParameter(0), fPol3->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.4e +/- %.4e", 
                   fPol3->GetParameter(1), fPol3->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("p_{2} = %.4e +/- %.4e",
                   fPol3->GetParameter(2), fPol3->GetParError(2)));
    text.DrawLatex(0.2, 0.60, Form("p_{3} = %.4e +/- %.4e",
                   fPol3->GetParameter(3), fPol3->GetParError(3)));
    text.DrawLatex(0.2, 0.55, Form("#chi^{2} = %.2f", fPol3->GetChisquare()));

    TCanvas* can_fun = new TCanvas("pr_fun", "pr_fun", 700, 500);

    TH1F* htemp = new TH1F("hPol3", "hPol3", 100, -1, 1);
    gPad->SetGrid(1, 1);
    htemp->Draw("func");
    fPol3->Draw("same");
    htemp->SetStats(0);
    htemp->GetXaxis()->SetRangeUser(-1, 1);
    htemp->GetYaxis()->SetRangeUser(fPol3->Eval(-1), fPol3->Eval(1));
    htemp->GetYaxis()->SetTitle("source position [mm]");
    htemp->GetXaxis()->SetTitle("ln(M_{LR})");
    htemp->SetTitle(Form("Calibration curve for position reconstruction S%i", seriesNo));

    //----- saving
    TString fname       = Form("posres_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in posres.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_posreco_dist->Write();
    can_posreco->Write();
    can_spec->Write();
    can_posres->Write();
    can_att->Write();
    can_fun->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "POSITION_RESOLUTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES, "
                    "POSITION_RES_ERR) VALUES(%i, '%s', %f, %f)",
                    table.Data(), seriesNo, fname_full.Data(), 
                    results->GetValue(SFResultTypeNum::kPositionRes),
                    results->GetUncertainty(SFResultTypeNum::kPositionRes)); 

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- posres writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- posres writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    for (auto h : hPosReco)
        delete h;
    
    for (auto h : spec)
        delete h;
    
    for (auto pf : peakFinAv)
        delete pf;
    
    delete data;
    delete posres;

    return 0;
}
