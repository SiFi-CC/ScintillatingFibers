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
#include <TLegend.h>

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
    //int anaGroup = data->GetAnalysisGroup();
    std::vector<double> positions = data->GetPositions();
    std::vector<int> ID = data->GetMeasurementsIDs();
    TString colimator   = data->GetCollimator();
    TString testBench   = data->GetTestBench();

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
/*
    DistributionContext ctx;
    ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

    ctx.dim   = DIM1;
    ctx.x.min = 0;
    ctx.x.max = 100;
    ctx.y.min = 0;
    ctx.y.max = 1000;
*/
    posres->AnalyzePositionRes();

    std::vector<SFResults*> results = posres->GetResults();
    results[0]->Print(); //pol0
    results[1]->Print(); //pol3
    posres->Print();
    
    TGraphErrors* gPosRecoVsPosPol1 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPosResVsPosPol1 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualsPol1   = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kResidualGraph);
    TGraphErrors* gPosRecoDiffPol1 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::/*kPositionDiff*/kResidualGraph); 
    TH1D*         hPosRecoAllPol1  = (TH1D*)results[0]->GetObject(SFResultTypeObj::kPositionAllHist);

    TGraphErrors* gAttenuation_pol3 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosVsMLRGraph);
    
    TGraphErrors* gPosRecoVsPosPol3 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPosResVsPosPol3 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualsPol3   = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kResidualGraph);
    TH1D*         hPosRecoAllPol3  = (TH1D*)results[1]->GetObject(SFResultTypeObj::kPositionAllHist);

    int col_pol1 = kAzure - 3;
    int col_pol3 = kPink - 8;
    
    //------
    double f                 = 2 * sqrt(2 * log(2));
    double pos_all_pol3      = hPosRecoAllPol3->GetFunction("fun_gauss_pol3")->GetParameter(1);
    double pos_all_pol3_err  = hPosRecoAllPol3->GetFunction("fun_gauss_pol3")->GetParError(1);
    double fwhm_all_pol3     = hPosRecoAllPol3->GetFunction("fun_gauss_pol3")->GetParameter(2) * f;
    double fwhm_all_pol3_err = hPosRecoAllPol3->GetFunction("fun_gauss_pol3")->GetParError(2) * f;
    
    TGraphErrors *gPosResAllPol3 = new TGraphErrors(1);
    gPosResAllPol3->SetMarkerStyle(5);
    gPosResAllPol3->SetMarkerColor(col_pol3 + 1);
    gPosResAllPol3->SetLineColor(col_pol3 + 1);
    
    gPosResAllPol3->SetPoint(0, 100, fwhm_all_pol3);
    gPosResAllPol3->SetPointError(0, SFTools::GetPosError(colimator, testBench), fwhm_all_pol3_err);
    
    double pos_all_pol1      = hPosRecoAllPol1->GetFunction("fun_gauss_pol1")->GetParameter(1);
    double pos_all_pol1_err  = hPosRecoAllPol1->GetFunction("fun_gauss_pol1")->GetParError(1);
    double fwhm_all_pol1     = hPosRecoAllPol1->GetFunction("fun_gauss_pol1")->GetParameter(2) * f;
    double fwhm_all_pol1_err = hPosRecoAllPol1->GetFunction("fun_gauss_pol1")->GetParError(2) * f;
    
    TGraphErrors *gPosResAllPol1 = new TGraphErrors(1);
    gPosResAllPol1->SetMarkerStyle(5);
    gPosResAllPol1->SetMarkerColor(col_pol1 - 1);
    gPosResAllPol1->SetLineColor(col_pol1 - 1);
    
    gPosResAllPol1->SetPoint(0, 100, fwhm_all_pol1);
    gPosResAllPol1->SetPointError(0, SFTools::GetPosError(colimator, testBench), fwhm_all_pol1_err);
    //------
    
    double* posResAllPol3    = gPosResVsPosPol3->GetY();
    double* posResAllPol3Err = gPosResVsPosPol3->GetEY();
    double* posRecoPol3      = gPosRecoVsPosPol3->GetY();
    double* posRecoPol3Err   = gPosRecoVsPosPol3->GetEY();
    
    double* posResAllPol1    = gPosResVsPosPol1->GetY();
    double* posResAllPol1Err = gPosResVsPosPol1->GetEY();
    double* posRecoPol1      = gPosRecoVsPosPol1->GetY();
    double* posRecoPol1Err   = gPosRecoVsPosPol1->GetEY();

    std::vector<TH1D*> spec         = posres->GetSpectra();
    std::vector<TH1D*> hPosRecoPol3 = posres->GetPositionRecoDist("pol3");
    std::vector<TH1D*> hPosRecoPol1 = posres->GetPositionRecoDist("pol1");

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

    //----- position distributions & spectra
    
    TCanvas* can_posreco_dist_pol3 = new TCanvas("pr_posreco_dist_pol3", "pr_posreco_dist_pol3", 2000, 1200);
    can_posreco_dist_pol3->DivideSquare(npoints + 1);

    TCanvas* can_posreco_dist_pol1 = new TCanvas("pr_posreco_dist_pol1", "pr_posreco_dist_pol1", 2000, 1200);
    can_posreco_dist_pol1->DivideSquare(npoints + 1);

    TCanvas* can_spec = new TCanvas("pr_spec", "pr_spec", 2000, 1200);
    can_spec->DivideSquare(npoints);
    
    double min_yaxis = 0.;
    double max_yaxis =  SFTools::FindMaxYaxis(spec[4]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(spec[0]);

    TString fun_name;
    
    for (int i = 0; i < npoints; i++)
    {
        can_posreco_dist_pol3->cd(i + 1);
        gPad->SetGrid(1, 1);
        hPosRecoPol3[i]->Draw();
        text.DrawLatex(0.2, 0.8, Form("#mu = (%.2f +/- %.2f) mm", posRecoPol3[i], posRecoPol3Err[i]));
        text.DrawLatex(0.2, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", posResAllPol3[i], posResAllPol3Err[i]));

        can_posreco_dist_pol1->cd(i + 1);
        gPad->SetGrid(1, 1);
        hPosRecoPol1[i]->Draw();
        text.DrawLatex(0.2, 0.8, Form("#mu = (%.2f +/- %.2f) mm", posRecoPol1[i], posRecoPol1Err[i]));
        text.DrawLatex(0.2, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", posResAllPol1[i], posResAllPol1Err[i]));
        
        can_spec->cd(i + 1);
        gPad->SetGrid(1, 1);
        spec[i]->SetStats(false);
        //ctx.configureFromJson("hSpecAv");
        //ctx.print();
        //spec[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //spec[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        //line.DrawLine(xmin[i], ctx.y.min, xmin[i], ctx.y.max);
        //line.DrawLine(xmax[i], ctx.y.min, xmax[i], ctx.y.max);
        spec[i]->SetTitle(Form("Average PE Spectrum S%i %.1f mm", seriesNo, positions[i]));
        spec[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        spec[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis);
        spec[i]->Draw();
        
        line.DrawLine(xmin[i], min_yaxis, xmin[i], max_yaxis);
        line.DrawLine(xmax[i], min_yaxis, xmax[i], max_yaxis);
        
        fun_name = Form("f_S%i_pos%.1f_ID%i_PEAverage", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.65, 0.75, Form("#mu = %.2f +/- %.2f", 
                       spec[i]->GetFunction(fun_name)->GetParameter(1),
                       spec[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.65, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       spec[i]->GetFunction(fun_name)->GetParameter(2),
                       spec[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.65, 0.60, Form("#chi^{2}/NDF = %.3f",
                       spec[i]->GetFunction(fun_name)->GetChisquare() / 
                       spec[i]->GetFunction(fun_name)->GetNDF()));
    }
    
    can_posreco_dist_pol3->cd(npoints + 1);
    gPad->SetGrid(1, 1);
    gPad->SetFillColor(kGreen - 10);
    hPosRecoAllPol3->SetLineColor(kBlack);
    hPosRecoAllPol3->Draw();
    text.DrawLatex(0.2, 0.8, Form("#mu = (%.2f +/- %.2f) mm", pos_all_pol3, pos_all_pol3_err));
    text.DrawLatex(0.2, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", fwhm_all_pol3, fwhm_all_pol3_err));
    
    can_posreco_dist_pol1->cd(npoints + 1);
    gPad->SetGrid(1, 1);
    gPad->SetFillColor(kGreen - 10);
    hPosRecoAllPol1->SetLineColor(kBlack);
    hPosRecoAllPol1->Draw();
    text.DrawLatex(0.2, 0.8, Form("#mu = (%.2f +/- %.2f) mm", pos_all_pol1, pos_all_pol1_err));
    text.DrawLatex(0.2, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", fwhm_all_pol1, fwhm_all_pol1_err));

    //----- position reconstruction
    
    TCanvas* can_posreco = new TCanvas("pr_posreco", "pr_posreco", 1500, 700);
    can_posreco->Divide(2, 1);
    
    can_posreco->cd(1);
    TPad* pad_posreco_pol1 = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0);
    TPad* pad_res_pol1     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);

    pad_posreco_pol1->Draw();
    pad_posreco_pol1->cd();
    pad_posreco_pol1->SetGrid(1, 1);
    gPosRecoVsPosPol1->Draw("AP");
    // gPosRecoVsPosPol3->GetXaxis()->SetLabelSize(0.035);
    // gPosRecoVsPosPol3->GetYaxis()->SetLabelSize(0.035);
    // gPosRecoVsPosPol3->GetYaxis()->SetTitleOffset(1.4);

    TF1* funpol1 = gPosRecoVsPosPol1->GetFunction("funpol1_p1");

    text.SetTextSize(0.03);
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.2f +/- %.2f", 
                   funpol1->GetParameter(0), funpol1->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.2f +/- %.2f", 
                   funpol1->GetParameter(1), funpol1->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   funpol1->GetChisquare() / funpol1->GetNDF()));

    can_posreco->cd(1);
    pad_res_pol1->Draw();
    pad_res_pol1->cd();
    pad_res_pol1->SetGrid(1, 1);
    gResidualsPol1->Draw("AP");
    
    can_posreco->cd(2);
    TPad* pad_posreco_pol3 = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0);
    TPad* pad_res_pol3     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);

    pad_posreco_pol3->Draw();
    pad_posreco_pol3->cd();
    pad_posreco_pol3->SetGrid(1, 1);
    gPosRecoVsPosPol3->Draw("AP");
    // gPosRecoVsPosPol3->GetXaxis()->SetLabelSize(0.035);
    // gPosRecoVsPosPol3->GetYaxis()->SetLabelSize(0.035);
    // gPosRecoVsPosPol3->GetYaxis()->SetTitleOffset(1.4);

    funpol1 = gPosRecoVsPosPol3->GetFunction("funpol1_p3");

    text.SetTextSize(0.03);
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.2f +/- %.2f", 
                   funpol1->GetParameter(0), funpol1->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.2f +/- %.2f", 
                   funpol1->GetParameter(1), funpol1->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   funpol1->GetChisquare() / funpol1->GetNDF()));

    can_posreco->cd(2);
    pad_res_pol3->Draw();
    pad_res_pol3->cd();
    pad_res_pol3->SetGrid(1, 1);
    gResidualsPol3->Draw("AP");

    //----- position resolution
    
    TCanvas* can_posres = new TCanvas("pr_posres", "pr_posres", 1500, 700);
    can_posres->Divide(2, 1);
    
    can_posres->cd(1);
    gPad->SetGrid(1, 1);
    
    TH1D *htemp_pol1 = new TH1D("htemp_pol1", gPosResVsPosPol1->GetTitle(), 105, 0, 105);
    htemp_pol1->GetXaxis()->SetTitle(gPosResVsPosPol1->GetXaxis()->GetTitle());
    htemp_pol1->GetYaxis()->SetTitle(gPosResVsPosPol1->GetYaxis()->GetTitle());
    htemp_pol1->SetStats(0);
    
    text.SetTextSize(0.04);
    htemp_pol1->Draw();
    gPosResVsPosPol1->SetMarkerColor(col_pol1);
    gPosResVsPosPol1->SetLineColor(col_pol1);
    gPosResVsPosPol1->Draw("P");
    // gPosResVsPosPol3->GetXaxis()->SetLabelSize(0.035);
    // gPosResVsPosPol3->GetYaxis()->SetLabelSize(0.035);
    // gPosResVsPosPol3->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.2, 0.8, Form("PR = (%.2f +/- %.2f) mm",
                   results[0]->GetValue(SFResultTypeNum::kPositionRes),
                   results[0]->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    gPosResAllPol1->Draw("P");

    double  ymin    = TMath::Min(TMath::MinElement(npoints, gPosResVsPosPol1->GetY()), gPosResAllPol1->GetPointY(0));
    double  ymax    = TMath::Max(TMath::MaxElement(npoints, gPosResVsPosPol3->GetY()), gPosResAllPol1->GetPointY(0));
    
    htemp_pol1->GetYaxis()->SetRangeUser(ymin - 0.05 * fabs(ymin), ymax + 0.05 * fabs(ymax));
    
    text.SetNDC(0);
    text.SetTextSize(0.03);
    text.DrawLatex(75, fwhm_all_pol1 + 0.2, Form("PR = (%.2f +/- %.2f) mm", fwhm_all_pol1, fwhm_all_pol1_err));
    text.SetNDC(1);
    text.SetTextSize(0.04);
    
    can_posres->cd(2);
    gPad->SetGrid(1, 1);

    TH1D *htemp_pol3 = new TH1D("htemp_pol3", gPosResVsPosPol3->GetTitle(), 105, 0, 105);
    htemp_pol3->GetXaxis()->SetTitle(gPosResVsPosPol3->GetXaxis()->GetTitle());
    htemp_pol3->GetYaxis()->SetTitle(gPosResVsPosPol3->GetYaxis()->GetTitle());
    htemp_pol3->SetStats(0);
    
    text.SetTextSize(0.04);
    htemp_pol3->Draw();
    gPosResVsPosPol3->SetMarkerColor(col_pol3);
    gPosResVsPosPol3->SetLineColor(col_pol3);
    gPosResVsPosPol3->Draw("P");
    // gPosResVsPosPol3->GetXaxis()->SetLabelSize(0.035);
    // gPosResVsPosPol3->GetYaxis()->SetLabelSize(0.035);
    // gPosResVsPosPol3->GetYaxis()->SetTitleOffset(1.4);
    text.DrawLatex(0.2, 0.8, Form("PR = (%.2f +/- %.2f) mm",
                   results[1]->GetValue(SFResultTypeNum::kPositionRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    gPosResAllPol3->Draw("P");

    ymin = TMath::Min(TMath::MinElement(npoints, gPosResVsPosPol3->GetY()), gPosResAllPol3->GetPointY(0));
    ymax = TMath::Max(TMath::MaxElement(npoints, gPosResVsPosPol3->GetY()), gPosResAllPol3->GetPointY(0));
    
    htemp_pol3->GetYaxis()->SetRangeUser(ymin - 0.05 * fabs(ymin), ymax + 0.05 * fabs(ymax));
    
    text.SetNDC(0);
    text.SetTextSize(0.03);
    text.DrawLatex(75, fwhm_all_pol3 + 0.2, Form("PR = (%.2f +/- %.2f) mm", fwhm_all_pol3, fwhm_all_pol3_err));
    text.SetNDC(1);
    text.SetTextSize(0.04);
    
    //----- reconstruction position difference
   
    TCanvas *can_diff = new TCanvas("pr_diff_pol1", "pr_diff_pol1", 700, 500);
    gPad->SetGrid(1, 1);
    
    gPosRecoDiffPol1->Draw("AP");
    //-----
    
    //----- MLR curve
    
    TCanvas* can_att = new TCanvas("pr_att", "pr_att", 1500, 700);
    can_att->Divide(2,1);
    
    can_att->cd(1);
    gPad->SetGrid(1, 1);
    
    TGraphErrors *gAttenaution_pol1 = (TGraphErrors*)gAttenuation_pol3->Clone();
    gAttenaution_pol1->Draw("AP");
    gAttenaution_pol1->GetFunction("funpol3")->Delete();
    TF1* fPol1 = (TF1*)gAttenaution_pol1->GetFunction("funpol1");
    fPol1->SetLineColor(col_pol1);
    
    text.SetTextSize(0.025);
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.4e +/- %.4e", 
                   fPol1->GetParameter(0), fPol1->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.4e +/- %.4e", 
                   fPol1->GetParameter(1), fPol1->GetParError(1)));
    text.DrawLatex(0.2, 0.55, Form("#chi^{2}/NDF = %.2f", fPol1->GetChisquare() / 
                   fPol1->GetNDF()));
    
    can_att->cd(2);
    gPad->SetGrid(1, 1);
    
    gAttenuation_pol3->Draw("AP");
    TF1* fPol3 = (TF1*)gAttenuation_pol3->GetFunction("funpol3");
    gAttenuation_pol3->GetFunction("funpol1")->Delete();
    fPol3->SetLineColor(col_pol3);
    
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
    text.DrawLatex(0.2, 0.55, Form("#chi^{2}/NDF = %.2f", fPol3->GetChisquare() / 
                   fPol3->GetNDF()));

    //----- reconstruction function
    
    TCanvas* can_fun = new TCanvas("pr_fun", "pr_fun", 700, 500);

    TH1F* hPol3 = new TH1F("hPol3", "hPol3", 100, -1, 1);
    gPad->SetGrid(1, 1);
    hPol3->Draw("func");
    fPol3->SetLineColor(col_pol3);
    fPol1->SetLineColor(col_pol1);
    fPol3->Draw("same");
    fPol1->Draw("same");
    hPol3->SetStats(0);
    hPol3->GetXaxis()->SetRangeUser(-1, 1);
    hPol3->GetYaxis()->SetRangeUser(fPol3->Eval(-1), fPol3->Eval(1));
    hPol3->GetYaxis()->SetTitle("source position [mm]");
    hPol3->GetXaxis()->SetTitle("ln(M_{LR})");
    hPol3->SetTitle(Form("Calibration Curve for Position Reconstruction S%i", seriesNo));

    TLegend *leg = new TLegend(0.2, 0.8, 0.3, 0.9);
    leg->AddEntry(fPol1, "pol1 curve", "L");
    leg->AddEntry(fPol3, "pol3 curve", "L");
    leg->Draw();
    
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

    can_posreco_dist_pol3->Write();
    can_posreco_dist_pol1->Write();
    can_posreco->Write();
    can_spec->Write();
    can_posres->Write();
    can_att->Write();
    can_fun->Write();
    can_diff->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "POSITION_RESOLUTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES_POL3, "
                    "POSITION_RES_POL3_ERR, POSITION_RES_POL3_ALL, POSITION_RES_POL3_ALL_ERR, "
                    "POSITION_RES_POL1, POSITION_RES_POL1_ERR, POSITION_RES_POL1_ALL, "
                    "POSITION_RES_POL1_ALL_ERR) VALUES(%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f)", 
                    table.Data(), seriesNo, fname_full.Data(), 
                    results[1]->GetValue(SFResultTypeNum::kPositionRes),
                    results[1]->GetUncertainty(SFResultTypeNum::kPositionRes),
                    fwhm_all_pol3, fwhm_all_pol3_err,
                    results[0]->GetValue(SFResultTypeNum::kPositionRes),
                    results[0]->GetUncertainty(SFResultTypeNum::kPositionRes),
                    fwhm_all_pol1, fwhm_all_pol1_err);

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

    for (auto h : hPosRecoPol3)
        delete h;
    
    for (auto h : hPosRecoPol1)
        delete h;
    
    for (auto h : spec)
        delete h;
    
    for (auto pf : peakFinAv)
        delete pf;
    
    delete data;
    delete posres;

    return 0;
}
