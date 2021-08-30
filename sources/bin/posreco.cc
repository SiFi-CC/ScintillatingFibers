// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              posreco.cc               *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFPositionReco.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMathBase.h>

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
        std::cout << "to run type: ./reco seriesNo";
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
        std::cerr << "##### Exception in posreco.cc!" << std::endl;
        return 1;
    }

    data->Print();

    TString desc      = data->GetDescription();
    TString colimator = data->GetCollimator();
    TString testBench = data->GetTestBench();
    int     npoints   = data->GetNpoints();
    double  fiberLen  = data->GetFiberLength();

    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in posreco.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cout << "Description: " << desc << std::endl;
        return 1;
    }

    SFPositionReco* reco;
    try
    {
        reco = new SFPositionReco(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in posreco.cc!" << std::endl;
        return 1;
    }

    reco->Print();

    reco->CalculateMLR();
    reco->CalculatePosRecoCoefficients();
    reco->PositionReco();

    std::vector<SFResults*> results = reco->GetResults();
    results[0]->Print(); // experimental
    results[1]->Print(); // corrected

    //----- getting position reconstruction results
    TGraphErrors* gMLRExp        = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors* gPRecoExp      = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPResExp       = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualExp   = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kResidualGraph);
    TGraphErrors* gPosDiffExp    = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::/*kPositionDiff*/kResidualGraph);
    TH1D*         hPosRecoExpAll = (TH1D*)results[0]->GetObject(SFResultTypeObj::kPositionAllHist);
    
    std::vector<TH1D*> hPosDistExp = reco->GetPositionDistributions("experimental");
    
    TGraphErrors* gMLRCorr        = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors* gPRecoCorr      = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPResCorr       = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualCorr   = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kResidualGraph);
    TGraphErrors* gACoeff         = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAGraph);
    TGraphErrors* gRevMLR         = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosVsMLRGraph);
    TGraphErrors* gPosDiffCorr    = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::/*kPositionDiff*/kResidualGraph);
    TH1D*         hPosRecoCorrAll = (TH1D*)results[1]->GetObject(SFResultTypeObj::kPositionAllHist);
    
    std::vector<TH1D*> hPosDistCorr   = reco->GetPositionDistributions("corrected");
    std::vector<TH1D*> hPosUncertCorr = reco->GetErrorDistributions(); 
    //-----
    
    //----- 
    double f            = 2 * sqrt(2 * log(2));
    double fwhm_exp_all     = hPosRecoExpAll->GetFunction("fun_gauss_pol1")->GetParameter(2) * f;
    double fwhm_exp_all_err = hPosRecoExpAll->GetFunction("fun_gauss_pol1")->GetParError(2) * f;
    
    TGraphErrors *gPosResExpAll = new TGraphErrors(1);
    gPosResExpAll->SetMarkerStyle(5);
    gPosResExpAll->SetMarkerColor(kBlack);
    gPosResExpAll->SetLineColor(kBlack);
    
    gPosResExpAll->SetPoint(0, 100, fwhm_exp_all);
    gPosResExpAll->SetPointError(0, SFTools::GetPosError(colimator, testBench), fwhm_exp_all_err);

    double pos_corr_all      =  hPosRecoCorrAll->GetFunction("fun_gauss")->GetParameter(1);
    double pos_corr_all_err  =  hPosRecoCorrAll->GetFunction("fun_gauss")->GetParError(1);
    double fwhm_corr_all     = hPosRecoCorrAll->GetFunction("fun_gauss")->GetParameter(2) * f;
    double fwhm_corr_all_err = hPosRecoCorrAll->GetFunction("fun_gauss")->GetParError(2) * f;
    
    TGraphErrors *gPosResCorrAll = new TGraphErrors(1);
    gPosResCorrAll->SetMarkerStyle(5);
    gPosResCorrAll->SetMarkerColor(kBlack);
    gPosResCorrAll->SetLineColor(kBlack);
    
    gPosResCorrAll->SetPoint(0, 100, fwhm_corr_all);
    gPosResCorrAll->SetPointError(0, SFTools::GetPosError(colimator, testBench), fwhm_corr_all_err);
    //-----

    int col_exp  = kTeal + 4;
    int col_corr = kAzure - 5;
    
    TLatex text;
    text.SetNDC(true);
    
    //----- drawing MLR
    TCanvas* can_mlr = new TCanvas("preco_mlr", "preco_mlr", 700, 500);
    gPad->SetGrid(1, 1);
    
    gMLRExp->SetMarkerColor(col_exp);
    gMLRExp->SetLineColor(col_exp);
    gMLRExp->SetMarkerStyle(8);
    gMLRExp->GetFunction("fpol1")->SetLineColor(col_exp);
    gMLRExp->GetFunction("fpol3")->Delete();
    gMLRExp->Draw("AP");
    
    gMLRCorr->SetMarkerColor(col_corr);
    gMLRCorr->SetLineColor(col_corr);
    gMLRCorr->GetFunction("fpol1")->SetLineColor(col_corr);
    gMLRCorr->Draw("P");
    
    text.SetTextColor(col_exp);
    text.DrawLatex(0.15, 0.80, Form("Slope_{exp} = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kMLRSlope),
                   results[0]->GetUncertainty(SFResultTypeNum::kMLRSlope)));
    text.DrawLatex(0.15, 0.75, Form("Offsef_{exp} = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kMLROffset),
                   results[0]->GetUncertainty(SFResultTypeNum::kMLROffset)));
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.70, Form("Slope_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kMLRSlope),
                   results[1]->GetUncertainty(SFResultTypeNum::kMLRSlope)));
    text.DrawLatex(0.15, 0.65, Form("Offset_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kMLROffset),
                   results[1]->GetUncertainty(SFResultTypeNum::kMLROffset)));
    
    double min = TMath::Min(TMath::MinElement(npoints, gMLRExp->GetY()),
                            TMath::MinElement(npoints, gMLRCorr->GetY()));
    double max = TMath::Max(TMath::MaxElement(npoints, gMLRExp->GetY()),
                            TMath::MaxElement(npoints, gMLRCorr->GetY()));
    gMLRExp->GetYaxis()->SetRangeUser(min - 0.1 * fabs(min), max + 0.1 * max);
    
    TLegend* leg_1 = new TLegend(0.5, 0.15, 0.9, 0.25);
    leg_1->AddEntry(gMLRExp, "M_{LR} experimental", "PE");
    leg_1->AddEntry(gMLRCorr, "M_{LR} corrected", "PE");
    leg_1->Draw();
    //-----
    
    //----- drawing A graph and reversed MLR
    TCanvas *can_a = new TCanvas("preco_a", "preco_a", 1200, 600);
    can_a->Divide(2, 1);
    
    can_a->cd(1);
    gPad->SetGrid(1, 1);
    gACoeff->SetMarkerColor(col_corr);
    gACoeff->SetLineColor(col_corr);
    gACoeff->Draw("AP");
    
    gACoeff->GetFunction("fpol0")->SetLineColor(col_corr);
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.75, Form("A = %.3f +/- %.3f",
                   gACoeff->GetFunction("fpol0")->GetParameter(0),
                   gACoeff->GetFunction("fpol0")->GetParError(0)));
    text.DrawLatex(0.15, 0.70, Form("B = %.3f", fiberLen / 2.));
    
    can_a->cd(2);
    gPad->SetGrid(1, 1);
    gRevMLR->SetMarkerColor(col_corr);
    gRevMLR->SetLineColor(col_corr);
    gRevMLR->Draw("AP");
    
    gRevMLR->GetFunction("fpol1")->SetLineColor(col_corr);
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.2, 0.75, Form("A = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kACoeff),
                   results[1]->GetUncertainty(SFResultTypeNum::kACoeff)));
    text.DrawLatex(0.2, 0.70, Form("B = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kBCoeff),
                   results[1]->GetUncertainty(SFResultTypeNum::kBCoeff)));
    //-----
    
    //----- drawing position resolution
    TCanvas *can_pos_res = new TCanvas("preco_pos_res", "preco_pos_res", 1200, 600);
    can_pos_res->Divide(2, 1);
    
    TH1D* htemp_exp = new TH1D("htemo_exp", gPResExp->GetTitle(), 105, 0, 105);
    htemp_exp->SetStats(0);
    htemp_exp->GetXaxis()->SetTitle(gPResExp->GetXaxis()->GetTitle());
    htemp_exp->GetYaxis()->SetTitle(gPResExp->GetYaxis()->GetTitle());
    
    can_pos_res->cd(1);
    gPad->SetGrid(1, 1);
    htemp_exp->Draw();
    
    gPResExp->SetMarkerColor(col_exp);
    gPResExp->SetLineColor(col_exp);
    gPResExp->SetMarkerStyle(8);
    gPResExp->Draw("P");
    
    text.SetTextColor(col_exp);
    text.SetTextSize(0.04);
    text.DrawLatex(0.15, 0.75, Form("PR_{exp} = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kPositionRes),
                   results[0]->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    gPosResExpAll->SetLineColor(col_exp - 1);
    gPosResExpAll->SetMarkerColor(col_exp - 1);
    gPosResExpAll->Draw("P");
    
    double  ymin    = TMath::Min(TMath::MinElement(npoints, gPResExp->GetY()), gPosResExpAll->GetPointY(0));
    double  ymax    = TMath::Max(TMath::MaxElement(npoints, gPResExp->GetY()), gPosResExpAll->GetPointY(0));
    
    htemp_exp->GetYaxis()->SetRangeUser(ymin - 0.05 * fabs(ymin), ymax + 0.05 * fabs(ymax));
    
    text.SetNDC(0);
    text.SetTextSize(0.03);
    text.SetTextColor(col_exp - 1);
    text.DrawLatex(75, fwhm_exp_all + 0.2, Form("PR = (%.2f +/- %.2f) mm", fwhm_exp_all, fwhm_exp_all_err));
    text.SetNDC(1);
    text.SetTextSize(0.04);
    
    TH1D* htemp_corr = new TH1D("htemo_corr", gPResCorr->GetTitle(), 105, 0, 105);
    htemp_corr->SetStats(0);
    htemp_corr->GetXaxis()->SetTitle(gPResCorr->GetXaxis()->GetTitle());
    htemp_corr->GetYaxis()->SetTitle(gPResCorr->GetYaxis()->GetTitle());
    
    can_pos_res->cd(2);
    gPad->SetGrid(1, 1);
    htemp_exp->Draw();
    gPResCorr->SetMarkerColor(col_corr);
    gPResCorr->SetLineColor(col_corr);
    gPResCorr->Draw("P");
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.75, Form("PR_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kPositionRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    gPosResCorrAll->SetLineColor(col_corr - 1);
    gPosResCorrAll->SetMarkerColor(col_corr - 1);
    gPosResCorrAll->Draw("P");
    
    //ymin    = TMath::Min(TMath::MinElement(npoints, gPResCorr->GetY()), gPosResCorrAll->GetPointY(0));
    //ymax    = TMath::Max(TMath::MaxElement(npoints, gPResCorr->GetY()), gPosResCorrAll->GetPointY(0));
    
    htemp_exp->GetYaxis()->SetRangeUser(ymin - 0.05 * fabs(ymin), ymax + 0.05 * fabs(ymax));
    
    text.SetNDC(0);
    text.SetTextSize(0.03);
    text.SetTextColor(col_corr - 1);
    text.DrawLatex(75, fwhm_corr_all + 0.2, Form("PR = (%.2f +/- %.2f) mm", fwhm_corr_all, fwhm_corr_all_err));
    text.SetNDC(1);
    text.SetTextSize(0.04);
    //-----
    
    //----- drawing position reconstructed
    TCanvas* can_pos_reco = new TCanvas("preco_pos_reco", "preco_pos_reco", 1200, 600);
    can_pos_reco->Divide(2, 1);
    
    TPad* pad_posreco_1 = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0, 0);
    TPad* pad_res_1     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);
    
    TLine line;
    line.SetLineStyle(2);
    line.SetLineColor(kRed);
    line.SetLineWidth(2);
    
    text.SetTextColor(kBlack);
    text.SetTextSize(0.030);
    
    can_pos_reco->cd(1);
    pad_posreco_1->Draw();
    pad_posreco_1->cd();
    gPad->SetGrid(1, 1);
    gPRecoExp->SetMarkerColor(col_exp);
    gPRecoExp->SetLineColor(col_exp);
    gPRecoExp->SetMarkerStyle(8);
    gPRecoExp->GetFunction("funpol1_p1")->SetLineColor(col_exp);
    gPRecoExp->Draw("AP");
    
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.3f +/- %.3f",  
                   gPRecoExp->GetFunction("funpol1_p1")->GetParameter(0),
                   gPRecoExp->GetFunction("funpol1_p1")->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.3f +/- %.3f",  
                   gPRecoExp->GetFunction("funpol1_p1")->GetParameter(1),
                   gPRecoExp->GetFunction("funpol1_p1")->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.3f",  
                   gPRecoExp->GetFunction("funpol1_p1")->GetChisquare() /
                   gPRecoExp->GetFunction("funpol1_p1")->GetNDF()));
    
    can_pos_reco->cd(1);
    pad_res_1->Draw();
    pad_res_1->cd();
    gPad->SetGrid(1, 1);
    gResidualExp->SetMarkerColor(col_exp);
    gResidualExp->SetLineColor(col_exp);
    gResidualExp->SetMarkerStyle(8);
    gResidualExp->Draw("AP");
    line.DrawLine(0, gResidualExp->GetPointX(0),
                  100, gResidualExp->GetPointX(npoints - 1));

    TPad* pad_posreco_2 = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0);
    TPad* pad_res_2     = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0, 0);
    
    can_pos_reco->cd(2);
    pad_posreco_2->Draw();
    pad_posreco_2->cd();
    gPad->SetGrid(1, 1);
    gPRecoCorr->SetMarkerColor(col_corr);
    gPRecoCorr->SetLineColor(col_corr);
    gPRecoCorr->GetFunction("funpol1")->SetLineColor(col_corr);
    gPRecoCorr->Draw("AP");
    
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.3f +/- %.3f",  
                   gPRecoCorr->GetFunction("funpol1")->GetParameter(0),
                   gPRecoCorr->GetFunction("funpol1")->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.3f +/- %.3f",  
                   gPRecoCorr->GetFunction("funpol1")->GetParameter(1),
                   gPRecoCorr->GetFunction("funpol1")->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.3f",  
                   gPRecoCorr->GetFunction("funpol1")->GetChisquare() /
                   gPRecoCorr->GetFunction("funpol1")->GetNDF()));
    
    can_pos_reco->cd(2);
    pad_res_2->Draw();
    pad_res_2->cd();
    gPad->SetGrid(1, 1);
    gResidualCorr->SetMarkerColor(col_corr);
    gResidualCorr->SetLineColor(col_corr);
    gResidualCorr->Draw("AP");
    line.DrawLine(0, gResidualCorr->GetPointX(0), 
                  100, gResidualCorr->GetPointX(npoints - 1));
    //----- 
    
    //----- drawing position diff
    
    TCanvas *can_diff = new TCanvas("preco_pos_diff", "preco_pos_diff", 700, 500);
    gPad->SetGrid(1, 1);
    
    gPosDiffCorr->GetYaxis()->SetRangeUser(-50, 50);
    gPosDiffCorr->SetMarkerColor(col_corr);
    gPosDiffCorr->SetLineColor(col_corr);
    gPosDiffExp->SetMarkerColor(col_exp);
    gPosDiffExp->SetLineColor(col_exp);
    
    gPosDiffCorr->Draw("AP");
    gPosDiffExp->Draw("P");
    
    TLegend *leg_2 = new TLegend();
    leg_2->AddEntry(gPosDiffCorr, "corrected data", "PEL");
    leg_2->AddEntry(gPosDiffExp, "experimental data", "PEL");
    leg_2->Draw();

    //-----
    
    //----- drawing position distributions
    TCanvas *can_pos_dist = new TCanvas("preco_pos_dist", "preco_pos_dist", 2000, 1200);
    can_pos_dist->DivideSquare(npoints + 1);
    
    TCanvas *can_pos_uncert = new TCanvas("preco_pos_uncert", "preco_pos_uncert", 2000, 1200);
    can_pos_uncert->DivideSquare(npoints);
    
    text.SetTextColor(kBlack);
    text.SetTextSize(0.040);
    
    for (int i = 0; i < npoints; i++)
    {
        can_pos_dist->cd(i + 1);
        gPad->SetGrid(1, 1);
        
        //hPosDistCorr[i]->SetStats(0);
        hPosDistCorr[i]->GetYaxis()->SetTitle("counts");
        hPosDistCorr[i]->GetXaxis()->SetTitle("reconstructed position [mm]");
        hPosDistCorr[i]->Draw();
        
        text.DrawLatex(0.2, 0.80, Form("#mu = (%.3f +/- %.3f) mm",
                       gPRecoCorr->GetPointY(i),
                       gPRecoCorr->GetErrorY(i)));
        text.DrawLatex(0.2, 0.75, Form("FWHM = (%.3f +/- %.3f) mm",
                       gPResCorr->GetPointY(i),
                       gPResCorr->GetErrorY(i)));
        
        can_pos_uncert->cd(i + 1);
        gPad->SetGrid(1, 1);
        
        hPosUncertCorr[i]->GetYaxis()->SetTitle("counts");
        hPosUncertCorr[i]->GetXaxis()->SetTitle("reconstructed position uncertainty [mm]");
        hPosUncertCorr[i]->Draw();
    }
    
    can_pos_dist->cd(npoints + 1);
    gPad->SetGrid(1, 1);    gPad->SetFillColor(kGreen - 10);
    hPosRecoCorrAll->SetLineColor(kBlack);
    hPosRecoCorrAll->Draw();
    text.DrawLatex(0.2, 0.8, Form("#mu = (%.2f +/- %.2f) mm", pos_corr_all, pos_corr_all_err));
    text.DrawLatex(0.2, 0.7, Form("FWHM = (%.2f +/- %.2f) mm", fwhm_corr_all, fwhm_corr_all_err));
    
    //-----
    
    //----- saving
    TString fname       = Form("posreco_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in posreco.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_mlr->Write();
    can_a->Write();
    can_pos_res->Write();
    can_pos_reco->Write();
    can_diff->Write();
    can_pos_dist->Write();
    can_pos_uncert->Write();
    file->Close();

    //-----writing results to the data base
    TString table = "POSITION_RECONSTRUCTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, A_COEFF, "
                         "A_COEFF_ERR, B_COEFF, MLR_SLOPE, MLR_SLOPE_ERR, MLR_OFFSET, " 
                         "MLR_OFFSET_ERR, MLR_SLOPE_EXP, MLR_SLOPE_EXP_ERR, MLR_OFFSET_EXP, "
                         "MLR_OFFSET_EXP_ERR, POSITION_RES, POSITION_RES_ERR, POSITION_RES_ALL, "
                         "POSITION_RES_ALL_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, "
                         "%f, %f, %f, %f,%f, %f, %f)", table.Data(), seriesNo, fname_full.Data(),
                         results[1]->GetValue(SFResultTypeNum::kACoeff),
                         results[1]->GetUncertainty(SFResultTypeNum::kACoeff),
                         results[1]->GetValue(SFResultTypeNum::kBCoeff),
                         results[1]->GetValue(SFResultTypeNum::kMLRSlope),
                         results[1]->GetUncertainty(SFResultTypeNum::kMLRSlope),
                         results[1]->GetValue(SFResultTypeNum::kMLROffset),
                         results[1]->GetUncertainty(SFResultTypeNum::kMLROffset),
                         results[0]->GetValue(SFResultTypeNum::kMLRSlope),
                         results[0]->GetUncertainty(SFResultTypeNum::kMLRSlope),
                         results[0]->GetValue(SFResultTypeNum::kMLROffset),
                         results[0]->GetUncertainty(SFResultTypeNum::kMLROffset),
                         results[1]->GetValue(SFResultTypeNum::kPositionRes),
                         results[1]->GetUncertainty(SFResultTypeNum::kPositionRes),
                         fwhm_corr_all, fwhm_corr_all_err);
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;
    
    do
    {
        std::cout << "----- posreco writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
        
        if (stat) break;
        
        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- posreco writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    for (auto h : hPosDistCorr)
        delete h;
    
    delete data;
    delete reco;
    
    return 0;
}
