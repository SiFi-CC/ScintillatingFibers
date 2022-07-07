// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             energyreco.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFEnergyReco.hh"
#include "common_options.h"

#include <TPave.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMathBase.h>
#include <TMath.h>

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
        std::cout << "to run type: ./energyreco seriesNo";
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
        std::cerr << "##### Exception in energyreco.cc!" << std::endl;
        return 1;
    }

    data->Print();

    TString             desc      = data->GetDescription();
    int                 npoints   = data->GetNpoints();
    std::vector<double> positions = data->GetPositions();

    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in energyreco.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cout << "Description: " << desc << std::endl;
        return 1;
    }

    SFEnergyReco* reco = nullptr;
    try
    {
        reco = new SFEnergyReco(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in energyreco.cc!" << std::endl;
        return 1;
    }

    reco->Print();
    reco->CalculateAlpha();
    reco->EnergyReco();
    reco->EnergyRecoByEvent();

    std::vector<SFResults*> results = reco->GetResults();
    results[0]->Print(); // experimental
    results[1]->Print(); // corrected

    //----- getting energy reconstruction results
    TGraphErrors* gEnReco = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kEnergyRecoGraph);
    TGraphErrors* gAlpha = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kAlphaGraph);
    TGraphErrors* gEnRecoSpec = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kEnergyRecoSpecGraph);
    TGraphErrors* gEnRes = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kEnergyResGraph);
    TF1*  fEnReco = (TF1*)results[0]->GetObject(SFResultTypeObj::kEnergyRecoFun);
    
    TGraphErrors* gEnRecoCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyRecoGraph);
    TGraphErrors* gAlphaCorr  = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAlphaGraph);
    TGraphErrors* gEnRecoSpecCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyRecoSpecGraph);
    TGraphErrors* gEnResCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyResGraph);
    TF1*  fEnRecoCor  = (TF1*)results[1]->GetObject(SFResultTypeObj::kEnergyRecoFun);
    
    std::vector <TH1D*> hEnReco       = reco->GetEnergySpectra("experimental");
    std::vector <TH1D*> hEnRecoCorr   = reco->GetEnergySpectra("corrected");
    std::vector <TH1D*> hEnRecoUncert = reco->GetErrorDistributions(); 
    //-----
    
    //----- collective reconstructed energy histograms
    TH1D* hEnRecoAllSpec = (TH1D*)results[0]->GetObject(SFResultTypeObj::kEnergyAllHist);
    TH1D* hEnRecoAllCorrSpec = (TH1D*)results[1]->GetObject(SFResultTypeObj::kEnergyAllHist);
    
    SFPeakFinder *pf_ereco_all = new SFPeakFinder(hEnRecoAllSpec, false, true);
    pf_ereco_all->FindPeakFit();
    SFResults *results_ereco_all = pf_ereco_all->GetResults();
    double eres_all = results_ereco_all->GetValue(SFResultTypeNum::kPeakSigma) / 
                      results_ereco_all->GetValue(SFResultTypeNum::kPeakPosition) * 100;
    double eres_all_err = eres_all * 
                          sqrt(pow((results_ereco_all->GetUncertainty(SFResultTypeNum::kPeakPosition) /
                          results_ereco_all->GetValue(SFResultTypeNum::kPeakPosition)), 2) + pow((results_ereco_all->GetUncertainty(SFResultTypeNum::kPeakSigma) /
                          results_ereco_all->GetValue(SFResultTypeNum::kPeakSigma)), 2));
    
    TGraphErrors *gEnResAll = new TGraphErrors(1);
    gEnResAll->SetPoint(0, 12, eres_all);
    gEnResAll->SetPointError(0, 0, eres_all_err);
    
    SFPeakFinder *pf_ereco_corr_all = new SFPeakFinder(hEnRecoAllCorrSpec, false, true);
    pf_ereco_corr_all->FindPeakFit();
    SFResults *results_ereco_corr_all = pf_ereco_corr_all->GetResults();
    double eres_corr_all = results_ereco_corr_all->GetValue(SFResultTypeNum::kPeakSigma) / 
                           results_ereco_corr_all->GetValue(SFResultTypeNum::kPeakPosition) * 100;
    double eres_corr_all_err = eres_corr_all *
                               sqrt(pow((results_ereco_corr_all->GetUncertainty(SFResultTypeNum::kPeakPosition) /
                               results_ereco_corr_all->GetValue(SFResultTypeNum::kPeakPosition)), 2) + pow((results_ereco_corr_all->GetUncertainty(SFResultTypeNum::kPeakSigma) /
                               results_ereco_corr_all->GetValue(SFResultTypeNum::kPeakSigma)), 2));
    
    TGraphErrors *gEnResCorrAll = new TGraphErrors(1);
    gEnResCorrAll->SetPoint(0, 12, eres_corr_all);
    gEnResCorrAll->SetPointError(0, 0, eres_corr_all_err);
    //-----

    int col_exp  = kMagenta -3;
    int col_corr = kAzure - 5;
    
    //----- drawing energy reconstruction results
    TCanvas* can_energy_reco = new TCanvas("ereco", "ereco", 700, 500);
    gPad->SetGrid(1, 1);

    gEnReco->SetMarkerColor(col_exp);
    gEnReco->SetLineColor(col_exp);
    gEnReco->Draw("AP");
    fEnReco->Draw("same");
    fEnReco->SetLineColor(col_exp);

    gEnRecoCorr->SetMarkerColor(col_corr);
    gEnRecoCorr->SetLineColor(col_corr);
    gEnRecoCorr->Draw("P");
    gEnRecoCorr->GetFunction("funPol0")->SetLineColor(col_corr);

    TLegend* leg_1 = new TLegend(0.5, 0.75, 0.9, 0.9);
    leg_1->AddEntry(gEnReco, "reconstructed energy", "PE");
    leg_1->AddEntry(gEnRecoCorr, "reconstructed energy (corrected)", "PE");
    leg_1->Draw();

    double xmin = TMath::Min(TMath::MinElement(npoints, gEnReco->GetY()),
                             TMath::MinElement(npoints, gEnRecoCorr->GetY()));
    double xmax = TMath::Max(TMath::MaxElement(npoints, gEnReco->GetY()),
                             TMath::MaxElement(npoints, gEnRecoCorr->GetY()));

    gEnReco->GetYaxis()->SetRangeUser(xmin-0.05*xmin, xmax+0.05*xmax);

    TCanvas *can_alpha = new TCanvas("ereco_alpha", "ereco_alpha", 700, 500);
    gPad->SetGrid(1, 1);
    
    gAlpha->SetMarkerColor(col_exp);
    gAlpha->SetLineColor(col_exp);
    gAlpha->Draw("AP");
    gAlpha->GetFunction("fun_alpha")->SetLineColor(col_exp);
    
    gAlphaCorr->SetMarkerColor(col_corr);
    gAlphaCorr->SetLineColor(col_corr);
    gAlphaCorr->Draw("P");
    gAlphaCorr->GetFunction("fun_alpha_corr")->SetLineColor(col_corr);
    
    TLegend* leg_2 = new TLegend(0.5, 0.5, 0.90, 0.6);
    leg_2->AddEntry(gAlpha, "#alpha coefficient", "PE");
    leg_2->AddEntry(gAlphaCorr, "#alpha coefficient (corrected)", "PE");
    leg_2->Draw();
    
    xmin = TMath::Min(TMath::MinElement(npoints, gAlpha->GetY()),
                      TMath::MinElement(npoints, gAlphaCorr->GetY()));
    xmax = TMath::Max(TMath::MaxElement(npoints, gAlpha->GetY()),
                      TMath::MaxElement(npoints, gAlphaCorr->GetY()));
    
    gAlpha->GetYaxis()->SetRangeUser(xmin-0.05*xmin, xmax+0.05*xmax);
    
    TLatex text;
    text.SetTextSize(0.04);
    text.SetNDC(true);
    text.SetTextColor(col_exp);
    text.DrawLatex(0.15, 0.5, Form("#alpha = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kAlpha),
                   results[0]->GetUncertainty(SFResultTypeNum::kAlpha)));
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.45, Form("#alpha_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kAlpha),
                   results[1]->GetUncertainty(SFResultTypeNum::kAlpha)));
    
    TCanvas* can_erecospec = new TCanvas("ereco_erecospec", "ereco_erecospec", 700, 500);
    can_erecospec->SetGrid(1, 1);
    
    gEnRecoSpec->SetMarkerColor(col_exp);
    gEnRecoSpec->SetLineColor(col_exp);
    gEnRecoSpec->GetYaxis()->SetRangeUser(0.8, 1.2);
    gEnRecoSpec->Draw("AP");
    
    gEnRecoSpecCorr->SetMarkerColor(col_corr);
    gEnRecoSpecCorr->SetLineColor(col_corr);
    gEnRecoSpecCorr->Draw("P");
    
    xmin = TMath::Min(TMath::MinElement(npoints, gEnRecoSpec->GetY()),
                      TMath::MinElement(npoints, gEnRecoSpecCorr->GetY()));
    xmax = TMath::Max(TMath::MaxElement(npoints, gEnRecoSpec->GetY()),
                      TMath::MaxElement(npoints, gEnRecoSpecCorr->GetY()));
    
    gEnRecoSpec->GetYaxis()->SetRangeUser(xmin-0.1*xmin, xmax+0.1*xmax);
    
    TLegend *leg_3 = new TLegend(0.5, 0.75, 0.9, 0.9);
    leg_3->AddEntry(gEnRecoSpec, "from experimental data", "PE");
    leg_3->AddEntry(gEnRecoSpecCorr, "from corrected data", "PE");
    leg_3->Draw();
    
    TCanvas* can_eres = new TCanvas("ereco_eres", "ereco_eres", 1200, 600);
    can_eres->Divide(2, 1);
    
    can_eres->cd(1);
    gPad->SetGrid(1, 1);
    
    gEnRes->SetTitle(Form("Energy Resolution S%i (from experimental data)", seriesNo));
    gEnRes->SetMarkerColor(col_exp);
    gEnRes->SetLineColor(col_exp);
    gEnRes->Draw("AP");
    
    text.SetTextColor(col_exp);
    text.DrawLatex(0.5, 0.75, Form("ER = (%.3f +/- %.3f) %%",
                   results[0]->GetValue(SFResultTypeNum::kEnergyRes),
                   results[0]->GetUncertainty(SFResultTypeNum::kEnergyRes)));
    
    gEnResAll->SetMarkerColor(col_exp + 6);
    gEnResAll->SetLineColor(col_exp + 6);
    gEnResAll->SetMarkerStyle(23);
    gEnResAll->Draw("P");
    
    text.SetTextColor(col_exp + 6);
    text.SetTextSize(0.03);
    text.SetNDC(false);
    text.DrawLatex(gEnResAll->GetPointX(0) + 0.05, gEnResAll->GetPointY(0) + 0.05,
                   Form("(%.3f +/- %.3f) %%", gEnResAll->GetPointY(0), gEnResAll->GetErrorY(0)));
    text.SetTextSize(0.04);
    text.SetNDC(true);
    
    double min_exp = std::min(TMath::MinElement(npoints, gEnRes->GetY()), eres_all);
    double max_exp = std::max(TMath::MaxElement(npoints, gEnRes->GetY()), eres_all);
     
    gEnRes->GetYaxis()->SetRangeUser(min_exp - 0.05 * min_exp,
                                     max_exp + 0.05 * max_exp);
    
    can_eres->cd(2);
    gPad->SetGrid(1, 1);
    
    gEnResCorr->SetTitle(Form("Energy Resolution S%i (from corrected data)", seriesNo));
    gEnResCorr->SetMarkerColor(col_corr);
    gEnResCorr->SetLineColor(col_corr);
    gEnResCorr->Draw("AP");
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.5, 0.75, Form("ER = (%.3f +/- %.3f) %%",
                   results[1]->GetValue(SFResultTypeNum::kEnergyRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kEnergyRes)));
    
    gEnResCorrAll->SetMarkerColor(col_corr + 5);
    gEnResCorrAll->SetLineColor(col_corr + 5);
    gEnResCorrAll->SetMarkerStyle(23);
    gEnResCorrAll->Draw("P");
    
    text.SetTextColor(col_corr + 5);
    text.SetTextSize(0.03);
    text.SetNDC(false);
    text.DrawLatex(gEnResCorrAll->GetPointX(0) + 0.05, gEnResCorrAll->GetPointY(0) + 0.05,
                   Form("(%.3f +/- %.3f) %%", gEnResCorrAll->GetPointY(0), gEnResCorrAll->GetErrorY(0)));
    text.SetTextSize(0.04);
    text.SetNDC(true);
    
    double min_corr = std::min(TMath::MinElement(npoints, gEnResCorr->GetY()), eres_corr_all);
    double max_corr = std::max(TMath::MaxElement(npoints, gEnResCorr->GetY()), eres_corr_all);
    
    gEnResCorr->GetYaxis()->SetRangeUser(min_corr - 0.05 * min_corr,
                                         max_corr + 0.05 * max_corr);
    
    //-----
    
    //----- drawing energy reconstruction event by event
    TCanvas *can_ereco_byevent_exp = new TCanvas("ereco_byevent_exp", "ereco_byevent_exp",
                                                 2000, 1200);
    can_ereco_byevent_exp->DivideSquare(npoints + 1);
    
    TCanvas *can_ereco_byevent_corr = new TCanvas("ereco_byevent_corr", "ereco_byevent_corr",
                                                  2000, 1200);
    can_ereco_byevent_corr->DivideSquare(npoints + 1);
    
    TCanvas *can_ereco_uncert_corr = new TCanvas("ereco_uncert_corr", "ereco_uncert_corr",
                                                  2000, 1200);
    can_ereco_uncert_corr->DivideSquare(npoints);
    
    text.SetTextColor(kBlack);

    TString fun_name;
    
    double miny = 0.;
    double maxy = maxy = SFTools::FindMaxYaxis(hEnReco[4]);
    
    for (int i = 0; i < npoints; i++)
    {
        can_ereco_byevent_exp->cd(i + 1);
        gPad->SetGrid(1, 1);
        hEnReco[i]->Draw();
        hEnReco[i]->GetXaxis()->SetTitle("reconstructed energy [keV]");
        hEnReco[i]->GetYaxis()->SetTitle("counts");
        hEnReco[i]->SetStats(0);
        fun_name = Form("f_S%i_hEnergyRecoExp_pos%.1f", seriesNo, positions[i]);
        text.DrawLatex(0.5, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gEnRes->GetPointY(i), gEnRes->GetErrorY(i)));
        text.DrawLatex(0.5, 0.75, Form("#mu_{Ereco} = %.2f +/- %.2f", 
                       hEnReco[i]->GetFunction(fun_name)->GetParameter(1),
                       hEnReco[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.5, 0.70, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
                       hEnReco[i]->GetFunction(fun_name)->GetParameter(2),
                       hEnReco[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.5, 0.60, Form("#chi^{2}/NDF = %.3f",
                       hEnReco[i]->GetFunction(fun_name)->GetChisquare() / 
                       hEnReco[i]->GetFunction(fun_name)->GetNDF()));
        hEnReco[i]->GetYaxis()->SetRangeUser(miny, maxy);
        
        can_ereco_byevent_corr->cd(i + 1);
        gPad->SetGrid(1, 1);
        hEnRecoCorr[i]->Draw();
        hEnRecoCorr[i]->GetXaxis()->SetTitle("reconstructed energy [keV]");
        hEnRecoCorr[i]->GetYaxis()->SetTitle("counts");
        hEnRecoCorr[i]->SetStats(0);
        fun_name = Form("f_S%i_hEnergyRecoCorr_pos%.1f", seriesNo, positions[i]);
        text.DrawLatex(0.5, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gEnResCorr->GetPointY(i), gEnResCorr->GetErrorY(i)));
        text.DrawLatex(0.5, 0.75, Form("#mu_{Ereco} = %.2f +/- %.2f", 
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetParameter(1),
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.5, 0.70, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetParameter(2),
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.5, 0.60, Form("#chi^{2}/NDF = %.3f",
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetChisquare() / 
                       hEnRecoCorr[i]->GetFunction(fun_name)->GetNDF()));
        hEnRecoCorr[i]->GetYaxis()->SetRangeUser(miny, maxy);
        
        can_ereco_uncert_corr->cd(i + 1);
        gPad->SetGrid(1, 1);
        hEnRecoUncert[i]->GetYaxis()->SetTitle("counts");
        hEnRecoUncert[i]->GetXaxis()->SetTitle("reconstructed energy uincertainty [keV]");
        hEnRecoUncert[i]->Draw();
    }
    
    can_ereco_byevent_exp->cd(npoints + 1);
    gPad->SetGrid(1, 1);
    gPad->SetFillColor(kGreen - 10);
    hEnRecoAllSpec->Draw();
    hEnRecoAllSpec->SetLineColor(kBlack);
    hEnRecoAllSpec->SetStats(0);
    fun_name = Form("f_S%i_hEnergyRecoAllExp", seriesNo);
    text.DrawLatex(0.55, 0.80, Form("ER = (%.2f +/- %.2f) %%", eres_all, eres_all_err));
    text.DrawLatex(0.55, 0.75, Form("#mu_{Ereco} = %.2f +/- %.2f", 
                   hEnRecoAllSpec->GetFunction(fun_name)->GetParameter(1),
                   hEnRecoAllSpec->GetFunction(fun_name)->GetParError(1)));
    text.DrawLatex(0.55, 0.70, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
                   hEnRecoAllSpec->GetFunction(fun_name)->GetParameter(2),
                   hEnRecoAllSpec->GetFunction(fun_name)->GetParError(2)));
    text.DrawLatex(0.55, 0.60, Form("#chi^{2}/NDF = %.3f",
                   hEnRecoAllSpec->GetFunction(fun_name)->GetChisquare() / 
                   hEnRecoAllSpec->GetFunction(fun_name)->GetNDF()));
    maxy = SFTools::FindMaxYaxis(hEnRecoAllSpec);
    hEnRecoAllSpec->GetYaxis()->SetRangeUser(miny, maxy);
    
    can_ereco_byevent_corr->cd(npoints + 1);
    gPad->SetGrid(1, 1);
    gPad->SetFillColor(kGreen - 10);
    hEnRecoAllCorrSpec->Draw();
    hEnRecoAllCorrSpec->SetLineColor(kBlack);
    hEnRecoAllCorrSpec->SetStats(0);
    fun_name = Form("f_S%i_hEnergyRecoAllCorr", seriesNo);
    text.DrawLatex(0.55, 0.8, Form("ER = (%.2f +/- %.2f) %%", eres_corr_all, eres_corr_all_err));
    text.DrawLatex(0.55, 0.75, Form("#mu_{Ereco} = %.2f +/- %.2f", 
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetParameter(1),
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetParError(1)));
    text.DrawLatex(0.55, 0.70, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetParameter(2),
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetParError(2)));
    text.DrawLatex(0.55, 0.60, Form("#chi^{2}/NDF = %.3f",
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetChisquare() / 
                   hEnRecoAllCorrSpec->GetFunction(fun_name)->GetNDF()));
    maxy = SFTools::FindMaxYaxis(hEnRecoAllCorrSpec);
    hEnRecoAllCorrSpec->GetYaxis()->SetRangeUser(miny, maxy);
    //----- 

    //----- saving
    TString fname       = Form("enreco_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in energyreco.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_energy_reco->Write();
    can_alpha->Write();
    can_erecospec->Write();
    can_eres->Write();
    can_ereco_byevent_corr->Write();
    can_ereco_byevent_exp->Write();
    can_ereco_uncert_corr->Write();
    file->Close();

    //-----writing results to the data base
    TString table = "ENERGY_RECONSTRUCTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ALPHA_EXP, "
                         "ALPHA_EXP_ERR, ALPHA_CORR, ALPHA_CORR_ERR, ERES_EXP, ERES_EXP_ERR, "
                         "ERES_CORR, ERES_CORR_ERR, ERES_ALL, ERES_ALL_ERR, ERES_CORR_ALL, "
                         "ERES_CORR_ALL_ERR) VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, " 
                         "%f, %f, %f, %f)", table.Data(), seriesNo, fname_full.Data(),
                         results[0]->GetValue(SFResultTypeNum::kAlpha),
                         results[0]->GetUncertainty(SFResultTypeNum::kAlpha),
                         results[1]->GetValue(SFResultTypeNum::kAlpha),
                         results[1]->GetUncertainty(SFResultTypeNum::kAlpha),
                         results[0]->GetValue(SFResultTypeNum::kEnergyRes),
                         results[0]->GetUncertainty(SFResultTypeNum::kEnergyRes),
                         results[1]->GetValue(SFResultTypeNum::kEnergyRes),
                         results[1]->GetUncertainty(SFResultTypeNum::kEnergyRes),
                         eres_all, eres_all_err, eres_corr_all, eres_corr_all_err);
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;
    
    do
    {
        std::cout << "----- enres writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
        
        if (stat) break;
        
        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- enres writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    for (auto h : hEnReco)
        delete h;
    
    for (auto h : hEnRecoCorr)
        delete h;
    
    delete data;
    delete reco;
    
    return 0;
}
