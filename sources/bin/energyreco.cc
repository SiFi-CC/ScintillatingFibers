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

#include <DistributionContext.h>

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
    double              fiberLen  = data->GetFiberLength();
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
    TGraphErrors* gAlpha  = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kAlphaGraph);
    //TGraphErrors* gEnDiff = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kEnergyDiffGraph);
    TF1*          fEnReco = (TF1*)results[0]->GetObject(SFResultTypeObj::kEnergyRecoFun);
    
    TGraphErrors* gEnRecoCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyRecoGraph);
    TGraphErrors* gAlphaCorr  = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAlphaGraph);
    //TGraphErrors* gEnDiffCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyDiffGraph);
    TF1*          fEnRecoCor  = (TF1*)results[1]->GetObject(SFResultTypeObj::kEnergyRecoFun); 
    
    std::vector <TH1D*> hEnReco     = reco->GetEnergyRecoHistograms("experimental");
    std::vector <TH1D*> hEnRecoCorr = reco->GetEnergyRecoHistograms("corrected");
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

    TLegend* leg_2 = new TLegend(0.5, 0.75, 0.9, 0.9);
    leg_2->AddEntry(gEnReco, "reconstructed energy", "PE");
    leg_2->AddEntry(gEnRecoCorr, "reconstructed energy (corrected)", "PE");
    leg_2->Draw();

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
    
    TLegend* leg_3 = new TLegend(0.5, 0.5, 0.90, 0.6);
    leg_3->AddEntry(gAlpha, "#alpha coefficient", "PE");
    leg_3->AddEntry(gAlphaCorr, "#alpha coefficient (corrected)", "PE");
    leg_3->Draw();
    
    xmin = TMath::Min(TMath::MinElement(npoints, gAlpha->GetY()),
                      TMath::MinElement(npoints, gAlphaCorr->GetY()));
    xmax = TMath::Max(TMath::MaxElement(npoints, gAlpha->GetY()),
                      TMath::MaxElement(npoints, gAlphaCorr->GetY()));
    
    gAlpha->GetYaxis()->SetRangeUser(xmin-0.05*xmin, xmax+0.05*xmax);
    
    TLatex text;
    text.SetNDC(true);
    text.SetTextColor(col_exp);
    text.DrawLatex(0.15, 0.5, Form("#alpha = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kAlpha),
                   results[0]->GetUncertainty(SFResultTypeNum::kAlpha)));
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.45, Form("#alpha_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kAlpha),
                   results[1]->GetUncertainty(SFResultTypeNum::kAlpha)));
    
//     TCanvas* can_energy_diff = new TCanvas("ereco_energy_diff", "ereco_energy_diff", 700, 500);
//     can_energy_diff->SetGrid(1, 1);
//     
//     gEnDiff->SetMarkerColor(col_exp);
//     gEnDiff->SetLineColor(col_exp);
//     gEnDiff->Draw("AP");
//     
//     gEnDiffCorr->SetMarkerColor(col_corr);
//     gEnDiffCorr->SetLineColor(col_corr);
//     gEnDiffCorr->Draw("P");
//     
//     TLegend *leg_4 = new TLegend(0.5, 0.75, 0.9, 0.9);
//     leg_4->AddEntry(gEnDiff, "energy difference", "PE");
//     leg_4->AddEntry(gEnDiffCorr, "energy difference (corrected)", "PE");
//     leg_4->Draw();
//     
//     TLine line;
//     line.SetLineColor(kGray + 2);
//     line.SetLineStyle(9);
//     line.DrawLine(positions[0], 0, positions[npoints - 1], 0);
    //-----
    
    //----- drawing energy reconstruction event by event
    TCanvas *can_ereco_byevent_exp = new TCanvas("ereco_byevent_exp", "ereco_byevent_exp",
                                                 2000, 1200);
    can_ereco_byevent_exp->DivideSquare(npoints);
    
    TCanvas *can_ereco_byevent_corr = new TCanvas("ereco_byevent_corr", "ereco_byevent_corr",
                                                  2000, 1200);
    can_ereco_byevent_corr->DivideSquare(npoints);
    
    text.SetTextColor(kBlack);
    
    for (int i = 0; i < npoints; i++)
    {
        can_ereco_byevent_exp->cd(i + 1);
        gPad->SetGrid(1, 1);
        hEnReco[i]->Draw();
        hEnReco[i]->GetXaxis()->SetTitle("reconstructed energy [keV]");
        hEnReco[i]->GetYaxis()->SetTitle("counts");
        hEnReco[i]->SetStats(0);
        //text.DrawLatex(0.2, 0.8, Form("#mu_{Ereco} = %.2f +/- %.2f", 
        //               hEnReco[i]->GetFunction("fGauss")->GetParameter(1),
        //               hEnReco[i]->GetFunction("fGauss")->GetParError(1)));
        //text.DrawLatex(0.2, 0.75, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
        //               hEnReco[i]->GetFunction("fGauss")->GetParameter(2),
        //               hEnReco[i]->GetFunction("fGauss")->GetParError(2)));
        
        can_ereco_byevent_corr->cd(i + 1);
        gPad->SetGrid(1, 1);
        hEnRecoCorr[i]->Draw();
        hEnRecoCorr[i]->GetXaxis()->SetTitle("reconstructed energy [keV]");
        hEnRecoCorr[i]->GetYaxis()->SetTitle("counts");
        hEnRecoCorr[i]->SetStats(0);
        //text.DrawLatex(0.2, 0.8, Form("#mu_{Ereco} = %.2f +/- %.2f", 
        //               hEnRecoCorr[i]->GetFunction("fGauss")->GetParameter(1),
        //               hEnRecoCorr[i]->GetFunction("fGauss")->GetParError(1)));
        //text.DrawLatex(0.2, 0.75, Form("#sigma_{Ereco} = %.2f +/- %.2f", 
        //               hEnRecoCorr[i]->GetFunction("fGauss")->GetParameter(2),
        //               hEnRecoCorr[i]->GetFunction("fGauss")->GetParError(2)));
    }
    
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
    //can_energy_diff->Write();
    can_ereco_byevent_corr->Write();
    can_ereco_byevent_exp->Write();
    file->Close();

    //-----writing results to the data base
    TString table = "ENERGY_RECONSTRUCTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ALPHA_EXP, ALPHA_EXP_ERR, "
                           "ALPHA_CORR, ALPHA_CORR_ERR) VALUES (%i, '%s', %f, %f, %f, %f)",
                           table.Data(), seriesNo, fname_full.Data(),
                           results[0]->GetValue(SFResultTypeNum::kAlpha),
                           results[0]->GetUncertainty(SFResultTypeNum::kAlpha),
                           results[1]->GetValue(SFResultTypeNum::kAlpha),
                           results[1]->GetUncertainty(SFResultTypeNum::kAlpha));
    
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
    
    delete data;
    delete reco;
    
    return 0;
}
