// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               model.cc                *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2020              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFAttenuationModel.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

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
        std::cout << "to run type: ./model seriesNo";
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
        std::cerr << "##### Exception in model.cc!" << std::endl;
        return 1;
    }

    data->Print();

    TString desc = data->GetDescription();

    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in model.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cout << "Description: " << desc << std::endl;
        return 1;
    }

    SFAttenuationModel* model = nullptr;
    try
    {
        model = new SFAttenuationModel(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in model.cc!" << std::endl;
        return 1;
    }

    model->Print();
    model->FitModel();

    SFResults* results = model->GetResults();
    results->Print();

    //----- getting model results
    TGraphErrors* gCh0 = (TGraphErrors*)results->GetObject(SFResultTypeObj::kSlVsPosGraph);
    TGraphErrors* gCh1 = (TGraphErrors*)results->GetObject(SFResultTypeObj::kSrVsPosGraph);
    
    gCh0->GetFunction("funCh0")->Delete();
    gCh1->GetFunction("funCh1")->Delete();

    TGraphErrors* gCh0Corr = (TGraphErrors*)results->GetObject(SFResultTypeObj::kPlVsPosGraph);
    TGraphErrors* gCh1Corr = (TGraphErrors*)results->GetObject(SFResultTypeObj::kPrVsPosGraph);

    TF1* funPl = (TF1*)results->GetObject(SFResultTypeObj::kPlFun);
    TF1* funPr = (TF1*)results->GetObject(SFResultTypeObj::kPrFun);
    TF1* funRl = (TF1*)results->GetObject(SFResultTypeObj::kRlFun);
    TF1* funRr = (TF1*)results->GetObject(SFResultTypeObj::kRrFun);
    TF1* funSl = (TF1*)results->GetObject(SFResultTypeObj::kSlFun);
    TF1* funSr = (TF1*)results->GetObject(SFResultTypeObj::kSrFun);
    //-----

    //----- drawing model results
    TCanvas* can_mod_ch = new TCanvas("mod_ch", "mod_ch", 1400, 500);
    can_mod_ch->Divide(2, 1);
    can_mod_ch->cd(1);
    gPad->SetGrid(1, 1);
    
    gCh0->SetMarkerColor(kBlue - 3);
    gCh0->SetLineColor(kBlue - 3);
    gCh0Corr->SetMarkerColor(kBlue - 3);
    gCh0Corr->SetLineColor(kBlue - 3);
    funPl->SetLineColor(kBlue - 3);
    funPl->SetLineStyle(2);
    funRl->SetLineColor(kBlue - 3);
    funRl->SetLineStyle(9);
    funSl->SetLineColor(kBlue - 3);

    gCh1->SetMarkerColor(kTeal - 7);
    gCh1->SetLineColor(kTeal - 7);
    gCh1Corr->SetMarkerColor(kTeal - 7);
    gCh1Corr->SetLineColor(kTeal - 7);
    funPr->SetLineColor(kTeal - 7);
    funPr->SetLineStyle(2);
    funRr->SetLineColor(kTeal - 7);
    funRr->SetLineStyle(9);
    funSr->SetLineColor(kTeal - 7);

    gCh0->Draw("AP");
    gCh1->Draw("P");
    gCh0Corr->Draw("P");
    gCh1Corr->Draw("P");

    TString hname = Form("Channels 0 & 1 Attenuation Curves S%i", seriesNo);
    gCh0->GetYaxis()->SetRangeUser(0, funPl->Eval(0) + (0.3 * funPl->Eval(0)));
    gCh0->SetTitle(hname);

    funPl->Draw("same");
    funPr->Draw("same");
    funRl->Draw("same");
    funRr->Draw("same");
    funSl->Draw("same");
    funSr->Draw("same");

    can_mod_ch->cd(2);

    TLegend* leg_1 = new TLegend(0.1, 0.4, 0.9, 0.9);
    leg_1->AddEntry(gCh0, "experimental data Ch0", "PE");
    leg_1->AddEntry(gCh1, "experimental data Ch1", "PE");
    leg_1->AddEntry(gCh0Corr, "recalculated data Ch0", "PE");
    leg_1->AddEntry(gCh1Corr, "recalculated data Ch1", "PE");
    leg_1->AddEntry(funPl, "direct light Ch0", "L");
    leg_1->AddEntry(funPr, "direct light Ch1", "L");
    leg_1->AddEntry(funRl, "reflected light Ch0", "L");
    leg_1->AddEntry(funRr, "reflected light Ch1", "L");
    leg_1->AddEntry(funSl, "total signal Ch0", "L");
    leg_1->AddEntry(funSr, "total signal Ch1", "L");
    leg_1->Draw();

    TLatex text;
    text.DrawLatex(0.2, 0.35, Form("S_{0} = %.2f +/- %.2f",
                   results->GetValue(SFResultTypeNum::kS0),
                   results->GetUncertainty(SFResultTypeNum::kS0)));
    
    text.DrawLatex(0.2, 0.30, Form("#lambda = %.2f +/- %.2f mm", 
                   results->GetValue(SFResultTypeNum::kLambda),
                   results->GetUncertainty(SFResultTypeNum::kLambda)));
    
    if (results->GetValue(SFResultTypeNum::kEtaR) < 0 ||
        results->GetValue(SFResultTypeNum::kEtaR) > 1)
        text.SetTextColor(kRed);
    
    text.DrawLatex(0.2, 0.25, Form("#eta_{R} = %.3f +/- %.3f", 
                   results->GetValue(SFResultTypeNum::kEtaR),
                   results->GetUncertainty(SFResultTypeNum::kEtaR)));
    text.SetTextColor(kBlack);
    
    if (results->GetValue(SFResultTypeNum::kEtaL) < 0 ||
        results->GetValue(SFResultTypeNum::kEtaL) > 1)
        text.SetTextColor(kRed);
    
    text.DrawLatex(0.2, 0.20, Form("#eta_{L} = %.3f +/- %.3f",
                   results->GetValue(SFResultTypeNum::kEtaL),
                   results->GetUncertainty(SFResultTypeNum::kEtaL)));
    text.SetTextColor(kBlack);
    
    text.DrawLatex(0.2, 0.15, Form("#xi = %.3f +/- %.3f",
                   results->GetValue(SFResultTypeNum::kKsi),
                   results->GetUncertainty(SFResultTypeNum::kKsi)));
    
    text.DrawLatex(0.2, 0.10, Form("L = %.2f mm (fixed)",
                   results->GetValue(SFResultTypeNum::kLength)));
    
    text.DrawLatex(0.2, 0.05, Form("#chi^{2}/NDF = %.3f",
                   results->GetValue(SFResultTypeNum::kChi2NDF)));
    //-----

    //----- saving
    TString fname       = Form("model_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in reconstruction.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_mod_ch->Write();
    file->Close();

    //-----writing results to the data base
    TString table = "ATTENUATION_MODEL";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, S0, S0_ERR, "
                           "LAMBDA, LAMBDA_ERR, ETAR, ETAR_ERR, ETAL, ETAL_ERR, KSI, KSI_ERR, CHI2NDF) "
                           "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)",
                           table.Data(), seriesNo, fname_full.Data(),
                           results->GetValue(SFResultTypeNum::kS0),
                           results->GetUncertainty(SFResultTypeNum::kS0),
                           results->GetValue(SFResultTypeNum::kLambda),
                           results->GetUncertainty(SFResultTypeNum::kLambda),
                           results->GetValue(SFResultTypeNum::kEtaR),
                           results->GetUncertainty(SFResultTypeNum::kEtaR),
                           results->GetValue(SFResultTypeNum::kEtaL),
                           results->GetUncertainty(SFResultTypeNum::kEtaL),
                           results->GetValue(SFResultTypeNum::kKsi),
                           results->GetUncertainty(SFResultTypeNum::kKsi),
                           results->GetValue(SFResultTypeNum::kChi2NDF));
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;
    
    do
    {
        std::cout << "----- model writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);
        
        if (stat) break;
        
        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- model writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    delete data;
    delete model;
    
    return 0;
}
