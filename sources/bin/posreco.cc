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

    TString desc     = data->GetDescription();
    int     npoints  = data->GetNpoints();

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
    TGraphErrors* gMLRExp = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kMLRvsPosGraph);
    TGraphErrors* gPRecoExp = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPResExp = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualExp = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kResidualGraph);
    
    std::vector<TH1D*> hPosDistExp = reco->GetPositionDistributions("experimental");
    
    TGraphErrors* gMLRCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kMLRvsPosGraph);
    TGraphErrors* gPRecoCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors* gPResCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors* gResidualCorr = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kResidualGraph);
    TGraphErrors* gACoeff = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAGraph);
    
    std::vector<TH1D*> hPosDistCorr = reco->GetPositionDistributions("corrected");
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
    
    gMLRExp->GetYaxis()->SetRangeUser(-0.5, 0.5);
    
    TLegend* leg_1 = new TLegend(0.5, 0.15, 0.9, 0.25);
    leg_1->AddEntry(gMLRExp, "M_{LR} experimental", "PE");
    leg_1->AddEntry(gMLRCorr, "M_{LR} corrected", "PE");
    leg_1->Draw();
    //-----
    
    //----- drawing A graph
    TCanvas *can_a = new TCanvas("preco_a", "preco_a", 700, 500);
    gPad->SetGrid(1, 1);
    
    gACoeff->SetMarkerColor(col_corr);
    gACoeff->SetLineColor(col_corr);
    gACoeff->Draw("AP");
    
    gACoeff->GetFunction("fpol0")->SetLineColor(col_corr);
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.75, Form("A = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kACoeff),
                   results[1]->GetUncertainty(SFResultTypeNum::kACoeff)));
    text.DrawLatex(0.15, 0.70, Form("B = %.3f",
                   results[1]->GetValue(SFResultTypeNum::kBCoeff)));
    //-----
    
    //----- drawing position resolution
    TCanvas *can_pos_res = new TCanvas("preco_pos_res", "preco_pos_res", 1200, 600);
    can_pos_res->Divide(2, 1);
    
    can_pos_res->cd(1);
    gPad->SetGrid(1, 1);
    gPResExp->SetMarkerColor(col_exp);
    gPResExp->SetLineColor(col_exp);
    gPResExp->Draw("AP");
    
    text.SetTextColor(col_exp);
    text.DrawLatex(0.15, 0.75, Form("PR_{exp} = %.3f +/- %.3f",
                   results[0]->GetValue(SFResultTypeNum::kPositionRes),
                   results[0]->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    can_pos_res->cd(2);
    gPad->SetGrid(1, 1);
    gPResCorr->SetMarkerColor(col_corr);
    gPResCorr->SetLineColor(col_corr);
    gPResCorr->Draw("AP");
    
    text.SetTextColor(col_corr);
    text.DrawLatex(0.15, 0.75, Form("PR_{corr} = %.3f +/- %.3f",
                   results[1]->GetValue(SFResultTypeNum::kPositionRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kPositionRes)));
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
    gPRecoExp->GetFunction("funpol1")->SetLineColor(col_exp);
    gPRecoExp->Draw("AP");
    
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.3f +/- %.3f",  
                   gPRecoExp->GetFunction("funpol1")->GetParameter(0),
                   gPRecoExp->GetFunction("funpol1")->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.3f +/- %.3f",  
                   gPRecoExp->GetFunction("funpol1")->GetParameter(1),
                   gPRecoExp->GetFunction("funpol1")->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.3f",  
                   gPRecoExp->GetFunction("funpol1")->GetChisquare() /
                   gPRecoExp->GetFunction("funpol1")->GetNDF()));
    
    can_pos_reco->cd(1);
    pad_res_1->Draw();
    pad_res_1->cd();
    gPad->SetGrid(1, 1);
    gResidualExp->SetMarkerColor(col_exp);
    gResidualExp->SetLineColor(col_exp);
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

    //----- drawing position distributions
    TCanvas *can_pos_dist = new TCanvas("preco_pos_dist", "preco_pos_dist", 2000, 1200);
    can_pos_dist->DivideSquare(npoints);
    
    text.SetTextColor(kBlack);
    text.SetTextSize(0.035);
    
    for (int i = 0; i < npoints; i++)
    {
        can_pos_dist->cd(i + 1);
        gPad->SetGrid(1, 1);
        
        hPosDistCorr[i]->SetStats(0);
        hPosDistCorr[i]->GetYaxis()->SetTitle("counts");
        hPosDistCorr[i]->GetXaxis()->SetTitle("reconstructed position [mm]");
        hPosDistCorr[i]->Draw();
        
        text.DrawLatex(0.2, 0.80, Form("#mu = (%.3f +/- %.3f) mm",
                       gPRecoCorr->GetPointY(i),
                       gPRecoCorr->GetErrorY(i)));
        text.DrawLatex(0.2, 0.75, Form("FWHM = (%.3f +/- %.3f) mm",
                       gPResCorr->GetPointY(i),
                       gPResCorr->GetErrorY(i)));
    }
    
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
    can_pos_dist->Write();
    file->Close();

    //-----writing results to the data base
    TString table = "POSITION_RECONSTRUCTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, A_COEFF, "
                         "A_COEFF_ERR, B_COEFF, MLR_SLOPE, MLR_SLOPE_ERR, MLR_OFFSET, " 
                         "MLR_OFFSET_ERR, MLR_SLOPE_EXP, MLR_SLOPE_EXP_ERR, MLR_OFFSET_EXP, "
                         "MLR_OFFSET_EXP_ERR, POSITION_RES, POSITION_RES_ERR) "
                         "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, "
                         "%f, %f, %f)",
                         table.Data(), seriesNo, fname_full.Data(),
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
                         results[1]->GetUncertainty(SFResultTypeNum::kPositionRes));
    
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
