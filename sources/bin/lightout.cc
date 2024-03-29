// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             lightout.cc               *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *         Katarzyna Rusiecka            *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFData.hh"
#include "SFLightOutput.hh"
#include "SFTools.hh"
#include "common_options.h"

#include <TCanvas.h>
#include <TLatex.h>
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
        std::cout << "to run type: ./lightout seriesNo";
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
        std::cerr << "##### Exception in lightout.cc!" << std::endl;
        return 1;
    }

    TString desc = data->GetDescription();
    if (!desc.Contains("Regular series"))
    {
        std::cerr << "##### Error in lightout.cc! This is not regular series!" << std::endl;
        std::cerr << "Series number: " << seriesNo << std::endl;
        std::cerr << "Description: " << desc << std::endl;
        return 1;
    }

    int                 npoints    = data->GetNpoints();
    //int                 anaGroup   = data->GetAnalysisGroup();
    TString             collimator = data->GetCollimator();
    std::vector<double> positions  = data->GetPositions();
    std::vector<int>    ID         = data->GetMeasurementsIDs();
    data->Print();

    SFLightOutput* lout;
    try
    {
        lout = new SFLightOutput(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in lightout.cc!" << std::endl;
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
    TCanvas* can = lout->GetInputData();

    //----- light output
    lout->CalculateLightOut(0);
    lout->CalculateLightOut(1);
    lout->CalculateLightOut();

    std::vector<SFResults*> LOresults = lout->GetLOResults();
    LOresults[0]->Print(); // channel 0
    LOresults[1]->Print(); // channel 1
    LOresults[2]->Print(); // summed

    TGraphErrors* gLightOutCh0 = (TGraphErrors*)LOresults[0]->GetObject(SFResultTypeObj::kLightGraph);
    TGraphErrors* gLightOutCh1 = (TGraphErrors*)LOresults[1]->GetObject(SFResultTypeObj::kLightGraph);
    TGraphErrors* gLightOut    = (TGraphErrors*)LOresults[2]->GetObject(SFResultTypeObj::kLightGraph);

    std::vector<TH1D*> specCh0 = lout->GetSpectra(0);
    std::vector<TH1D*> specCh1 = lout->GetSpectra(1);

    //----- light collection
    lout->CalculateLightCol(0);
    lout->CalculateLightCol(1);
    lout->CalculateLightCol();

    std::vector<SFResults*> LCresults = lout->GetLCResults();
    LCresults[0]->Print(); // channel 0
    LCresults[1]->Print(); // channel 1
    LCresults[2]->Print(); // summed

    TGraphErrors* gLightColCh0 = (TGraphErrors*)LCresults[0]->GetObject(SFResultTypeObj::kLightGraph);
    TGraphErrors* gLightColCh1 = (TGraphErrors*)LCresults[1]->GetObject(SFResultTypeObj::kLightGraph);
    TGraphErrors* gLightCol    = (TGraphErrors*)LCresults[2]->GetObject(SFResultTypeObj::kLightGraph);

    //----- drawing
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.04);

    //----- light output channels 0 and 1
    TCanvas* can_lout_ch = new TCanvas("lo_lout_ch", "lo_lout_ch", 1200, 600);
    can_lout_ch->Divide(2, 1);

    can_lout_ch->cd(1);
    gPad->SetGrid(1, 1);
    gLightOutCh0->GetYaxis()->SetMaxDigits(2);
    gLightOutCh0->GetXaxis()->SetLabelSize(0.035);
    gLightOutCh0->GetYaxis()->SetLabelSize(0.035);
    gLightOutCh0->GetYaxis()->SetTitleOffset(1.4);
    gLightOutCh0->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", 
                  LOresults[0]->GetValue(SFResultTypeNum::kLight),
                  LOresults[0]->GetUncertainty(SFResultTypeNum::kLight)));

    can_lout_ch->cd(2);
    gPad->SetGrid(1, 1);
    gLightOutCh1->GetYaxis()->SetMaxDigits(2);
    gLightOutCh1->GetXaxis()->SetLabelSize(0.035);
    gLightOutCh1->GetYaxis()->SetLabelSize(0.035);
    gLightOutCh1->GetYaxis()->SetTitleOffset(1.4);
    gLightOutCh1->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", 
                   LOresults[1]->GetValue(SFResultTypeNum::kLight),
                   LOresults[1]->GetUncertainty(SFResultTypeNum::kLight)));

    //----- light output summed
    TCanvas* can_lout = new TCanvas("lo_lout", "lo_lout", 700, 500);
    can_lout->cd();
    gPad->SetGrid(1, 1);
    gLightOut->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LO = (%.2f +/- %.2f) PE/MeV", 
                   LOresults[2]->GetValue(SFResultTypeNum::kLight),
                   LOresults[2]->GetUncertainty(SFResultTypeNum::kLight)));

    //----- light collection channels 0 and 1
    TCanvas* can_lcol_ch = new TCanvas("lo_lcol_ch", "lo_lcol_ch", 1200, 600);
    can_lcol_ch->Divide(2, 1);

    can_lcol_ch->cd(1);
    gPad->SetGrid(1, 1);
    gLightColCh0->GetYaxis()->SetMaxDigits(2);
    gLightColCh0->GetXaxis()->SetLabelSize(0.035);
    gLightColCh0->GetYaxis()->SetLabelSize(0.035);
    gLightColCh0->GetYaxis()->SetTitleOffset(1.4);
    gLightColCh0->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LC = (%.2f +/- %.2f) PE/MeV", 
                   LCresults[0]->GetValue(SFResultTypeNum::kLight),
                   LCresults[0]->GetUncertainty(SFResultTypeNum::kLight)));

    can_lcol_ch->cd(2);
    gPad->SetGrid(1, 1);
    gLightColCh1->GetYaxis()->SetMaxDigits(2);
    gLightColCh1->GetXaxis()->SetLabelSize(0.035);
    gLightColCh1->GetYaxis()->SetLabelSize(0.035);
    gLightColCh1->GetYaxis()->SetTitleOffset(1.4);
    gLightColCh1->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LC = (%.2f +/- %.2f) PE/MeV", 
                   LCresults[1]->GetValue(SFResultTypeNum::kLight),
                   LCresults[1]->GetUncertainty(SFResultTypeNum::kLight)));

    //----- light collection summed
    TCanvas* can_lcol = new TCanvas("lo_lcol", "lo_lcol", 700, 500);
    can_lcol->cd();
    gPad->SetGrid(1, 1);
    gLightCol->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("LC = (%.2f +/- %.2f) PE/MeV", 
                   LCresults[2]->GetValue(SFResultTypeNum::kLight),
                   LCresults[2]->GetUncertainty(SFResultTypeNum::kLight)));

    //----- spectra
    TCanvas* can_spec_ch0 = new TCanvas("lo_spec_ch0", "lo_spec_ch0", 2000, 1200);
    can_spec_ch0->DivideSquare(npoints);

    TCanvas* can_spec_ch1 = new TCanvas("lo_spec_ch1", "lo_spec_ch1", 2000, 1200);
    can_spec_ch1->DivideSquare(npoints);
    
    double min_yaxis = 0.;
    double max_yaxis_0 =  SFTools::FindMaxYaxis(specCh0[npoints-1]);
    double max_yaxis_1 =  SFTools::FindMaxYaxis(specCh1[0]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(specCh0[0]);
    
    TString fun_name;
    
    for (int i = 0; i < npoints; i++)
    {
        can_spec_ch0->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh0[i]->SetStats(false);
        specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specCh0[i]->GetYaxis()->SetTitle("counts");
        specCh0[i]->SetTitle(Form("Cherge spectrum Ch0, position %.2f mm", positions[i]));
        //ctx.configureFromJson("hSpec");
        //ctx.print();
        //specCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
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

        can_spec_ch1->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh1[i]->SetStats(false);
        specCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specCh1[i]->GetYaxis()->SetTitle("counts");
        specCh1[i]->SetTitle(Form("Cherge spectrum Ch1, position %.2f mm", positions[i]));
        //specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
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
    }

    //----- saving
    TString fname       = Form("lightout_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in lightout.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_lout_ch->Write();
    can_lout->Write();
    can_lcol_ch->Write();
    can_lcol->Write();
    can_spec_ch0->Write();
    can_spec_ch1->Write();
    can->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "LIGHT_OUTPUT";
    TString query = Form(
        "INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, "
        "LOUT, LOUT_ERR, LOUT_CH0, LOUT_CH0_ERR, LOUT_CH1, LOUT_CH1_ERR, "
        "LCOL, LCOL_ERR, LCOL_CH0, LCOL_CH0_ERR, LCOL_CH1, LCOL_CH1_ERR) "
        "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)",
        table.Data(), seriesNo, fname_full.Data(), 
        LOresults[2]->GetValue(SFResultTypeNum::kLight),
        LOresults[2]->GetUncertainty(SFResultTypeNum::kLight),
        LOresults[0]->GetValue(SFResultTypeNum::kLight),
        LOresults[0]->GetUncertainty(SFResultTypeNum::kLight),
        LOresults[1]->GetValue(SFResultTypeNum::kLight),
        LOresults[1]->GetUncertainty(SFResultTypeNum::kLight),
        LCresults[2]->GetValue(SFResultTypeNum::kLight), 
        LCresults[2]->GetUncertainty(SFResultTypeNum::kLight),
        LCresults[0]->GetValue(SFResultTypeNum::kLight),
        LCresults[0]->GetUncertainty(SFResultTypeNum::kLight),
        LCresults[1]->GetValue(SFResultTypeNum::kLight),
        LCresults[1]->GetUncertainty(SFResultTypeNum::kLight));

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- lightout writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- lightout writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    delete lout;

    return 0;
}
