// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *             energyres.cc              *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFEnergyRes.hh"
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
        std::cout << "to run type: ./energyres seriesNo ";
        std::cout << "-out path/to/output -db database" << std::endl;
        return 1;
    }

    //----- accessing results of energy resolution analysis
    SFData* data;
    try
    {
        data = new SFData(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cout << "##### Exception in energyres.cc!" << std::endl;
        return 1;
    }

    data->Print();
    int                 npoints    = data->GetNpoints();
    //int                 anaGroup   = data->GetAnalysisGroup();
    TString             collimator = data->GetCollimator();
    std::vector<double> positions  = data->GetPositions();
    std::vector<int>    ID         = data->GetMeasurementsIDs();

    SFEnergyRes* enres;
    try
    {
        enres = new SFEnergyRes(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in energyres.cc!" << std::endl;
        return 1;
    }

    enres->CalculateEnergyRes(0);
    enres->CalculateEnergyRes(1);
    enres->CalculateEnergyRes();

    std::vector<SFResults*> results = enres->GetResults();
    results[0]->Print(); // channel 0
    results[1]->Print(); // channel 1
    results[2]->Print(); // averaged

    TGraphErrors* gEnResAve = (TGraphErrors*)results[2]->GetObject(SFResultTypeObj::kEnergyResGraph);
    TGraphErrors* gEnResCh0 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kEnergyResGraph);
    TGraphErrors* gEnResCh1 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kEnergyResGraph);

    std::vector<TH1D*> specAve = enres->GetSpectra();
    std::vector<TH1D*> specCh0 = enres->GetSpectra(0);
    std::vector<TH1D*> specCh1 = enres->GetSpectra(1);

    //----- accessing json file
/*    
   DistributionContext ctx;
    ctx.findJsonFile("./", Form(".configAG%i.json", anaGroup));

    ctx.dim   = DIM1;
    ctx.x.min = 0;
    ctx.x.max = 100;
    ctx.y.min = 0;
    ctx.y.max = 1000;
*/
    //----- drawing energy resolution graphs
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.04);

    //----- averaged channels
    TCanvas* can_ave = new TCanvas("er_ave", "er_ave", 700, 500);
    can_ave->cd();
    gPad->SetGrid(1, 1);
    gEnResAve->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                  results[2]->GetValue(SFResultTypeNum::kEnergyRes),
                  results[2]->GetUncertainty(SFResultTypeNum::kEnergyRes)));

    //----- channel 0
    TCanvas* can_ch0 = new TCanvas("er_ch0", "er_ch0", 700, 500);
    can_ch0->cd();
    gPad->SetGrid(1, 1);
    gEnResCh0->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                   results[0]->GetValue(SFResultTypeNum::kEnergyRes),
                   results[0]->GetUncertainty(SFResultTypeNum::kEnergyRes)));

    //----- channel 1
    TCanvas* can_ch1 = new TCanvas("er_ch1", "er_ch1", 700, 500);
    can_ch1->cd();
    gPad->SetGrid(1, 1);
    gEnResCh1->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                   results[1]->GetValue(SFResultTypeNum::kEnergyRes),
                   results[1]->GetUncertainty(SFResultTypeNum::kEnergyRes)));

    //----- drawing spectra
    TCanvas* can_spec_ave = new TCanvas("er_spec_ave", "er_spec_ave", 2000, 1200);
    can_spec_ave->DivideSquare(npoints);

    TCanvas* can_spec_ch0 = new TCanvas("er_spec_ch0", "er_spec_ch0", 2000, 1200);
    can_spec_ch0->DivideSquare(npoints);

    TCanvas* can_spec_ch1 = new TCanvas("er_spec_ch1", "er_spec_ch1", 2000, 1200);
    can_spec_ch1->DivideSquare(npoints);
    
    double min_yaxis = 0.;
    double max_yaxis_0  =  SFTools::FindMaxYaxis(specCh0[npoints-1]);
    double max_yaxis_1  =  SFTools::FindMaxYaxis(specCh1[0]);
    double max_yaxis_av =  SFTools::FindMaxYaxis(specAve[4]); 
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(specCh0[0]);

    TString fun_name;
    
    for (int i = 0; i < npoints; i++)
    {
        can_spec_ave->cd(i + 1);
        gPad->SetGrid(1, 1);
        specAve[i]->SetStats(false);
        specAve[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specAve[i]->GetYaxis()->SetTitle("counts");
        specAve[i]->SetTitle(Form("Average PE spectrum, position %.2f mm", positions[i]));
        //ctx.configureFromJson("hSpecAv");
        //ctx.print();
        //specAve[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specAve[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        specAve[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specAve[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_av);
        specAve[i]->Draw();
        fun_name = Form("f_S%i_pos%.1f_ID%i_PEAverage", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gEnResAve->GetPointY(i), gEnResAve->GetErrorY(i)));
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       specAve[i]->GetFunction(fun_name)->GetParameter(1),
                       specAve[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       specAve[i]->GetFunction(fun_name)->GetParameter(2),
                       specAve[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       specAve[i]->GetFunction(fun_name)->GetChisquare() / 
                       specAve[i]->GetFunction(fun_name)->GetNDF()));

        can_spec_ch0->cd(i + 1);
        gPad->SetGrid(1, 1);
        specCh0[i]->SetStats(false);
        specCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        specCh0[i]->GetYaxis()->SetTitle("counts");
        specCh0[i]->SetTitle(Form("Ch0 PE spectrum, position %.2f mm", positions[i]));
        //ctx.configureFromJson("hSpec");
        //ctx.print();
        //specCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        specCh0[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh0[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        specCh0[i]->Draw();
        fun_name = Form("f_S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gEnResCh0->GetPointY(i), gEnResCh0->GetErrorY(i)));
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
        specCh1[i]->SetTitle(Form("Ch1 PE spectrum, position %.2f mm", positions[i]));
        //specCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //specCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        specCh1[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        specCh1[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        specCh1[i]->Draw();
        fun_name = Form("f_S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]);
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gEnResCh1->GetPointY(i), gEnResCh1->GetErrorY(i)));
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

    //----- saving ROOT file
    TString fname       = Form("enres_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in energyres.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_ave->Write();
    can_ch0->Write();
    can_ch1->Write();
    can_spec_ave->Write();
    can_spec_ch0->Write();
    can_spec_ch1->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "ENERGY_RESOLUTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ENRES_AV, "
                         "ENRES_AV_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR) "
                         "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f)",
                         table.Data(), seriesNo, fname_full.Data(),
                         results[2]->GetValue(SFResultTypeNum::kEnergyRes),
                         results[2]->GetUncertainty(SFResultTypeNum::kEnergyRes),
                         results[0]->GetValue(SFResultTypeNum::kEnergyRes),
                         results[0]->GetUncertainty(SFResultTypeNum::kEnergyRes),
                         results[1]->GetValue(SFResultTypeNum::kEnergyRes),
                         results[1]->GetUncertainty(SFResultTypeNum::kEnergyRes));

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- energyres writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- energyres writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    delete enres;

    return 0;
}
