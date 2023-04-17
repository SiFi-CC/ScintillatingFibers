// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            attenuation.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2018              *
// *                                       *
// *****************************************

#include "SFAttenuation.hh"
#include "SFData.hh"
#include "common_options.h"

//#include <DistributionContext.h>

#include <TCanvas.h>
#include <TLatex.h>

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
        std::cout << "to run type: ./attenuation seriesNo";
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

    int                 npoints    = data->GetNpoints();
    //int               anaGroup   = data->GetAnalysisGroup();
    TString             collimator = data->GetCollimator();
    TString             sipm       = data->GetSiPM();
    std::vector<double> positions  = data->GetPositions();
    std::vector<int>    ID         = data->GetMeasurementsIDs();

    SFAttenuation* att;
    try
    {
        att = new SFAttenuation(seriesNo);
    }
    catch (const char* message)
    {
        std::cerr << message << std::endl;
        std::cerr << "##### Exception in attenuation.cc!" << std::endl;
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
    //----- analysis
    att->AttCombinedCh();
    att->AttSeparateCh(0);
    att->AttSeparateCh(1);
    att->FitSimultaneously();

    //----- accessing results
    std::vector<SFResults*> results = att->GetResults();
    results[0]->Print(); // channel 0
    results[1]->Print(); // channel 1
    results[2]->Print(); // combined channels/pol1
    results[3]->Print(); // combined channels/pol3
    results[4]->Print(); // simple model sim. fitting

    TGraphErrors* attGraphPol1 = (TGraphErrors*)results[2]->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors* attGraphPol3 = (TGraphErrors*)attGraphPol1->Clone("attGraphPol3");
    TGraphErrors* sigGraph     = (TGraphErrors*)results[2]->GetObject(SFResultTypeObj::kMLRSigmaGraph);

    TGraphErrors* attGraphCh0 = (TGraphErrors*)results[0]->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors* attGraphCh1 = (TGraphErrors*)results[1]->GetObject(SFResultTypeObj::kAttGraph);
    
    TF1 *funSl = (TF1*)results[4]->GetObject(SFResultTypeObj::kSlFun);
    TF1 *funSr = (TF1*)results[4]->GetObject(SFResultTypeObj::kSrFun);

    std::vector<TH1D*> attRatios  = att->GetRatios();
    std::vector<TH1D*> spectraCh1 = att->GetSpectra(1);
    std::vector<TH1D*> spectraCh0 = att->GetSpectra(0);

    //-----drawing averaged channels
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.03);

    TCanvas* can_averaged_ch = new TCanvas("att_averaged_ch", "att_averaged_ch", 1200, 600);
    can_averaged_ch->Divide(2, 1);

    can_averaged_ch->cd(1);
    gPad->SetGrid(1, 1);
    attGraphPol1->SetTitle(Form("Attenuation Curve S%i", seriesNo));
    attGraphPol1->GetFunction("fpol3")->Delete();
    attGraphPol1->GetYaxis()->SetTitleOffset(1.4);
    attGraphPol1->GetXaxis()->SetLabelSize(0.035);
    attGraphPol1->GetYaxis()->SetLabelSize(0.035);
    attGraphPol1->Draw("AP");
    text.DrawLatex(0.2, 0.85, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.80, Form("p_{0} = %.3e +/- %.3e", 
                   attGraphPol1->GetFunction("fpol1")->GetParameter(0),
                   attGraphPol1->GetFunction("fpol1")->GetParError(0)));
    text.DrawLatex(0.2, 0.75, Form("p_{1} = %.3e +/- %.3e", 
                   attGraphPol1->GetFunction("fpol1")->GetParameter(1),
                   attGraphPol1->GetFunction("fpol1")->GetParError(1)));
    text.DrawLatex(0.2, 0.70, Form("L_{att} = (%.2f +/- %.2f) mm",
                   results[2]->GetValue(SFResultTypeNum::kLambda),
                   results[2]->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   results[2]->GetValue(SFResultTypeNum::kChi2NDF)));

    can_averaged_ch->cd(2);
    gPad->SetGrid(1, 1);
    attGraphPol3->SetTitle(Form("Attenuation Curve S%i", seriesNo));
    attGraphPol3->GetFunction("fpol1")->Delete();
    attGraphPol3->GetYaxis()->SetTitleOffset(1.2);
    attGraphPol3->Draw("AP");
    TF1* funPol3 = (TF1*)attGraphPol3->GetFunction("fpol3");
    text.DrawLatex(0.2, 0.85, "y = p_{0} + p_{1}[x + p_{2}x^{3}]");
    text.DrawLatex(0.2, 0.8, Form("p_{0} = %.4e +/- %.4e", 
                   funPol3->GetParameter(0), funPol3->GetParError(0)));
    text.DrawLatex(0.2, 0.75, Form("p_{1} = %.4e +/- %.4e", 
                   funPol3->GetParameter(1), funPol3->GetParError(1)));
    text.DrawLatex(0.2, 0.70, Form("p_{2} = %.4e +/- %.4e", 
                   funPol3->GetParameter(2), funPol3->GetParError(2)));
    text.DrawLatex(0.2, 0.60, Form("L_{att} = (%.2f +/- %.2f) mm",
                   results[3]->GetValue(SFResultTypeNum::kLambda),
                   results[3]->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.2, 0.55, Form("#chi^{2}/NDF = %.2f", 
                   results[3]->GetValue(SFResultTypeNum::kChi2NDF)));

    TCanvas* can_sigma = new TCanvas("att_sigma", "att_sigma", 700, 500);
    gPad->SetGrid(1, 1);
    sigGraph->SetTitle(Form("Sigma of M_{LR} Distribution S%i", seriesNo));
    sigGraph->GetYaxis()->SetTitleOffset(1.3);
    sigGraph->Draw("AP");

    TCanvas* can_ratios = new TCanvas("att_ratios", "att_ratios", 2000, 1200);
    can_ratios->Divide(3, 3);
    TString htitle;
    text.SetTextSize(0.02);

    TF1* fun;
    TF1* fun_clone = new TF1("fun_clone", "gaus", -1.5, 1.5);
    TF1* fthin     = new TF1("fthin", "gaus", -1.5, 1.5);
    TF1* fthick    = new TF1("fthick", "gaus", -1.5, 1.5);
    text.SetTextSize(0.045);

    for (int i = 0; i < npoints; i++)
    {
        can_ratios->cd(i + 1);
        gPad->SetGrid(1, 1);
        attRatios[i]->SetTitle(Form("M_{LR}, S%i %.2f mm", seriesNo, positions[i]));
        attRatios[i]->GetXaxis()->SetRangeUser(-1.5, 1.5);
        attRatios[i]->GetXaxis()->SetTitle("M_{LR}");
        attRatios[i]->Draw();
        if ((collimator == "Lead") || (collimator == "Electronic" && sipm == "SensL"))
        {
            fun = attRatios[i]->GetFunction("fDGauss");
            fthin->SetParameters(fun->GetParameter(0), fun->GetParameter(1), fun->GetParameter(2));
            fthick->SetParameters(fun->GetParameter(3), fun->GetParameter(4), fun->GetParameter(5));
            fthin->SetLineColor(kMagenta);
            fthick->SetLineColor(kMagenta - 10);
            fthin->DrawClone("same");
            fthick->DrawClone("same");
            text.DrawLatex(0.2, 0.80,
                           Form("sig/bg = %.2f", fun->GetParameter(0) / fun->GetParameter(3)));
            text.DrawLatex(0.2, 0.75, Form("sig = %.4f", fun->GetParameter(1)));
            text.DrawLatex(0.2, 0.70, Form("bg = %.4f", fun->GetParameter(4)));
        }
        else if (collimator == "Electronic" && sipm == "Hamamatsu")
        {
            fun = attRatios[i]->GetFunction("fGauss");
            fun_clone->SetParameters(fun->GetParameter(0), fun->GetParameter(1),
                                     fun->GetParameter(2));
            fun_clone->SetLineColor(kMagenta - 10);
            fun_clone->DrawClone("same");
        }
    }

    //----- drawing separate channels
    text.SetTextSize(0.035);
    TCanvas* can_separate_ch = new TCanvas("att_separate_ch", "att_separate_ch", 1200, 600);
    can_separate_ch->Divide(2,1);
    
    can_separate_ch->cd(1);
    gPad->SetGrid(1, 1);

    TString hname = Form("Channels 0 & 1 Attenuation Curves S%i (Separate Fitting)", seriesNo);
    TH1D *h = new TH1D("h", hname, 100, 0, 100);
    h->GetXaxis()->SetTitle("source position [mm]");
    h->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    h->SetStats(false);
    h->Draw();
    
    attGraphCh0->SetTitle(Form("Channel 0 Attenuation Curve S%i", seriesNo));
    attGraphCh0->GetYaxis()->SetTitleSize(0.03);
    attGraphCh0->SetMarkerColor(kPink - 8);
    attGraphCh0->SetLineColor(kPink - 8);
    attGraphCh0->Draw("P");
    attGraphCh0->GetFunction("funCh0")->SetLineColor(kPink - 8);
    text.SetTextColor(kPink - 8);
    text.DrawLatex(0.3, 0.8, Form("L_{att Ch0} = (%.2f +/- %.2f) mm",
                   results[0]->GetValue(SFResultTypeNum::kLambda),
                   results[0]->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.65, Form("#chi^{2}/NDF_{0} = %.2f", 
                   results[0]->GetValue(SFResultTypeNum::kChi2NDF)));

    attGraphCh1->SetTitle(Form("Channel 1 Attenuation Curve S%i", seriesNo));
    attGraphCh1->GetYaxis()->SetTitleSize(0.03);
    attGraphCh1->SetMarkerColor(kAzure - 6);
    attGraphCh1->SetLineColor(kAzure - 6);
    attGraphCh1->Draw("P");
    attGraphCh1->GetFunction("funCh1")->SetLineColor(kAzure - 6);
    text.SetTextColor(kAzure - 6);
    text.DrawLatex(0.3, 0.75, Form("L_{att Ch1} = (%.2f +/- %.2f) mm",
                   results[1]->GetValue(SFResultTypeNum::kLambda),
                   results[1]->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.60, Form("#chi^{2}/NDF_{1} = %.2f", 
                   results[1]->GetValue(SFResultTypeNum::kChi2NDF)));

    double* yCh0    = attGraphCh0->GetY();
    double* yCh1    = attGraphCh1->GetY();
    double  yminCh0 = TMath::MinElement(npoints, yCh0);
    double  yminCh1 = TMath::MinElement(npoints, yCh1);
    double  ymin    = TMath::Min(yminCh0, yminCh1);
    double  ymaxCh0 = TMath::MaxElement(npoints, yCh0);
    double  ymaxCh1 = TMath::MaxElement(npoints, yCh1);
    double  ymax    = TMath::Max(ymaxCh0, ymaxCh1);

    h->GetYaxis()->SetRangeUser(ymin - 0.2 * ymin, ymax + 0.1 * ymax);
    
    can_separate_ch->cd(2);
    gPad->SetGrid(1, 1);

    hname = Form("Channels 0 & 1 Attenuation Curves S%i (Simultaneous Fitting)", seriesNo);
    TH1D *hh = new TH1D("hh", hname, 100, 0, 100);
    hh->GetXaxis()->SetTitle("source position [mm]");
    hh->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    hh->SetStats(false);
    hh->Draw();
    
    TGraphErrors *attGraphCh0Clone = (TGraphErrors*)attGraphCh0->Clone();
    attGraphCh0Clone->GetFunction("funCh0")->Delete();
    attGraphCh0Clone->Draw("P");
    funSl->SetLineColor(kPink - 8);
    funSl->Draw("same");
    
    TGraphErrors *attGraphCh1Clone = (TGraphErrors*)attGraphCh1->Clone();
    attGraphCh1Clone->GetFunction("funCh1")->Delete();
    attGraphCh1Clone->Draw("P");
    funSr->SetLineColor(kAzure - 6);
    funSr->Draw("same");
    
    text.SetTextColor(kBlack);
    text.DrawLatex(0.3, 0.75, Form("L_{att} = (%.2f +/- %.2f) mm",
                   results[4]->GetValue(SFResultTypeNum::kLambda),
                   results[4]->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.60, Form("#chi^{2}/NDF_{1} = %.2f", 
                   results[4]->GetValue(SFResultTypeNum::kChi2NDF)));
    
    hh->GetYaxis()->SetRangeUser(ymin - 0.2 * ymin, ymax + 0.1 * ymax);
    
    TCanvas* can_spectra_ch0 = new TCanvas("att_spectra_ch0", "att_spectra_ch0", 2000, 1200);
    can_spectra_ch0->DivideSquare(npoints);

    TCanvas* can_spectra_ch1 = new TCanvas("att_spectra_ch1", "att_spectra_ch1", 2000, 1200);
    can_spectra_ch1->DivideSquare(npoints);

    double min_yaxis = 0.;
    double max_yaxis_0 =  SFTools::FindMaxYaxis(spectraCh0[npoints-1]);
    double max_yaxis_1 =  SFTools::FindMaxYaxis(spectraCh1[0]);
    
    double min_xaxis = 10.;
    double max_xaxis = SFTools::FindMaxXaxis(spectraCh0[0]);
    
    text.SetTextColor(kBlack);
    text.SetTextSize(0.040);
    TString fun_name;
    
    for (int i = 0; i < npoints; i++)
    {
        can_spectra_ch0->cd(i + 1);
        gPad->SetGrid(1, 1);
        spectraCh0[i]->SetStats(false);
        spectraCh0[i]->GetXaxis()->SetRangeUser(0, 1200);
        spectraCh0[i]->SetTitle(Form("PE spectrum S%i Ch0, source position %.2f mm",
                                     seriesNo, positions[i]));
        spectraCh0[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        spectraCh0[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        spectraCh0[i]->GetXaxis()->SetTitle("charge [P.E.]");
        spectraCh0[i]->GetYaxis()->SetTitle("counts");
        //ctx.configureFromJson("hSpec");
        //ctx.print();
        //spectraCh0[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //spectraCh0[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        spectraCh0[i]->GetYaxis()->SetMaxDigits(2);
        spectraCh0[i]->Draw();
        
        fun_name = Form("f_S%i_ch0_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]); 
        TF1 *fun_tmp_ch0 = spectraCh0[i]->GetFunction(fun_name);
        
        if (fun_tmp_ch0)
        {
            text.DrawLatex(0.65, 0.75, Form("#chi^{2}/NDF = %.3f", 
                           fun_tmp_ch0->GetChisquare() / fun_tmp_ch0->GetNDF()));
            text.DrawLatex(0.65, 0.70, Form("#mu = %.3f +/- %.3f", 
                           fun_tmp_ch0->GetParameter(1), fun_tmp_ch0->GetParError(1)));
            text.DrawLatex(0.65, 0.65, Form("#sigma = %.3f +/- %.3f", 
                           fun_tmp_ch0->GetParameter(2), fun_tmp_ch0->GetParError(2)));
        }

        can_spectra_ch1->cd(i + 1);
        gPad->SetGrid(1, 1);
        spectraCh1[i]->SetStats(false);
        spectraCh1[i]->GetXaxis()->SetRangeUser(0, 1200);
        spectraCh1[i]->SetTitle(Form("PE spectrum S%i Ch1, source position %.2f mm",
                                     seriesNo, positions[i]));
        spectraCh1[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
        spectraCh1[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        spectraCh1[i]->GetXaxis()->SetTitle("charge [P.E.]");
        spectraCh1[i]->GetYaxis()->SetTitle("counts");
        //spectraCh1[i]->GetXaxis()->SetRangeUser(ctx.x.min, ctx.x.max);
        //spectraCh1[i]->GetYaxis()->SetRangeUser(ctx.y.min, ctx.y.max);
        spectraCh1[i]->GetYaxis()->SetMaxDigits(2);
        spectraCh1[i]->Draw();
        
        fun_name = Form("f_S%i_ch1_pos%.1f_ID%i_PE", seriesNo, positions[i], ID[i]); 
        TF1* fun_tmp_ch1 = spectraCh1[i]->GetFunction(fun_name);
        
        if (fun_tmp_ch1)
        {
            text.DrawLatex(0.65, 0.75, Form("#chi^{2}/NDF = %.3f", 
                           fun_tmp_ch1->GetChisquare() / fun_tmp_ch1->GetNDF()));
            text.DrawLatex(0.65, 0.70, Form("#mu = %.3f +/- %.3f", 
                           fun_tmp_ch1->GetParameter(1), fun_tmp_ch1->GetParError(1)));
            text.DrawLatex(0.65, 0.65, Form("#sigma = %.3f +/- %.3f", 
                           fun_tmp_ch1->GetParameter(2), fun_tmp_ch1->GetParError(2)));
        }
    }

    //----- saving
    TString fname       = Form("attenuation_series%i.root", seriesNo);
    TString fname_full  = outdir + fname;
    TString dbname_full = outdir + dbase;

    TFile* file = new TFile(fname_full, "RECREATE");

    if (!file->IsOpen())
    {
        std::cerr << "##### Error in attenuation.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << fname_full << std::endl;
        return 1;
    }

    can_averaged_ch->Write();
    can_sigma->Write();
    can_separate_ch->Write();
    can_ratios->Write();
    can_spectra_ch0->Write();
    can_spectra_ch1->Write();
    file->Close();

    //----- writing results to the data base
    TString table = "ATTENUATION_LENGTH";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ATT_CH0, "
                         "ATT_CH0_ERR, CHI2NDF_CH0, ATT_CH1, ATT_CH1_ERR, CHI2NDF_CH1, "
                         "ATT_COMB, ATT_COMB_ERR, CHI2NDF_COMB, ATT_COMB_POL3, ATT_COMB_POL3_ERR, "
                         "CHI2NDF_POL3, ATT_SIM, ATT_SIM_ERR, CHI2NDF_SIM) VALUES "
                         "(%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, "
                         "%f, %f, %f)", table.Data(), seriesNo, fname_full.Data(), 
                         results[0]->GetValue(SFResultTypeNum::kLambda),
                         results[0]->GetUncertainty(SFResultTypeNum::kLambda), 
                         results[0]->GetValue(SFResultTypeNum::kChi2NDF),
                         results[1]->GetValue(SFResultTypeNum::kLambda),
                         results[1]->GetUncertainty(SFResultTypeNum::kLambda),
                         results[1]->GetValue(SFResultTypeNum::kChi2NDF),
                         results[2]->GetValue(SFResultTypeNum::kLambda),
                         results[2]->GetUncertainty(SFResultTypeNum::kLambda),
                         results[2]->GetValue(SFResultTypeNum::kChi2NDF),
                         results[3]->GetValue(SFResultTypeNum::kLambda),
                         results[3]->GetValue(SFResultTypeNum::kLambda),
                         results[3]->GetValue(SFResultTypeNum::kChi2NDF),
                         results[4]->GetValue(SFResultTypeNum::kLambda),
                         results[4]->GetUncertainty(SFResultTypeNum::kLambda),
                         results[4]->GetValue(SFResultTypeNum::kChi2NDF));
     

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- attenuation writing, try number " << (max_tries - i_try) + 1
                  << std::endl;
        stat = SFTools::SaveResultsDB(dbname_full, table, query, seriesNo);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- attenuation writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);

    delete data;
    delete att;

    return 0;
}
