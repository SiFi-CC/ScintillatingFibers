// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *          attenuation_krk.cc           *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFInfo.hh"
#include "SFPeakFinder.hh"
#include "SFTools.hh"
#include "SFAttenuation.hh"
#include "SFAttenuationModel.hh"
#include "common_options_krk.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>

#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{
    //----- parsing arguments
    
    TString log_path;
    TString data_path;
    TString outdir;
    TString dbase;
    int m = 0;
    int l =0;
    int f = 0;
    
    if (argc < 2)
    {
        std::cerr << "to run type: ./attenuation_krk path/to/data/log.txt ";
        std::cerr << "-m 0 -l 0 -f 0 -out path/to/results/ ";
        std::cerr << "-db path/to/database.db" << std::endl;
        return 1;
    }
    
    int ret = parse_common_options(argc, argv, log_path, data_path, outdir, dbase, m, l, f);
    if (ret != 0)
        exit(ret);
    
    TString tmp    = TString(gSystem->DirName(outdir));
    char* abs_path = realpath((std::string(tmp)).c_str(), NULL);
    
    //----- setting measurement properties
    
    std::ifstream logfile(log_path);
    
    if(!logfile.is_open())
    {
        std::cerr << "##### Error in attenuation_krk!" << std::endl;
        std::cerr << "Couldn't open logfile: " << log_path << std::endl;
        return 1;
    }
    
    int npoints = 0;
    std::vector<TString> meas_names;
    std::vector<double> positions;
    std::vector<int> times;
    std::vector<int> start;
    std::vector<int> stop;
    
    std::string line = " ";
    TString name_tmp = "";
    double pos_tmp   = 0;
    double time_tmp  = 0;
    double start_tmp = 0;
    double stop_tmp  = 0;
    
    while (true)
    {
        getline(logfile, line);
        logfile >> line >> line >> start_tmp;
        logfile >> line >> line >> stop_tmp;
        logfile >> name_tmp;
        logfile >> line >> line >> pos_tmp >> line;
        getline(logfile, line);
        getline(logfile, line);
 
        time_tmp = ((stop_tmp - start_tmp)/60) - 1;
        
        std::cout << name_tmp << "\t" << pos_tmp << " mm \t" << time_tmp << " min" << std::endl;
        meas_names.push_back(data_path + name_tmp);
        positions.push_back(pos_tmp);   // milimeters
        times.push_back(time_tmp * 60); // seconds        
        npoints++;
        
        if(logfile.eof())
            break;
    }
    
    std::cout << "\nnpoints: " << npoints << std::endl;
    logfile.close();
    
    SFInfo *info = new SFInfo();
    info->SetNpoints(npoints);
    info->SetFiberLength(100);
    info->SetFiber("LYSO:Ce");
    info->SetCollimator("Electronic long");
    info->SetTestBench("PRO KRK");
    info->SetSiPM("Ketek");
    info->SetLogFile(log_path);
    info->SetNames(meas_names);
    info->SetPositions(positions);
    info->SetTimes(times);
    info->SetStartTimes(start);
    info->SetStopTimes(stop);
    
    double pos_uncert = 2.; // mm // TODO check
    
    //----- getting histograms
        
    std::vector<TH1D*> hQL;
    std::vector<TH1D*> hQR;
    std::vector<TH1D*> hMLR;
    
    for (int i=0; i<npoints; i++)
    {
        TString full_root_path = meas_names[i] + "/sifi_results_HISTOS.root";
        TString dir_name = Form("/Module%i/Layer%i/Fiber%i/", m, l, f);
        std::cout << "\nOpening file: " << full_root_path << std::endl;
        
        auto file = std::make_unique<TFile>(full_root_path, "READ");
        
        if (!file.get()->IsOpen())
        {
            std::cerr << "##### Error in attenuation_krk! Could not open file: " << std::endl;
            std::cerr << full_root_path << std::endl;
            std::abort();
        }
        
        //----- spectrum - left
        TString hname = Form("Q_M%iL%iF%iL", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);
        auto htemp_L = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_L)
        {
            std::cerr << "##### Error in attenuation_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_L = (TH1D*)htemp_L->Clone(hname);
        
        if (!clone_L)
        {
            std::cerr << "##### Error in attenuation_pmi! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }
        
        clone_L->SetDirectory(nullptr);
        hQL.push_back(clone_L);
        
        //----- spectrum - right
        hname = Form("Q_M%iL%iF%iR", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);

        auto htemp_R = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_R)
        {
            std::cerr << "##### Error in attenuation_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_R = (TH1D*)htemp_R->Clone(hname);
        
        if (!clone_R)
        {
            std::cerr << "##### Error in attenuation_pmi! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_R->SetDirectory(nullptr);
        hQR.push_back(clone_R);
        
        //----- spectrum - average
        hname = Form("MLR_M%iL%iF%i", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);

        auto htemp_MLR = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_MLR) 
        {
            std::cerr << "##### Error in attenuation_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_MLR = (TH1D*)htemp_MLR->Clone(hname);
        
        if (!clone_MLR)
        {
            std::cerr << "##### Error in attenuation_pmi! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_MLR->SetDirectory(nullptr);
        hMLR.push_back(clone_MLR);
    }
    
    //----- attenuation analysis MLR
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.03);
    
    SFResults *results_mlr = SFAttenuation::AttCombinedCh("pol1", info, pos_uncert,
                                                          positions, hMLR);
    
    TGraphErrors* gMLR = (TGraphErrors*)results_mlr->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors* gSigma = (TGraphErrors*)results_mlr->GetObject(SFResultTypeObj::kMLRSigmaGraph);
    
    TCanvas *can_mlr = new TCanvas("can_mlr", "can_mlr", 1200, 600);
    can_mlr->Divide(2, 1);
    
    can_mlr->cd(1);
    gPad->SetGrid(1, 1); 
    gMLR->SetTitle(Form("M_{LR} vs. source position M%iL%iF%i", m, l, f));
    gMLR->GetYaxis()->SetTitleOffset(1.4);
    gMLR->GetXaxis()->SetLabelSize(0.035);
    gMLR->GetYaxis()->SetLabelSize(0.035);
    gMLR->Draw("AP");
    text.DrawLatex(0.5, 0.85, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.5, 0.80, Form("p_{0} = %.3e +/- %.3e", 
                   gMLR->GetFunction("fpol1")->GetParameter(0),
                   gMLR->GetFunction("fpol1")->GetParError(0)));
    text.DrawLatex(0.5, 0.75, Form("p_{1} = %.3e +/- %.3e", 
                   gMLR->GetFunction("fpol1")->GetParameter(1),
                   gMLR->GetFunction("fpol1")->GetParError(1)));
    text.DrawLatex(0.5, 0.70, Form("L_{att} = (%.2f +/- %.2f) mm",
                   results_mlr->GetValue(SFResultTypeNum::kLambda),
                   results_mlr->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.5, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   results_mlr->GetValue(SFResultTypeNum::kChi2NDF)));
    
    can_mlr->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetGrid(1, 1);
    
    gSigma->SetTitle(Form("#sigma M_{LR} vs. source position M%iL%iF%i", m, l, f));
    gSigma->Draw("AP");
    
    TCanvas *can_ratios = new TCanvas("can_ratios", "can_ratios", 2000, 1200);
    can_ratios->DivideSquare(npoints);
    
    TCanvas *can_l = new TCanvas("can_l", "can_l", 2000, 1200);
    can_l->DivideSquare(npoints);
    
    TCanvas *can_r = new TCanvas("can_r", "can_r", 2000, 1200);
    can_r->DivideSquare(npoints);
    
    for (int i=0; i<npoints; i++)
    {
        can_ratios->cd(i+1);
        gPad->SetGrid(1, 1);
        hMLR[i]->Draw();
        
        can_l->cd(i+1);
        gPad->SetGrid(1, 1);
        hQL[i]->Draw();
        
        can_r->cd(i+1);
        gPad->SetGrid(1, 1);
        hQR[i]->Draw();
    }
    
    //----- Attenuation analysis ELA (separate fit)
    int col_l = kPink-8;
    int col_r = kAzure-6;
    
    SFResults *results_L = SFAttenuation::AttSeparateCh('l', pos_uncert, info,
                                                        positions, hQL, meas_names);
    
    SFResults *results_R = SFAttenuation::AttSeparateCh('r', pos_uncert, info,
                                                        positions, hQR, meas_names);
    
    TGraphErrors *gL = (TGraphErrors*)results_L->GetObject(SFResultTypeObj::kAttGraph);
    TGraphErrors *gR = (TGraphErrors*)results_R->GetObject(SFResultTypeObj::kAttGraph);
    
    gL->SetMarkerColor(col_l);
    gL->SetLineColor(col_l);
    
    gR->SetMarkerColor(col_r);
    gR->SetLineColor(col_r);
    
    bool valid_ela_l = true;
    double lambda_ela_l = 0;
    double lambda_ela_l_err = 0;
    double chi2ndf_ela_l = 0;
    
    if(results_L->GetValue(SFResultTypeNum::kLambda) != 
       results_L->GetValue(SFResultTypeNum::kLambda) ||
       results_L->GetUncertainty(SFResultTypeNum::kLambda) !=
       results_L->GetUncertainty(SFResultTypeNum::kLambda) ||
       results_L->GetValue(SFResultTypeNum::kLambda) < 0 ||
       results_L->GetValue(SFResultTypeNum::kLambda) > 5000)
    {
        std::cout << "\n\nWarning in attenuation_krk!" << std::endl;
        std::cout << "ELA model (left) - invalid results!\n\n" << std::endl;
        valid_ela_l = false;
        lambda_ela_l = -1;
        lambda_ela_l_err = -1;
        chi2ndf_ela_l = -1;
    }
    else
    {
        lambda_ela_l = results_L->GetValue(SFResultTypeNum::kLambda);
        lambda_ela_l_err = results_L->GetUncertainty(SFResultTypeNum::kLambda);
        chi2ndf_ela_l = results_L->GetValue(SFResultTypeNum::kChi2NDF);
    }
        
    
    bool valid_ela_r = true;
    double lambda_ela_r = 0;
    double lambda_ela_r_err = 0;
    double chi2ndf_ela_r = 0;
    
    if(results_R->GetValue(SFResultTypeNum::kLambda) != 
       results_R->GetValue(SFResultTypeNum::kLambda) ||
       results_R->GetUncertainty(SFResultTypeNum::kLambda) !=
       results_R->GetUncertainty(SFResultTypeNum::kLambda) ||
       results_R->GetValue(SFResultTypeNum::kLambda) < 0 ||
       results_R->GetValue(SFResultTypeNum::kLambda) > 5000)
    {
        std::cout << "Warning in attenuation_krk!" << std::endl;
        std::cout << "ELA model (right) - invalid results!" << std::endl;
        valid_ela_r = false;
        lambda_ela_r = -1;
        lambda_ela_r = -1;
        chi2ndf_ela_r = -1;
    }
    else
    {
        lambda_ela_r = results_R->GetValue(SFResultTypeNum::kLambda);
        lambda_ela_r_err = results_R->GetUncertainty(SFResultTypeNum::kLambda);
        chi2ndf_ela_r = results_R->GetValue(SFResultTypeNum::kChi2NDF);
    }
    
    //------ Attenuation analysis ELA (simultaneous fit)
    
    TGraphErrors *gL_sim = (TGraphErrors*)gL->Clone(TString(gL->GetName())+"_sim");
   
    TF1 *tmp_fun_1 = gL_sim->GetFunction("fun_l");
    if(tmp_fun_1)
        gL_sim->GetFunction("fun_l")->Delete();
    
    TGraphErrors *gR_sim = (TGraphErrors*)gR->Clone(TString(gR->GetName())+"_sim");
    
    TF1 *tmp_fun_2 = gR_sim->GetFunction("fun_r");
    if(tmp_fun_2)
        gR_sim->GetFunction("fun_r")->Delete();
    
    SFResults *results_sim = SFAttenuation::AttSimultaneousFit(gL_sim, gR_sim);
    
    bool valid_sim = true;
    double lambda_sim = 0;
    double lambda_sim_err = 0;
    double chi2ndf_sim = 0;
    
    if(results_sim->GetValue(SFResultTypeNum::kLambda) != 
       results_sim->GetValue(SFResultTypeNum::kLambda) ||
       results_sim->GetUncertainty(SFResultTypeNum::kLambda) !=
       results_sim->GetUncertainty(SFResultTypeNum::kLambda) ||
       results_sim->GetValue(SFResultTypeNum::kLambda) < 0 ||
       results_sim->GetValue(SFResultTypeNum::kLambda) > 5000)
    {
        std::cout << "Warning in attenuation_krk!" << std::endl;
        std::cout << "ELA model (simultaneous) - invalid results!" << std::endl;
        valid_sim = false;
        lambda_sim = -1;
        lambda_sim_err = -1;
        chi2ndf_sim = -1;
    }
    else
    {
        lambda_sim = results_sim->GetValue(SFResultTypeNum::kLambda);
        lambda_sim_err = results_sim->GetUncertainty(SFResultTypeNum::kLambda);
        chi2ndf_sim = results_sim->GetValue(SFResultTypeNum::kChi2NDF);
    }
    
    //----- ELA drawing
    
    TString hname = Form("ELA Separate Fitting M%iL%iF%i", m, l, f);
    TH1D *htmp = new TH1D("h", hname, 100, 0, 100);
    htmp->GetXaxis()->SetTitle("source position [mm]");
    htmp->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    htmp->SetStats(false);
    
    hname = Form("ELA Simulataneus Fitting M%iL%iF%i", m, l, f);
    TH1D *htmp_2 = new TH1D("h2", hname, 100, 0, 100);
    htmp_2->GetXaxis()->SetTitle("source position [mm]");
    htmp_2->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    htmp_2->SetStats(false);
    
    TCanvas *can_ela = new TCanvas("can_ela", "can_ela", 1200, 600);
    can_ela->Divide(2, 1);
    
    can_ela->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.01);
    gPad->SetGrid(1, 1);
    
    htmp->Draw();
    gL->Draw("P");
    gL->GetFunction("fun_l")->SetLineColor(col_l);
    text.SetTextColor(col_l);
    if(!valid_ela_l)
        text.DrawLatex(0.3, 0.85, "INVALID");
    text.DrawLatex(0.3, 0.8, Form("L_{att l} = (%.2f +/- %.2f) mm",
                   results_L->GetValue(SFResultTypeNum::kLambda),
                   results_L->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.75, Form("#chi^{2}/NDF_{l} = %.2f", 
                   results_L->GetValue(SFResultTypeNum::kChi2NDF)));
    
    gR->Draw("P");
    gR->GetFunction("fun_r")->SetLineColor(col_r);
    text.SetTextColor(col_r);
    if(!valid_ela_r)
        text.DrawLatex(0.3, 0.65, "INVALID");
    text.DrawLatex(0.3, 0.60, Form("L_{att r} = (%.2f +/- %.2f) mm",
                   results_R->GetValue(SFResultTypeNum::kLambda),
                   results_R->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.55, Form("#chi^{2}/NDF_{r} = %.2f", 
                   results_R->GetValue(SFResultTypeNum::kChi2NDF)));
    
//     double min_l = gL->GetMinimum();
//     double max_l = gL->GetMaximum();
//     double min_r = gR->GetMinimum();
//     double max_r = gR->GetMaximum();
//     
//     double min = min_l < min_r ? min_l : min_r;
//     double max = max_l > max_r ? max_l : max_r;
//     
//     htmp->GetYaxis()->SetRangeUser(min-20, max+20);
    
    htmp->GetYaxis()->SetRangeUser(0, 2000);
    
    can_ela->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.01);
    gPad->SetGrid(1, 1);
    
    htmp_2->Draw();
    htmp_2->GetYaxis()->SetRangeUser(0, 2000);
    gL_sim->Draw("P");
    gR_sim->Draw("P");
    
    text.SetTextColor(kBlack);
    if(!valid_sim)
        text.DrawLatex(0.3, 0.80, "INVALID");
    text.DrawLatex(0.3, 0.75, Form("L_{att} = (%.2f +/- %.2f) mm",
                   results_sim->GetValue(SFResultTypeNum::kLambda),
                   results_sim->GetUncertainty(SFResultTypeNum::kLambda)));
    text.DrawLatex(0.3, 0.70, Form("#chi^{2}/NDF_{1} = %.2f", 
                   results_sim->GetValue(SFResultTypeNum::kChi2NDF)));
    
    //----- attenuation analysis ELAR
    
    TGraphErrors *gL_elar = (TGraphErrors*)gL->Clone(TString(gL->GetName())+"_elar");
    
    TF1* tmp_fun_3 = gL_elar->GetFunction("fun_l");
    if(tmp_fun_3)
        gL_elar->GetFunction("fun_l")->Delete();
    
    TGraphErrors *gR_elar = (TGraphErrors*)gR->Clone(TString(gR->GetName())+"_elar");
    
    TF1* tmp_fun_4 = gR_elar->GetFunction("fun_r");
    if(tmp_fun_4)
        gR_elar->GetFunction("fun_r")->Delete();
    
    SFResults *results_elar = SFAttenuationModel::FitModel(gL_elar, gR_elar, pos_uncert, 
                                                           info->GetFiberLength());
    
    TF1* fun_Sl = (TF1*)results_elar->GetObject(SFResultTypeObj::kSlFun);
    TF1* fun_Sr = (TF1*)results_elar->GetObject(SFResultTypeObj::kSrFun);
    TF1* fun_Pl = (TF1*)results_elar->GetObject(SFResultTypeObj::kPlFun);
    TF1* fun_Pr = (TF1*)results_elar->GetObject(SFResultTypeObj::kPrFun);
    TF1* fun_Rl = (TF1*)results_elar->GetObject(SFResultTypeObj::kRlFun);
    TF1* fun_Rr = (TF1*)results_elar->GetObject(SFResultTypeObj::kRrFun);
    
    TGraphErrors *gL_corr = (TGraphErrors*)results_elar->GetObject(SFResultTypeObj::kPlVsPosGraph);
    TGraphErrors *gR_corr = (TGraphErrors*)results_elar->GetObject(SFResultTypeObj::kPrVsPosGraph);
    
    TCanvas *can_elar = new TCanvas("can_elar", "can_elar", 1200, 600);
    can_elar->Divide(2, 1);
       
    can_elar->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetGrid(1, 1);
    
    hname = Form("ELAR Attenuation Curves M%iL%iF%i", m, l, f);
    TH1D *htmp_3 = new TH1D("h3", hname, 100, 0, 100);
    htmp_3->GetXaxis()->SetTitle("source position [mm]");
    htmp_3->GetYaxis()->SetTitle("511 keV peak position [P.E.]");
    htmp_3->SetStats(false);
//     htmp_3->GetYaxis()->SetRangeUser(min-20, max+20);
    htmp_3->GetYaxis()->SetRangeUser(0, 2000);
    htmp_3->Draw();
    
    gL_corr->SetMarkerColor(col_l);
    gL_corr->SetLineColor(col_l);
    gR_corr->SetMarkerColor(col_r);
    gR_corr->SetLineColor(col_r);
    gL_elar->Draw("P");
    gR_elar->Draw("P");
    gL_corr->Draw("P");
    gR_corr->Draw("P");
    
    fun_Sl->SetLineColor(col_l);
    fun_Sr->SetLineColor(col_r);
    fun_Pl->SetLineColor(col_l);
    fun_Pl->SetLineStyle(2);
    fun_Pr->SetLineColor(col_r);
    fun_Pr->SetLineStyle(2);
    fun_Rl->SetLineColor(col_l);
    fun_Rl->SetLineStyle(9);
    fun_Rr->SetLineColor(col_r);
    fun_Rr->SetLineStyle(9);
    
    fun_Sl->Draw("same");
    fun_Sr->Draw("same");
    fun_Pl->Draw("same");
    fun_Pr->Draw("same");
    fun_Rl->Draw("same");
    fun_Rr->Draw("same");
    
    can_elar->cd(2);

    TLegend* leg = new TLegend(0.1, 0.4, 0.9, 0.9);
    leg->AddEntry(gL_elar, "experimental data L", "PE");
    leg->AddEntry(gR_elar, "experimental data R", "PE");
    leg->AddEntry(gL_corr, "recalculated data L", "PE");
    leg->AddEntry(gR_corr, "recalculated data R", "PE");
    leg->AddEntry(fun_Pl, "direct light L", "L");
    leg->AddEntry(fun_Pr, "direct light R", "L");
    leg->AddEntry(fun_Rl, "reflected light L", "L");
    leg->AddEntry(fun_Rr, "reflected light R", "L");
    leg->AddEntry(fun_Sl, "total signal L", "L");
    leg->AddEntry(fun_Sr, "total signal R", "L");
    leg->Draw();

    text.DrawLatex(0.2, 0.35, Form("S_{0} = %.2f +/- %.2f",
                   results_elar->GetValue(SFResultTypeNum::kS0),
                   results_elar->GetUncertainty(SFResultTypeNum::kS0)));
    
    text.DrawLatex(0.2, 0.30, Form("#lambda = %.2f +/- %.2f mm", 
                   results_elar->GetValue(SFResultTypeNum::kLambda),
                   results_elar->GetUncertainty(SFResultTypeNum::kLambda)));
    
    if (results_elar->GetValue(SFResultTypeNum::kEtaR) < 0 ||
        results_elar->GetValue(SFResultTypeNum::kEtaR) > 1)
        text.SetTextColor(kRed);
    
    text.DrawLatex(0.2, 0.25, Form("#eta_{R} = %.3f +/- %.3f", 
                   results_elar->GetValue(SFResultTypeNum::kEtaR),
                   results_elar->GetUncertainty(SFResultTypeNum::kEtaR)));
    text.SetTextColor(kBlack);
    
    if (results_elar->GetValue(SFResultTypeNum::kEtaL) < 0 ||
        results_elar->GetValue(SFResultTypeNum::kEtaL) > 1)
        text.SetTextColor(kRed);
    
    text.DrawLatex(0.2, 0.20, Form("#eta_{L} = %.3f +/- %.3f",
                   results_elar->GetValue(SFResultTypeNum::kEtaL),
                   results_elar->GetUncertainty(SFResultTypeNum::kEtaL)));
    text.SetTextColor(kBlack);
    
    text.DrawLatex(0.2, 0.15, Form("#xi = %.3f +/- %.3f",
                   results_elar->GetValue(SFResultTypeNum::kKsi),
                   results_elar->GetUncertainty(SFResultTypeNum::kKsi)));
    
    text.DrawLatex(0.2, 0.10, Form("L = %.2f mm (fixed)",
                   results_elar->GetValue(SFResultTypeNum::kLength)));
    
    text.DrawLatex(0.2, 0.05, Form("#chi^{2}/NDF = %.3f",
                   results_elar->GetValue(SFResultTypeNum::kChi2NDF)));
    
    //----- saving to ROOT file
    TString fout_name = Form("/attenuation_M%iL%iF%i.root", m, l, f);
    TFile *out_file = new TFile(std::string(abs_path) + fout_name, "RECREATE");
    
    if (!out_file->IsOpen())
    {
        std::cerr << "##### Error in attenuation_krk.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << std::string(abs_path) + fout_name << std::endl;
        return 1;
    }
    
    can_mlr->Write();
    can_ela->Write();
    can_elar->Write();
    can_ratios->Write();
    can_l->Write();
    can_r->Write();
    
    out_file->mkdir("histograms");
    out_file->cd("histograms");
    
    for (int i=0; i<npoints; i++)
    {
        hQL[i]->Write();
        hQR[i]->Write();
        hMLR[i]->Write();
    }
    
    out_file->Close();
    
    //----- saving to the database
    
    //----- saving ELA
    TString address_str = Form("%i%i%i", m, l, f);
    int address = atoi(address_str);
    
    TString table = "ATTENUATION_LENGTH";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ATT_CH0, "
                         "ATT_CH0_ERR, CHI2NDF_CH0, ATT_CH1, ATT_CH1_ERR, CHI2NDF_CH1, "
                         "ATT_COMB, ATT_COMB_ERR, CHI2NDF_COMB, ATT_COMB_POL3, ATT_COMB_POL3_ERR, "
                         "CHI2NDF_POL3, ATT_SIM, ATT_SIM_ERR, CHI2NDF_SIM) VALUES "
                         "(%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, "
                         "%f, %f, %f)", table.Data(), address, 
                         (std::string(abs_path) + fout_name).Data(), 
                         lambda_ela_l, lambda_ela_l_err, chi2ndf_ela_l,
                         lambda_ela_r, lambda_ela_r_err, chi2ndf_ela_r,
                         results_mlr->GetValue(SFResultTypeNum::kLambda),
                         results_mlr->GetUncertainty(SFResultTypeNum::kLambda),
                         results_mlr->GetValue(SFResultTypeNum::kChi2NDF),
                         -1.0, -1.0, -1.0,
                         lambda_sim, lambda_sim_err, chi2ndf_sim);
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- attenuation_krk writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbase, table, query, address);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- attenuation writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    //----- saving ELAR
    
    table = "ATTENUATION_MODEL";
    query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, S0, S0_ERR, "
                 "LAMBDA, LAMBDA_ERR, ETAR, ETAR_ERR, ETAL, ETAL_ERR, KSI, KSI_ERR, CHI2NDF) "
                 "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)",
                 table.Data(), address, (std::string(abs_path) + fout_name).Data(),
                 results_elar->GetValue(SFResultTypeNum::kS0),
                 results_elar->GetUncertainty(SFResultTypeNum::kS0),
                 results_elar->GetValue(SFResultTypeNum::kLambda),
                 results_elar->GetUncertainty(SFResultTypeNum::kLambda),
                 results_elar->GetValue(SFResultTypeNum::kEtaR),
                 results_elar->GetUncertainty(SFResultTypeNum::kEtaR),
                 results_elar->GetValue(SFResultTypeNum::kEtaL),
                 results_elar->GetUncertainty(SFResultTypeNum::kEtaL),
                 results_elar->GetValue(SFResultTypeNum::kKsi),
                 results_elar->GetUncertainty(SFResultTypeNum::kKsi),
                 results_elar->GetValue(SFResultTypeNum::kChi2NDF));
    
    i_try = max_tries;
    wait  = 0;
    stat  = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- attenuation_krk writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbase, table, query, address);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- attenuation_krk writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);    
    return 0;
}