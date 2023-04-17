// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *           energyres_krk.cc            *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2022              *
// *                                       *
// *****************************************

#include "SFInfo.hh"
#include "SFPeakFinder.hh"
#include "SFTools.hh"
#include "SFEnergyRes.hh"
#include "common_options_krk.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>

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
        std::cerr << "to run type: ./energyres_krk path/to/data/log.txt ";
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
        std::cerr << "##### Error in energyres_krk!" << std::endl;
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
    std::vector<TH1D*> hQAve;
    
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
            std::cerr << "##### Error in attenuation_krk! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_R->SetDirectory(nullptr);
        hQR.push_back(clone_R);
        
        //----- spectrum - average
        hname = Form("QAve_M%iL%iF%i", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);

        auto htemp_QAve = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_QAve) 
        {
            std::cerr << "##### Error in energyres_krk! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_QAve = (TH1D*)htemp_QAve->Clone(hname);
        
        if (!clone_QAve)
        {
            std::cerr << "##### Error in energyres_krk! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_QAve->SetDirectory(nullptr);
        hQAve.push_back(clone_QAve);
    }
        
    //----- energy resolution analysis
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.04);
    
    TString suffix = Form("M%iL%iF%i", m, l, f);
    
    TCanvas *can_ave = new TCanvas("er_ave", "er_ave", 700, 500);
    gPad->SetGrid(1, 1);
    
    SFResults *res_ave = SFEnergyRes::CalculateEnergyResSeries(suffix, pos_uncert,
                                                               positions, hQAve,
                                                               meas_names);
    
    TGraphErrors *gERAve = (TGraphErrors*)res_ave->GetObject(SFResultTypeObj::kEnergyResGraph);
    double er_ave = res_ave->GetValue(SFResultTypeNum::kEnergyRes);
    double er_ave_err = res_ave->GetUncertainty(SFResultTypeNum::kEnergyRes);
    
    gERAve->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                  er_ave, er_ave_err));
    
    TCanvas *can_lr = new TCanvas("er_lr", "er_lr", 1400, 500);
    can_lr->Divide(2, 1);
    
    can_lr->cd(1);
    gPad->SetGrid(1, 1);
    
    SFResults *res_l = SFEnergyRes::CalculateEnergyResSeries(suffix, pos_uncert,
                                                             positions, hQL,
                                                             meas_names);
    
    TGraphErrors *gERL = (TGraphErrors*)res_l->GetObject(SFResultTypeObj::kEnergyResGraph);
    double er_l = res_l->GetValue(SFResultTypeNum::kEnergyRes);
    double er_l_err = res_l->GetUncertainty(SFResultTypeNum::kEnergyRes);
    
    gERL->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                  er_l, er_l_err));
    
    can_lr->cd(2);
    gPad->SetGrid(1, 1);
    
    SFResults *res_r = SFEnergyRes::CalculateEnergyResSeries(suffix, pos_uncert,
                                                             positions, hQR,
                                                             meas_names);
    
    TGraphErrors *gERR = (TGraphErrors*)res_r->GetObject(SFResultTypeObj::kEnergyResGraph);
    double er_r = res_r->GetValue(SFResultTypeNum::kEnergyRes);
    double er_r_err = res_r->GetUncertainty(SFResultTypeNum::kEnergyRes);
    
    gERR->Draw("AP");
    text.DrawLatex(0.2, 0.8, Form("ER = (%.2f +/- %.2f) %%", 
                  er_r, er_r_err));
    
    //----- drawing histograms
        
    TCanvas *can_hl = new TCanvas("er_spec_l", "er_spec_l", 2000, 1200);
    can_hl->DivideSquare(npoints);
    
    TCanvas *can_hr = new TCanvas("er_spec_r", "er_spec_r", 2000, 1200);
    can_hr->DivideSquare(npoints);
    
    TCanvas *can_have = new TCanvas("er_spec_ave", "er_spec_ave", 2000, 1200);
    can_have->DivideSquare(npoints);
    
    TString fun_name = "";
    
    for(int i=0; i<npoints; i++)
    {
        can_have->cd(i + 1);
        gPad->SetGrid(1, 1);
        hQAve[i]->SetStats(false);
        hQAve[i]->GetXaxis()->SetTitle("charge [a.u.]");
        hQAve[i]->GetYaxis()->SetTitle("counts");
        hQAve[i]->SetTitle(Form("Average PE spectrum, position %.2f mm", positions[i]));
//         hQAve[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
//         hQAve[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_av);
        hQAve[i]->Draw();
        fun_name = "f_" + TString(hQAve[i]->GetName());
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gERAve->GetPointY(i), gERAve->GetErrorY(i)));
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       hQAve[i]->GetFunction(fun_name)->GetParameter(1),
                       hQAve[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       hQAve[i]->GetFunction(fun_name)->GetParameter(2),
                       hQAve[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       hQAve[i]->GetFunction(fun_name)->GetChisquare() / 
                       hQAve[i]->GetFunction(fun_name)->GetNDF()));

        can_hl->cd(i + 1);
        gPad->SetGrid(1, 1);
        hQL[i]->SetStats(false);
        hQL[i]->GetXaxis()->SetTitle("charge [a.u.]");
        hQL[i]->GetYaxis()->SetTitle("counts");
        hQL[i]->SetTitle(Form("Charge spectrum L, position %.2f mm", positions[i]));
//         hQL[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
//         hQL[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_0);
        hQL[i]->Draw();
        fun_name = "f_" + TString(hQL[i]->GetName());
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gERL->GetPointY(i), gERL->GetErrorY(i)));
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       hQL[i]->GetFunction(fun_name)->GetParameter(1),
                       hQL[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       hQL[i]->GetFunction(fun_name)->GetParameter(2),
                       hQL[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       hQL[i]->GetFunction(fun_name)->GetChisquare() / 
                       hQL[i]->GetFunction(fun_name)->GetNDF()));

        can_hr->cd(i + 1);
        gPad->SetGrid(1, 1);
        hQR[i]->SetStats(false);
        hQR[i]->GetXaxis()->SetTitle("charge [P.E.]");
        hQR[i]->GetYaxis()->SetTitle("counts");
        hQR[i]->SetTitle(Form("Charge spectrum R, position %.2f mm", positions[i]));
//         hQR[i]->GetXaxis()->SetRangeUser(min_xaxis, max_xaxis);
//         hQR[i]->GetYaxis()->SetRangeUser(min_yaxis, max_yaxis_1);
        hQR[i]->Draw();
        fun_name = "f_" + TString(hQR[i]->GetName());
        text.DrawLatex(0.6, 0.80, Form("ER = (%.2f +/- %.2f) %%", 
                       gERR->GetPointY(i), gERR->GetErrorY(i)));
        text.DrawLatex(0.6, 0.75, Form("#mu = %.2f +/- %.2f", 
                       hQR[i]->GetFunction(fun_name)->GetParameter(1),
                       hQR[i]->GetFunction(fun_name)->GetParError(1)));
        text.DrawLatex(0.6, 0.70, Form("#sigma = %.2f +/- %.2f", 
                       hQR[i]->GetFunction(fun_name)->GetParameter(2),
                       hQR[i]->GetFunction(fun_name)->GetParError(2)));
        text.DrawLatex(0.6, 0.60, Form("#chi^{2}/NDF = %.3f",
                       hQR[i]->GetFunction(fun_name)->GetChisquare() / 
                       hQR[i]->GetFunction(fun_name)->GetNDF()));
    }
    
    //----- saving to ROOT file

    TString fout_name = Form("/energyres_M%iL%iF%i.root", m, l, f);
    TFile *out_file = new TFile(std::string(abs_path) + fout_name, "RECREATE");
    
    if (!out_file->IsOpen())
    {
        std::cerr << "##### Error in energyres_krk.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << std::string(abs_path) + fout_name << std::endl;
        return 1;
    }
    
    can_ave->Write();
    can_lr->Write();
    can_hl->Write();
    can_hr->Write();
    can_have->Write();
    
    out_file->mkdir("histograms");
    out_file->cd("histograms");
    
    for (int i=0; i<npoints; i++)
    {
        hQL[i]->Write();
        hQR[i]->Write();
        hQAve[i]->Write();
    }
    
    out_file->Close();    
    
    //----- saving to the data base
    
    TString address_str = Form("%i%i%i", m, l, f);
    int address = atoi(address_str);
        
    TString table = "ENERGY_RESOLUTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, ENRES_AV, "
                         "ENRES_AV_ERR, ENRES_CH0, ENRES_CH0_ERR, ENRES_CH1, ENRES_CH1_ERR) "
                         "VALUES (%i, '%s', %f, %f, %f, %f, %f, %f)",
                         table.Data(), address, (std::string(abs_path) + fout_name).Data(),
                         er_ave, er_ave_err, er_l, er_l_err, er_r, er_r_err);
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- energyres_krk writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbase, table, query, address);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- energyres_krk writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    return 0;
}
