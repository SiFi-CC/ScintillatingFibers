// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            peakfin_pmi.cc             *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFInfo.hh"
#include "SFPeakFinder.hh"
#include "SFTools.hh"
#include "common_options_pmi.h"

#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>

#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char** argv)
{
    //----- parsing arguments
    
    TString path;
    TString outdir;
    TString dbase;
    int m = 0;
    int l = 0;
    int f = 0;
    
    if (argc < 2)
    {
        std::cerr << "to run type: ./peakfin_pmi path/to/data/ ";
        std::cerr << "-m 0 -l 0 -f 0 -out path/to/results/ ";
        std::cerr << "-db path/to/database.db" << std::endl;
        return 1;
    }
    
    int ret = parse_common_options(argc, argv, path, outdir, dbase, m, l, f);
    if (ret != 0)
        exit(ret);
    
    TString logname = path + "regular_measurement_list.txt";
    
    //----- setting measurement properties
    
    std::ifstream logfile(logname);
    
    if(!logfile.is_open())
    {
        std::cerr << "##### Error in peakfin_pmi!" << std::endl;
        std::cerr << "Couldn't open logfile: " << logname << std::endl;
        return 1;
    }
    
    std::string line;
    getline(logfile, line);
    
    int npoints = 0;
    std::vector<TString> meas_names;
    std::vector<double> positions;
    std::vector<int> times;
    
    TString name_tmp;
    double pos_tmp = 0;
    double time_tmp = 0;
    
    while (true)
    {
        logfile >> name_tmp >> pos_tmp >> time_tmp;
        
        if (name_tmp == "" || name_tmp == " ")
            break;
        
        std::cout << name_tmp << "\t" << pos_tmp << "\t" << time_tmp << std::endl;
        meas_names.push_back(name_tmp);
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
    info->SetCollimator("FBC");
    info->SetTestBench("PMI");
    info->SetLogFile(logname);
    info->SetNames(meas_names);
    info->SetPositions(positions);
    info->SetTimes(times);
    
    //----- getting histograms
    
    std::vector<TH1D*> hQL;
    std::vector<TH1D*> hQR;
    std::vector<TH1D*> hQAve;

    for (int i=0; i<npoints; i++)
    {
        TString full_root_path = path + meas_names[i] + "/sifi_results_HISTOS.root";
        TString dir_name = Form("/Module%i/Layer%i/Fiber%i/", m, l, f);
        std::cout << "\nOpening file: " << full_root_path << std::endl;
        
        auto file = std::make_unique<TFile>(full_root_path, "READ");
        
        if (!file.get()->IsOpen())
        {
            std::cerr << "##### Error in peakfin_pmi! Could not open file: " << std::endl;
            std::cerr << "full_root_path" << std::endl;
            std::abort();
        }
        
        //----- spectrum - left
        TString hname = Form("Q_M%iL%iF%iL", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);
        auto htemp_L = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_L)
        {
            std::cerr << "##### Error in peakfin_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_L = (TH1D*)htemp_L->Clone(hname);
        
        if (!clone_L)
        {
            std::cerr << "##### Error in peakfin_pmi! Cloning went wrong: " << std::endl;
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
            std::cerr << "##### Error in peakfin_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_R = (TH1D*)htemp_R->Clone(hname);
        
        if (!clone_R)
        {
            std::cerr << "##### Error in peakfin_pmi! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_R->SetDirectory(nullptr);
        hQR.push_back(clone_R);
        
        //----- spectrum - average
        hname = Form("QAve_M%iL%iF%i", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        file->cd(dir_name);

        auto htemp_A = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp_A) 
        {
            std::cerr << "##### Error in peakfin_pmi! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto clone_A = (TH1D*)htemp_A->Clone(hname);
        
        if (!clone_A)
        {
            std::cerr << "##### Error in peakfin_pmi! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        clone_A->SetDirectory(nullptr);
        hQAve.push_back(clone_A);
    }
    
    //----- fitting
    std::vector<SFResults*> params_L;
    std::vector<SFResults*> params_R;
    std::vector<SFResults*> params_Ave;
    
    for (int i=0; i<npoints; i++)
    {
        TString path_meas = path + meas_names[i];
        
        params_L.push_back(SFPeakFinder::FindPeakFit(hQL[i], path_meas, 1, 1));
        
        params_R.push_back(SFPeakFinder::FindPeakFit(hQR[i], path_meas, 1, 1));
        
        params_Ave.push_back(SFPeakFinder::FindPeakFit(hQAve[i], path_meas, 1, 1));
    }
    
    TCanvas *can_QL = new TCanvas("can_QL", "can_QL", 2000, 1200);
    TCanvas *can_QR = new TCanvas("can_QR", "can_QR", 2000, 1200);
    TCanvas *can_QAve = new TCanvas("can_QAve", "can_QAve", 2000, 1200);
    
    can_QL->DivideSquare(npoints);
    can_QR->DivideSquare(npoints);
    can_QAve->DivideSquare(npoints);
    
    TLatex text;
    text.SetNDC(true);
    
    double ymax = 0;
    
    for (int i=0; i<npoints; i++)
    {
        
        can_QL->cd(i+1);
        gPad->SetGrid(1, 1);
        ymax = hQL[i]->GetBinContent(hQL[i]->GetMaximumBin());
        hQL[i]->GetYaxis()->SetRangeUser(0, ymax + 0.25*ymax);
        hQL[i]->Draw();
        hQL[i]->SetStats(0);
        hQL[i]->SetTitle(Form("Charge spectrum (left) %.2f mm", positions[i]));
        
        text.DrawLatex(0.50, 0.75, Form("c = %.3f +/- %.3f",
                       params_L[i]->GetValue(SFResultTypeNum::kPeakConst),
                       params_L[i]->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.50, 0.70, Form("#mu = %.3f +/- %.3f",
                       params_L[i]->GetValue(SFResultTypeNum::kPeakPosition),
                       params_L[i]->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.50, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       params_L[i]->GetValue(SFResultTypeNum::kPeakSigma),
                       params_L[i]->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.50, 0.60, Form("#chi^{2}/NDF = %.3f",
                       params_L[i]->GetValue(SFResultTypeNum::kChi2NDF)));
        
        can_QR->cd(i+1);
        gPad->SetGrid(1, 1);
        ymax = hQR[i]->GetBinContent(hQR[i]->GetMaximumBin());
        hQR[i]->GetYaxis()->SetRangeUser(0, ymax + 0.25*ymax);
        hQR[i]->Draw();
        hQR[i]->SetStats(0);
        hQR[i]->SetTitle(Form("Charge spectrum (right) %.2f mm", positions[i]));
        
        text.DrawLatex(0.50, 0.75, Form("c = %.3f +/- %.3f",
                       params_R[i]->GetValue(SFResultTypeNum::kPeakConst),
                       params_R[i]->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.50, 0.70, Form("#mu = %.3f +/- %.3f",
                       params_R[i]->GetValue(SFResultTypeNum::kPeakPosition),
                       params_R[i]->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.50, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       params_R[i]->GetValue(SFResultTypeNum::kPeakSigma),
                       params_R[i]->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.50, 0.60, Form("#chi^{2}/NDF = %.3f",
                       params_R[i]->GetValue(SFResultTypeNum::kChi2NDF)));
        
        can_QAve->cd(i+1);
        gPad->SetGrid(1, 1);
        ymax = hQAve[i]->GetBinContent(hQAve[i]->GetMaximumBin());
        hQAve[i]->GetYaxis()->SetRangeUser(0, ymax + 0.25*ymax);
        hQAve[i]->Draw();
        hQAve[i]->SetStats(0);
        hQAve[i]->SetTitle(Form("Average charge spectrum %.2f mm", positions[i]));
        
        text.DrawLatex(0.50, 0.75, Form("c = %.3f +/- %.3f",
                       params_Ave[i]->GetValue(SFResultTypeNum::kPeakConst),
                       params_Ave[i]->GetUncertainty(SFResultTypeNum::kPeakConst)));
        text.DrawLatex(0.50, 0.70, Form("#mu = %.3f +/- %.3f",
                       params_Ave[i]->GetValue(SFResultTypeNum::kPeakPosition),
                       params_Ave[i]->GetUncertainty(SFResultTypeNum::kPeakPosition)));
        text.DrawLatex(0.50, 0.65, Form("#sigma = %.3f +/- %.3f", 
                       params_Ave[i]->GetValue(SFResultTypeNum::kPeakSigma),
                       params_Ave[i]->GetUncertainty(SFResultTypeNum::kPeakSigma)));
        text.DrawLatex(0.50, 0.60, Form("#chi^{2}/NDF = %.3f",
                       params_Ave[i]->GetValue(SFResultTypeNum::kChi2NDF)));
    }
    
    //---- saving in ROOT file
    TString fout_name = Form("peakfin_M%iL%iF%i.root", m, l, f);
    TFile *out_file = new TFile(outdir + fout_name, "RECREATE");
    
    if (!out_file->IsOpen())
    {
        std::cerr << "##### Error in peakfin_pmi.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << outdir + fout_name << std::endl;
        return 1;
    }
    
    can_QL->Write();
    can_QR->Write();
    can_QAve->Write();
    
    out_file->mkdir("histograms");
    out_file->cd("histograms");
    
    for (int i=0; i<npoints; i++)
    {
        hQL[i]->Write();
        hQR[i]->Write();
        hQAve[i]->Write();
    }
    
    out_file->Close();
    
    //----- saving in data base
    
    TString address_str = Form("%i%i%i", m, l, f);
    int address = atoi(address_str);
    
    TString table = "PEAK_FINDER";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE) VALUES(%i, '%s')",
                         table.Data(), address, (outdir + fout_name).Data());
    
    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- peakfin_pmi writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbase, table, query, address);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- peakfin_pmi writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    return 0;
}



