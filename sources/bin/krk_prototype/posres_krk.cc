// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *            posres_krk.cc              *
// *          Katarzyna Rusiecka           *
// * katarzyna.rusiecka@doctoral.uj.edu.pl *
// *          Created in 2021              *
// *                                       *
// *****************************************

#include "SFInfo.hh"
#include "SFPeakFinder.hh"
#include "SFTools.hh"
#include "SFAttenuation.hh"
#include "SFPositionReco.hh"
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
        std::cerr << "to run type: ./posres_krk path/to/data/log.txt ";
        std::cerr << "-m 0 -l 0 -f 0 -out -out path/to/results/ ";
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
        std::cerr << "##### Error in posres_krk!" << std::endl;
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
    info->SetTestBench("Prototype KRK");
    info->SetSiPM("Ketek");
    info->SetLogFile(log_path);
    info->SetNames(meas_names);
    info->SetPositions(positions);
    info->SetTimes(times);
    info->SetStartTimes(start);
    info->SetStopTimes(stop);
    
    double pos_uncert = 2.; // mm // TODO check
    double BL_sigma_cut = 5.;     // TODO check
    
    SFChAddr addr;
    addr.fModule = m;
    addr.fLayer  = l;
    addr.fFiber  = f;
    
    //----- getting histograms
    
    std::vector<TH1D*> hAvg;
    
    for (int i=0; i<npoints; i++)
    {
        //----- spectrum - average
        
        TString full_root_path = meas_names[i] + "/sifi_results_HISTOS.root";
        TString dir_name = Form("/Module%i/Layer%i/Fiber%i/", m, l, f);
        std::cout << "\nOpening file: " << full_root_path << std::endl;
        
        auto hfile = std::make_unique<TFile>(full_root_path, "READ");
        
        if (!hfile.get()->IsOpen())
        {
            std::cerr << "##### Error in posres_krk! Could not open file: " << std::endl;
            std::cerr << full_root_path << std::endl;
            std::abort();
        }
        
        TString hname = Form("QAve_M%iL%iF%i", m, l, f);
        std::cout << "\tGetting histogram: " << dir_name + hname << std::endl;
        hfile->cd(dir_name);
        auto htemp = (TH1D*)gDirectory->GetObjectChecked(hname, "TH1D");
        
        if (!htemp)
        {
            std::cerr << "##### Error in posres_krk! Could not get histogram: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }

        auto hclone = (TH1D*)htemp->Clone(hname);
        
        if (!hclone)
        {
            std::cerr << "##### Error in posres_krk! Cloning went wrong: " << std::endl;
            std::cerr << dir_name + hname << std::endl;
            std::abort();
        }
        
        hclone->SetDirectory(nullptr);
        hAvg.push_back(hclone);
    }
    
    //----- getting MLR curve (reversed)
    TString path_to_attenuation = std::string(abs_path) + 
                                  Form("/attenuation_M%iL%iF%i.root", m, l, f);
                                  
    std::cout << "\n\nabsolute path: " << abs_path << std::endl;
    std::cout << "path to attenuation: " << path_to_attenuation << std::endl;
                                  
    TFile *file_attenuation = new TFile(path_to_attenuation, "READ");
    
    if(!file_attenuation->IsOpen())
    {
        std::cerr << "Error in posres_krk! Could not open attenuation file!" << std::endl;
        std::cerr << path_to_attenuation << std::endl;
        return 1;
    }
    
    TCanvas *can_tmp = (TCanvas*)file_attenuation->FindObjectAny("can_mlr");
    TPad *pad_tmp = (TPad*)can_tmp->GetListOfPrimitives()->FindObject("can_mlr_1");
    
    TGraphErrors *gMLR = (TGraphErrors*)pad_tmp->FindObject("attenuation_mlr");
    
    if (!gMLR)
    {
        std::cerr << "##### Error in posres_krk! Could not find MLR attenuation graph!" << std::endl;
        return 1;
    }
    
    double *MLR_x = gMLR->GetX();
    double *MLR_y = gMLR->GetY();
    double *MLR_ex = gMLR->GetEX();
    double *MLR_ey = gMLR->GetEY();
    
    TGraphErrors *gMLR_rev = new TGraphErrors(npoints, MLR_y, MLR_x, MLR_ey, MLR_ex);
    gMLR_rev->SetTitle(Form("source position vs. M_{LR} (M%iL%iF%i)", m, l, f));
    gMLR_rev->GetXaxis()->SetTitle("M_{LR}");
    gMLR_rev->GetYaxis()->SetTitle("source position [mm]");
    
    TF1 *fMLR_rev = new TF1("fMLR_rev", "pol1", -2, 2);
    gMLR_rev->Fit(fMLR_rev, "QR");
    
    file_attenuation->Close();
    
    //----- position-wise position resolution
    auto tuple_results = SFPositionReco::ReconstructPositionDistAll(addr,
                                                                    hAvg,
                                                                    positions,
                                                                    meas_names,
                                                                    fMLR_rev,
                                                                    BL_sigma_cut,
                                                                    pos_uncert,
                                                                    info->GetCollimator());
    
    SFResults *results_single = std::get<0>(tuple_results);
    std::vector<TH1D*> hPosRecoAll = std::get<1>(tuple_results);
    
    TGraphErrors *gPosRecoVsPos = (TGraphErrors*)results_single->GetObject(SFResultTypeObj::kPosRecoVsPosGraph);
    TGraphErrors *gPosResVsPos = (TGraphErrors*)results_single->GetObject(SFResultTypeObj::kPosResVsPosGraph);
    TGraphErrors *gPosRecoDiff = (TGraphErrors*)results_single->GetObject(SFResultTypeObj::kPosDiffVsPosGraph);
    TGraphErrors* gResiduals = (TGraphErrors*)results_single->GetObject(SFResultTypeObj::kResidualGraph);
    
    //----- summed position resolution
    
    SFResults *results_summed = SFPositionReco::ReconstructPositionDistSum(addr,
                                                                           hAvg,
                                                                           positions,
                                                                           meas_names,
                                                                           fMLR_rev,
                                                                           BL_sigma_cut,
                                                                           info->GetCollimator());
    
    TH1D* hPosRecoSum = (TH1D*)results_summed->GetObject(SFResultTypeObj::kPositionDist);
    
    //----- drawing
    
    double sigma_to_fwhm = 2 * sqrt(2 * log(2));
    
    TLatex text;
    text.SetNDC(true);
    text.SetTextSize(0.03);
    
    //----- distribution of reconstructed positions
    
    TCanvas *can_pos_reco_all = new TCanvas("can_pos_reco_all", "can_pos_reco_all", 2000, 1200);
    can_pos_reco_all->DivideSquare(npoints + 1);
    
    for (int i=0; i<npoints; i++)
    {
        can_pos_reco_all->cd(i + 1);
        gPad->SetGrid(1, 1);
        hPosRecoAll[i]->Draw();
        TF1* fun_gaus_tmp = hPosRecoAll[i]->GetFunction("fun_gaus");
        double fwhm = fun_gaus_tmp->GetParameter(2) * sigma_to_fwhm;
        double fwhm_err = fun_gaus_tmp->GetParError(2) * sigma_to_fwhm;
        double pos = fun_gaus_tmp->GetParameter(1);
        double pos_err = fun_gaus_tmp->GetParError(1);
        text.DrawLatex(0.20, 0.75, Form("#mu = (%.2f +/- %.2f) mm", pos, pos_err));
        text.DrawLatex(0.20, 0.70, Form("FWHM = (%.2f =/- %.2f) mm", fwhm, fwhm_err));
    }
    
    can_pos_reco_all->cd(npoints + 1);
    gPad->SetGrid(1, 1);
    gPad->SetFillColor(kGreen - 10);
    hPosRecoSum->SetLineColor(kBlack);
    hPosRecoSum->Draw();
    
    TF1 *fun_gaus_tmp = hPosRecoSum->GetFunction("fun_gaus");
    double fwhm_sum = fun_gaus_tmp->GetParameter(2) * sigma_to_fwhm;
    double fwhm_sum_err = fun_gaus_tmp->GetParError(2) * sigma_to_fwhm;
    double pos_sum = fun_gaus_tmp->GetParameter(1);
    double pos_sum_err = fun_gaus_tmp->GetParError(1);
    
    text.DrawLatex(0.20, 0.75, Form("#mu = (%.2f +/- %.2f) mm", pos_sum, pos_sum_err));
    text.DrawLatex(0.20, 0.70, Form("FWHM = (%.2f +/- %.2f) mm", fwhm_sum, fwhm_sum_err));
    
    //----- average spectra
    
    TCanvas *can_spec_avg = new TCanvas("can_spec_avg", "can_spec_avg", 2000, 1200);
    can_spec_avg->DivideSquare(npoints);
    
    TLine rline;
    rline.SetLineColor(kRed);
    rline.SetLineWidth(1);
    
    for (int i=0; i<npoints; i++)
    {
        can_spec_avg->cd(i + 1);
        gPad->SetGrid(1, 1);
        TString new_title = TString(hAvg[i]->GetTitle()) + 
                            Form(" position %.1f mm", positions[i]);
        hAvg[i]->SetTitle(new_title);
        hAvg[i]->Draw();
        hAvg[i]->SetStats(false);
        
        TF1* fun_tmp = hAvg[i]->GetFunction(Form("f_QAve_M%iL%iF%i", m, l, f));
        
        if(!fun_tmp)
        {
            std::cerr << "##### Error in posres_krk! Could not get spectrum function!" << std::endl;
            return 1;
        }
        
        double max_value = hAvg[i]->GetBinContent(hAvg[i]->GetMaximumBin());
        double mean = fun_tmp->GetParameter(1);
        double sigma = fabs(fun_tmp->GetParameter(2));
        
        rline.DrawLine(mean - 2*sigma, 0, mean - 2*sigma, max_value);
        rline.DrawLine(mean + 2*sigma, 0, mean + 2*sigma, max_value);
    }
    
    //----- reversed MLR
    
    TCanvas *can_mlr_rev = new TCanvas("can_mlr_rev", "can_rev_mlr", 1200, 1000);
    gPad->SetGrid(1, 1);
    gMLR_rev->Draw("AP");
    
    text.DrawLatex(0.6, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.6, 0.75, Form("p_{0} = %.2f +/- %.2f",
                   fMLR_rev->GetParameter(0), fMLR_rev->GetParError(0)));
    text.DrawLatex(0.6, 0.70, Form("p_{1} = %.2f +/- %.2f",
                   fMLR_rev->GetParameter(1), fMLR_rev->GetParError(1)));
    
    //----- position reco graph
    
    TCanvas *can_posreco = new TCanvas("can_posreco", "can_posreco", 1200, 1000);
    TPad *pad_posreco = new TPad("pad_posreco", "pad_posreco", 0, 0.3, 1, 1, 10, 0);
    TPad *pad_res = new TPad("pad_res", "pad_res", 0, 0, 1, 0.3, 10, 0);
    
    pad_posreco->Draw();
    pad_posreco->cd();
    gPad->SetGrid(1, 1);
    gPosRecoVsPos->Draw("AP");
    
    TF1 *funpol1 = gPosRecoVsPos->GetFunction("fun_pol1");
    
    text.DrawLatex(0.2, 0.80, "y = p_{0} + p_{1}x");
    text.DrawLatex(0.2, 0.75, Form("p_{0} = %.2f +/- %.2f", 
                   funpol1->GetParameter(0), funpol1->GetParError(0)));
    text.DrawLatex(0.2, 0.70, Form("p_{1} = %.2f +/- %.2f", 
                   funpol1->GetParameter(1), funpol1->GetParError(1)));
    text.DrawLatex(0.2, 0.65, Form("#chi^{2}/NDF = %.2f", 
                   funpol1->GetChisquare() / funpol1->GetNDF()));
    
    can_posreco->cd();
    pad_res->Draw();
    pad_res->cd();
    gPad->SetGrid(1, 1);
    gResiduals->Draw("AP");
    
    //----- position resolution graph
    
    TCanvas *can_posres = new TCanvas("can_posres", "can_posres", 1200, 1000);
    gPad->SetGrid(1, 1);
    
    TGraphErrors *gPosResSum = new TGraphErrors(1);
    gPosResSum->SetMarkerStyle(5);
    gPosResSum->SetPoint(0, 100, fwhm_sum);
    gPosResSum->SetPointError(0, pos_uncert, fwhm_sum_err);
    
    TH1D *htemp = new TH1D("htemp", gPosResVsPos->GetTitle(), 105, 0, 105);
    htemp->GetXaxis()->SetTitle(gPosResVsPos->GetXaxis()->GetTitle());
    htemp->GetYaxis()->SetTitle(gPosResVsPos->GetYaxis()->GetTitle());
    htemp->SetStats(0);
    htemp->Draw();
    
    gPosResVsPos->Draw("P");
    gPosResSum->Draw("P");
    
    double  ymin = TMath::Min(TMath::MinElement(npoints, gPosResVsPos->GetY()), gPosResSum->GetPointY(0));
    double  ymax = TMath::Max(TMath::MaxElement(npoints, gPosResVsPos->GetY()), gPosResSum->GetPointY(0));
    
    htemp->GetYaxis()->SetRangeUser(ymin - 0.05 * fabs(ymin), ymax + 0.05 * fabs(ymax));
    
    text.DrawLatex(0.2, 0.8, Form("FWHM_{avg} = (%.2f +/- %.2f) mm",
                   results_single->GetValue(SFResultTypeNum::kPositionRes),
                   results_single->GetUncertainty(SFResultTypeNum::kPositionRes)));
    
    text.SetNDC(false);
    text.DrawLatex(70, fwhm_sum + 0.2, Form("FWHM_{sum} = (%.2f +/- %.2f) mm", fwhm_sum, fwhm_sum_err));
    
    //----- position reconstruction difference
    
    TCanvas *can_diff = new TCanvas("can_diff", "can_diff", 1200, 1000);
    gPad->SetGrid(1, 1);
    gPosRecoDiff->Draw("AP");
    
    //----- saving to ROOT file

    TString fout_name = Form("/posres_M%iL%iF%i.root", m, l, f);
    TFile *out_file = new TFile(std::string(abs_path) + fout_name, "RECREATE");
    
    if (!out_file->IsOpen())
    {
        std::cerr << "##### Error in posres_krk.cc!" << std::endl;
        std::cerr << "Couldn't open file: " << std::string(abs_path) + fout_name << std::endl;
        return 1;
    }
    
    can_pos_reco_all->Write();
    can_spec_avg->Write();
    can_posreco->Write();
    can_posres->Write();
    can_diff->Write();
    can_mlr_rev->Write();
    
    out_file->mkdir("histograms");
    out_file->cd("histograms");
    
    for (int i=0; i<npoints; i++)
    {
        hAvg[i]->Write();
        hPosRecoAll[i]->Write();
    }
    
    hPosRecoSum->Write();
    
    out_file->Close();
    
    //----- saving to the data base
    
    TString address_str = Form("%i%i%i", m, l, f);
    int address = atoi(address_str);
    
    TString table = "POSITION_RESOLUTION";
    TString query = Form("INSERT OR REPLACE INTO %s (SERIES_ID, RESULTS_FILE, POSITION_RES_POL3, "
                         "POSITION_RES_POL3_ERR, POSITION_RES_POL3_ALL, POSITION_RES_POL3_ALL_ERR, "
                         "POSITION_RES_POL1, POSITION_RES_POL1_ERR, POSITION_RES_POL1_ALL, "
                         "POSITION_RES_POL1_ALL_ERR) VALUES(%i, '%s', %f, %f, %f, %f, %f, %f, %f, %f)", 
                         table.Data(), address, (std::string(abs_path) + fout_name).Data(), 
                         -1.0, -1.0, -1.0, -1.0,
                         results_single->GetValue(SFResultTypeNum::kPositionRes),
                         results_single->GetUncertainty(SFResultTypeNum::kPositionRes),
                         fwhm_sum, fwhm_sum_err);

    const int max_tries = 20;
    int       i_try     = max_tries;
    float     wait      = 0;
    bool      stat      = false;

    srand(time(NULL));

    do
    {
        std::cout << "----- posres_krk writing, try number " << (max_tries - i_try) + 1 << std::endl;
        stat = SFTools::SaveResultsDB(dbase, table, query, address);

        if (stat) break;

        --i_try;
        wait = rand() % 10 + 1;
        std::cout << "----- posres_krk writing, database locked..." << std::endl;
        std::cout << "----- waiting " << wait << " s" << std::endl;
        sleep(wait);
    } while (i_try > 0);
    
    return 0;
}
    
