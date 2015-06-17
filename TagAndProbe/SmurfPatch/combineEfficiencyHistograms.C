
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"

#include <vector>
#include <iostream>

TH2D *new_histogram(std::vector<TString> fileNames, TString newName, TString hist);

void combineFakeRateHistograms()
{

    TString effDir = "/smurf/dlevans/FakeRates/V00-02-07_HCP_V0";
    TString name = "V00-02-07_HCP_V0";
    TFile outfile(effDir + "/summary.root", "RECREATE");
    TDirectory *outdir = (TDirectory*)outfile.GetDirectory("");

    gROOT->cd();
    std::vector<TString> fileNames;

    // muon selection
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet5_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold5_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold5_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold5_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet10_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold10_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold10_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold10_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet15_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold15_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold15_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold15_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet20_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold20_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold20_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold20_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet25_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold25_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold25_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold25_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet30_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold30_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold30_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold30_PtEta->SetDirectory(outdir);
    fileNames.push_back(effDir + "/HWW_Muon_"+name+"_FRJet50_IDISO/extra/smurf.root");
    TH2D *MuonFakeRate_M2_ptThreshold50_PtEta = new_histogram(fileNames, "MuonFakeRate_M2_ptThreshold50_PtEta", "Data");
    MuonFakeRate_M2_ptThreshold50_PtEta->SetDirectory(outdir);

    // electron selection

    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet15_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold15_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold15_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold15_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet20_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold20_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold20_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold20_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet25_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold25_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold25_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold25_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet30_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold30_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold30_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold30_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet35_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold35_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold35_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold35_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet40_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold40_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold40_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold40_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet45_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold45_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold45_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold45_PtEta->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+name+"_FRJet50_IDISO/extra/smurf.root");
    TH2D *ElectronFakeRate_V4_ptThreshold50_PtEta = new_histogram(fileNames, "ElectronFakeRate_V4_ptThreshold50_PtEta", "Data");
    ElectronFakeRate_V4_ptThreshold50_PtEta->SetDirectory(outdir);

    // write output
    outfile.Write();
    outfile.Close();

}


void combineEfficiencyHistograms()
{

    TString effDir = "/smurf/dlevans/Efficiencies/V00-02-07_trigNameFix_HCP_V1";
    TString nameEl = "V00-02-07_HCP_V0";
    TString nameMu = nameEl;
    TFile outfile(effDir + "/summary.root", "RECREATE");
    TDirectory *outdir = (TDirectory*)outfile.GetDirectory("");

    gROOT->cd();
    std::vector<TString> fileNames;

    // muon selection
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_ID/extra/smurf.root");
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_Iso/extra/smurf.root");
    TH2D *h2_results_muon_selection = new_histogram(fileNames, "h2_results_muon_selection", "SF");
    h2_results_muon_selection->SetDirectory(outdir);
   
    // electron selection
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_ID/extra/smurf.root");
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_Iso/extra/smurf.root");
    TH2D *h2_results_electron_selection = new_histogram(fileNames, "h2_results_electron_selection", "SF");
    h2_results_electron_selection->SetDirectory(outdir);

    // triggers

    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigSgl/extra/smurf.root");
    TH2D *h2_results_electron_single = new_histogram(fileNames, "h2_results_electron_single", "Data");
    h2_results_electron_single->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigLeadDbl/extra/smurf.root");
    TH2D *h2_results_electron_double_leadingleg = new_histogram(fileNames, "h2_results_electron_double_leadingleg", "Data");
    h2_results_electron_double_leadingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Electron_"+nameEl+"_NM1Eff_TrigTrailDbl/extra/smurf.root");
    TH2D *h2_results_electron_double_trailingleg = new_histogram(fileNames, "h2_results_electron_double_trailingleg", "Data");
    h2_results_electron_double_trailingleg->SetDirectory(outdir);

    std::cout << "changing namer for muons due to trigger matching fix..." << std::endl;
    nameMu = "V00-02-07_trigNameFix_HCP_V1";

    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigSgl/extra/smurf.root");
    TH2D *h2_results_muon_single = new_histogram(fileNames, "h2_results_muon_single", "Data");
    h2_results_muon_single->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigLeadDbl/extra/smurf.root");
    TH2D *h2_results_muon_double_leadingleg = new_histogram(fileNames, "h2_results_muon_double_leadingleg", "Data");
    h2_results_muon_double_leadingleg->SetDirectory(outdir);
    fileNames.clear();
    fileNames.push_back(effDir + "/HWW_Muon_"+nameMu+"_NM1Eff_TrigTrailDbl/extra/smurf.root");
    TH2D *h2_results_muon_double_trailingleg = new_histogram(fileNames, "h2_results_muon_double_trailingleg", "Data");
    h2_results_muon_double_trailingleg->SetDirectory(outdir);

    // write output
    outfile.Write();
    outfile.Close();

}

TH2D *new_histogram(std::vector<TString> fileNames, TString newName, TString hist)
{

    std::vector<TFile*> files;
    for (unsigned int i = 0; i < fileNames.size(); ++i)
        files.push_back(new TFile(fileNames[i], "READ"));

    gROOT->cd();
    files[0]->ls();

    TH2D *new_hist = 0;
    if (hist == "Data") 
        new_hist = (TH2D*)files[0]->Get("smurf_"+hist+"Eff")->Clone(newName);
    else 
        new_hist = (TH2D*)files[0]->Get("smurf_"+hist)->Clone(newName);
    
    if (new_hist == 0) {
        std::cout << "ERROR" << std::endl;
        return 0;
    }

    for (int ix = 1; ix <= new_hist->GetXaxis()->GetNbins() + 1; ++ix) {
        for (int iy = 1; iy <= new_hist->GetYaxis()->GetNbins() + 1; ++iy) {

            double sf       = 1.0;
            double err2h   = 0.0;
            double err2l    = 0.0;
            for (unsigned int f = 0; f < files.size(); ++f) {
                TH2D *temp_sf   = 0;
                if (hist == "Data") temp_sf = (TH2D*)files[f]->Get("smurf_"+hist+"Eff");
                else temp_sf = (TH2D*)files[f]->Get("smurf_"+hist);
                TH2D *temp_errh = (TH2D*)files[f]->Get("smurf_"+hist+"ErrHigh");
                TH2D *temp_errl = (TH2D*)files[f]->Get("smurf_"+hist+"ErrLow");
                sf      *= temp_sf->GetBinContent(ix, iy);
                err2h   += pow(temp_errh->GetBinContent(ix, iy)/temp_sf->GetBinContent(ix, iy), 2);
                err2l   += pow(temp_errl->GetBinContent(ix, iy)/temp_sf->GetBinContent(ix, iy), 2);
                delete temp_sf;
                delete temp_errl;
                delete temp_errh;
            }

            new_hist->SetBinContent(ix, iy, sf);
            new_hist->SetBinError(ix, iy, sf * (sqrt(err2h) + sqrt(err2l))/2.0);

        }
    }

    for (unsigned int i = 0; i < fileNames.size(); ++i) delete files[i];
    return new_hist;

}

