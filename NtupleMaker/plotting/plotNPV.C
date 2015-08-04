#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLatex.h>
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

void plotNPV(const TString infile) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();
  
  //
  // Declare variables to read in ntuple
  //
  Int_t npv;
  
  //
  // Declare histograms
  //
  //
  TH1D *hNPV  = new TH1D("hNPV","",100,0,100);
        hNPV->SetStats(0);
        hNPV->SetLineColor(1);
        hNPV->GetXaxis()->SetTitle("NPV");
        hNPV->GetXaxis()->SetTitleSize(0.05);
        hNPV->GetXaxis()->SetTitleOffset(0.75);
        hNPV->GetYaxis()->SetTitle("Events / 1.0 GeV");
        hNPV->GetYaxis()->SetTitleSize(0.055);
        hNPV->GetYaxis()->SetTitleOffset(0.9);


  //
  // Setup input ntuple
  //
  inputFile = new TFile(infile);
  inputTree = (TTree*)inputFile->Get("selectZmm/Events");
  
  inputTree->SetBranchAddress("npv", &npv);   // number of primary vertex
           
  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    //
    // Fill histograms
    //
    hNPV->Fill(npv);
  }

  TLegend *leg = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05084746);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(hNPV,"NPV","lp");


  //
  // Save plots
  //
  TCanvas *cmass = new TCanvas("cmass", "cmass",18,40,700,500);
  cmass->Range(18.78366,0.5855855,56.2702,0.9308669);
  cmass->SetFillColor(0);
  cmass->SetBorderMode(0);
  cmass->SetBorderSize(2);
  cmass->SetFrameBorderMode(0);
  cmass->cd();
  hNPV->Draw();
  leg->Draw("same");
  cmass->Print("npv.png");
}
