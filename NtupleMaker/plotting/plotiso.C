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

void plotiso(const TString infile) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();
  
  //
  // Declare variables to read in ntuple
  //
  Int_t   isTightMuon1, isTightMuon2;
  Float_t pfCombRelIso1, pfCombRelIso2;
  
  //
  // Declare histograms
  //
  //
  TH1D *h1_pfCombRelIso1  = new TH1D("h1_pfCombRelIso1","",100,0,10);
        h1_pfCombRelIso1->SetStats(0);
        h1_pfCombRelIso1->SetLineColor(1);
        h1_pfCombRelIso1->GetXaxis()->SetTitle("pfCombRelIso");
        h1_pfCombRelIso1->GetXaxis()->SetTitleSize(0.05);
        h1_pfCombRelIso1->GetXaxis()->SetTitleOffset(0.75);
        h1_pfCombRelIso1->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_pfCombRelIso1->GetYaxis()->SetTitleSize(0.055);
        h1_pfCombRelIso1->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_pfCombRelIso1_isTight = (TH1D*)h1_pfCombRelIso1->Clone("h1_pfCombRelIso1isTight");
        h1_pfCombRelIso1_isTight->SetLineColor(2);
  
  TH1D *h1_pfCombRelIso2  = new TH1D("h1_pfCombRelIso2","",100,0,10);
        h1_pfCombRelIso2->SetStats(0);
        h1_pfCombRelIso2->SetLineColor(1);
        h1_pfCombRelIso2->GetXaxis()->SetTitle("pfCombRelIso");
        h1_pfCombRelIso2->GetXaxis()->SetTitleSize(0.05);
        h1_pfCombRelIso2->GetXaxis()->SetTitleOffset(0.75);
        h1_pfCombRelIso2->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_pfCombRelIso2->GetYaxis()->SetTitleSize(0.055);
        h1_pfCombRelIso2->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_pfCombRelIso2_isTight = (TH1D*)h1_pfCombRelIso2->Clone("h1_pfCombRelIso2_isTight");
        h1_pfCombRelIso2_isTight->SetLineColor(2);

  //
  // Setup input ntuple
  //
  inputFile = new TFile(infile);
  inputTree = (TTree*)inputFile->Get("selectZmm/Events");
  
  inputTree->SetBranchAddress("isTightMuon1", &isTightMuon1);   
  inputTree->SetBranchAddress("isTightMuon2", &isTightMuon2);   
  inputTree->SetBranchAddress("pfCombRelIso1", &pfCombRelIso1);   
  inputTree->SetBranchAddress("pfCombRelIso2", &pfCombRelIso2);   
           
  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    //
    // Fill histograms
    //
    h1_pfCombRelIso1->Fill(pfCombRelIso1);
    h1_pfCombRelIso2->Fill(pfCombRelIso2);
    if(isTightMuon1 == true)
    {
      h1_pfCombRelIso1_isTight->Fill(pfCombRelIso1);
    }
    if(isTightMuon2 == true)
    {
      h1_pfCombRelIso2_isTight->Fill(pfCombRelIso2);
    }

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
  leg->AddEntry(h1_pfCombRelIso1,"pfCombRelIso1","lp");
  leg->AddEntry(h1_pfCombRelIso1_isTight,"pfCombRelIso1_isTight","lp");

  TLegend *leg2 = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg2->SetBorderSize(1);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.05084746);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(h1_pfCombRelIso2,"pfCombRelIso2","lp");
  leg2->AddEntry(h1_pfCombRelIso2_isTight,"pfCombRelIso2_isTight","lp");

  //
  // Save plots
  //
  TCanvas *c = new TCanvas("c", "c",18,40,700,500);
  c->Range(18.78366,0.5855855,56.2702,0.9308669);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetLogy();
  c->cd();
  h1_pfCombRelIso1->Draw();
  h1_pfCombRelIso1_isTight->Draw("same");
  leg->Draw("same");
  c->Print("pfCombRelIso1.png");

  TCanvas *c2 = new TCanvas("c2","c2",18,40,700,500);
  c2->Range(18.78366,0.5855855,56.2702,0.9308669);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  c2->cd();
  h1_pfCombRelIso2->Draw();
  h1_pfCombRelIso2_isTight->Draw("same");
  leg2->Draw("same");
  c2->Print("pfCombRelIso2.png");
}
