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

void plotmuNchi2(const TString infile) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();
  
  //
  // Declare variables to read in ntuple
  //
  Int_t   isTightMuon1, isTightMuon2;
  Float_t muNchi21, muNchi22;
  
  //
  // Declare histograms
  //
  //
  TH1D *h1_muNchi21  = new TH1D("h1_muNchi21","",100,0,100);
        h1_muNchi21->SetStats(0);
        h1_muNchi21->SetLineColor(1);
        h1_muNchi21->GetXaxis()->SetTitle("muNchi2");
        h1_muNchi21->GetXaxis()->SetTitleSize(0.05);
        h1_muNchi21->GetXaxis()->SetTitleOffset(0.75);
        h1_muNchi21->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_muNchi21->GetYaxis()->SetTitleSize(0.055);
        h1_muNchi21->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_muNchi21_isTight = (TH1D*)h1_muNchi21->Clone("h1_muNchi21_isTight");
        h1_muNchi21_isTight->SetLineColor(2);
  
  TH1D *h1_muNchi22  = new TH1D("h1_muNchi22","",100,0,100);
        h1_muNchi22->SetStats(0);
        h1_muNchi22->SetLineColor(1);
        h1_muNchi22->GetXaxis()->SetTitle("muNchi2");
        h1_muNchi22->GetXaxis()->SetTitleSize(0.05);
        h1_muNchi22->GetXaxis()->SetTitleOffset(0.75);
        h1_muNchi22->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_muNchi22->GetYaxis()->SetTitleSize(0.055);
        h1_muNchi22->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_muNchi22_isTight = (TH1D*)h1_muNchi22->Clone("h1_muNchi22_isTight");
        h1_muNchi22_isTight->SetLineColor(2);

  //
  // Setup input ntuple
  //
  inputFile = new TFile(infile);
  inputTree = (TTree*)inputFile->Get("selectZmm/Events");
  
  inputTree->SetBranchAddress("isTightMuon1", &isTightMuon1);   
  inputTree->SetBranchAddress("isTightMuon2", &isTightMuon2);   
  inputTree->SetBranchAddress("muNchi21", &muNchi21);   
  inputTree->SetBranchAddress("muNchi22", &muNchi22);   
           
  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    //
    // Fill histograms
    //
    h1_muNchi21->Fill(muNchi21);
    h1_muNchi22->Fill(muNchi22);
    if(isTightMuon1 == true)
    {
      h1_muNchi21_isTight->Fill(muNchi21);
    }
    if(isTightMuon2 == true)
    {
      h1_muNchi22_isTight->Fill(muNchi22);
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
  leg->AddEntry(h1_muNchi21,"muNchi21","lp");
  leg->AddEntry(h1_muNchi21_isTight,"muNchi21_isTight","lp");

  TLegend *leg2 = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg2->SetBorderSize(1);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.05084746);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(h1_muNchi22,"muNchi22","lp");
  leg2->AddEntry(h1_muNchi22_isTight,"muNchi22_isTight","lp");

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
  h1_muNchi21->Draw();
  h1_muNchi21_isTight->Draw("same");
  leg->Draw("same");
  c->Print("muNchi21.png");

  TCanvas *c2 = new TCanvas("c2","c2",18,40,700,500);
  c2->Range(18.78366,0.5855855,56.2702,0.9308669);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  c2->cd();
  h1_muNchi22->Draw();
  h1_muNchi22_isTight->Draw("same");
  leg2->Draw("same");
  c2->Print("muNchi22.png");
}
