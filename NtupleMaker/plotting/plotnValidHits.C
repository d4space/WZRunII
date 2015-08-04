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

void plotnValidHits(const TString infile) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();
  
  //
  // Declare variables to read in ntuple
  //
  Int_t   isTightMuon1, isTightMuon2;
  UInt_t  nValidHits1, nValidHits2;
  
  //
  // Declare histograms
  //
  //
  TH1D *h1_nValidHits1  = new TH1D("h1_nValidHits1","",10,0,10);
        h1_nValidHits1->SetStats(0);
        h1_nValidHits1->SetLineColor(1);
        h1_nValidHits1->GetXaxis()->SetTitle("nValidHits");
        h1_nValidHits1->GetXaxis()->SetTitleSize(0.05);
        h1_nValidHits1->GetXaxis()->SetTitleOffset(0.75);
        h1_nValidHits1->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_nValidHits1->GetYaxis()->SetTitleSize(0.055);
        h1_nValidHits1->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_nValidHits1_isTight = (TH1D*)h1_nValidHits1->Clone("h1_nValidHits1_isTight");
        h1_nValidHits1_isTight->SetLineColor(2);
  
  TH1D *h1_nValidHits2  = new TH1D("h1_nValidHits2","",10,0,10);
        h1_nValidHits2->SetStats(0);
        h1_nValidHits2->SetLineColor(1);
        h1_nValidHits2->GetXaxis()->SetTitle("nValidHits");
        h1_nValidHits2->GetXaxis()->SetTitleSize(0.05);
        h1_nValidHits2->GetXaxis()->SetTitleOffset(0.75);
        h1_nValidHits2->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_nValidHits2->GetYaxis()->SetTitleSize(0.055);
        h1_nValidHits2->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_nValidHits2_isTight = (TH1D*)h1_nValidHits2->Clone("h1_nValidHits2_isTight");
        h1_nValidHits2_isTight->SetLineColor(2);

  //
  // Setup input ntuple
  //
  inputFile = new TFile(infile);
  inputTree = (TTree*)inputFile->Get("selectZmm/Events");
  
  inputTree->SetBranchAddress("isTightMuon1", &isTightMuon1);   
  inputTree->SetBranchAddress("isTightMuon2", &isTightMuon2);   
  inputTree->SetBranchAddress("nValidHits1", &nValidHits1);   
  inputTree->SetBranchAddress("nValidHits2", &nValidHits2);   
           
  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    //
    // Fill histograms
    //
    h1_nValidHits1->Fill(nValidHits1);
    h1_nValidHits2->Fill(nValidHits2);
    if(isTightMuon1 == true)
    {
      h1_nValidHits1_isTight->Fill(nValidHits1);
    }
    if(isTightMuon2 == true)
    {
      h1_nValidHits2_isTight->Fill(nValidHits2);
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
  leg->AddEntry(h1_nValidHits1,"nValidHits1","lp");
  leg->AddEntry(h1_nValidHits1_isTight,"nValidHits1_isTight","lp");

  TLegend *leg2 = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg2->SetBorderSize(1);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.05084746);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(h1_nValidHits2,"nValidHits2","lp");
  leg2->AddEntry(h1_nValidHits2_isTight,"nValidHits2_isTight","lp");

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
  h1_nValidHits1->Draw();
  h1_nValidHits1_isTight->Draw("same");
  leg->Draw("same");
  c->Print("nValidHits1.png");

  TCanvas *c2 = new TCanvas("c2","c2",18,40,700,500);
  c2->Range(18.78366,0.5855855,56.2702,0.9308669);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  c2->cd();
  h1_nValidHits2->Draw();
  h1_nValidHits2_isTight->Draw("same");
  leg2->Draw("same");
  c2->Print("nValidHits2.png");
}
