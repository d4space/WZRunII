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

void plotdB(const TString infile) {

  TFile* inputFile = new TFile();
  TTree* inputTree = new TTree();
  
  //
  // Declare variables to read in ntuple
  //
  Int_t   isTightMuon1, isTightMuon2;
  Float_t dB1, dB2;
  
  //
  // Declare histograms
  //
  //
  TH1D *h1_dB1  = new TH1D("h1_dB1","",100,0,1);
        h1_dB1->SetStats(0);
        h1_dB1->SetLineColor(1);
        h1_dB1->GetXaxis()->SetTitle("dB");
        h1_dB1->GetXaxis()->SetTitleSize(0.05);
        h1_dB1->GetXaxis()->SetTitleOffset(0.75);
        h1_dB1->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_dB1->GetYaxis()->SetTitleSize(0.055);
        h1_dB1->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_dB1_isTight = (TH1D*)h1_dB1->Clone("h1_dB1isTight");
        h1_dB1_isTight->SetLineColor(2);
  
  TH1D *h1_dB2  = new TH1D("h1_dB2","",100,0,1);
        h1_dB2->SetStats(0);
        h1_dB2->SetLineColor(1);
        h1_dB2->GetXaxis()->SetTitle("dB");
        h1_dB2->GetXaxis()->SetTitleSize(0.05);
        h1_dB2->GetXaxis()->SetTitleOffset(0.75);
        h1_dB2->GetYaxis()->SetTitle("Events / 1.0 GeV");
        h1_dB2->GetYaxis()->SetTitleSize(0.055);
        h1_dB2->GetYaxis()->SetTitleOffset(0.9);
  TH1D *h1_dB2_isTight = (TH1D*)h1_dB2->Clone("h1_dB2_isTight");
        h1_dB2_isTight->SetLineColor(2);

  //
  // Setup input ntuple
  //
  inputFile = new TFile(infile);
  inputTree = (TTree*)inputFile->Get("selectZmm/Events");
  
  inputTree->SetBranchAddress("isTightMuon1", &isTightMuon1);   
  inputTree->SetBranchAddress("isTightMuon2", &isTightMuon2);   
  inputTree->SetBranchAddress("dB1", &dB1);   
  inputTree->SetBranchAddress("dB2", &dB2);   
           
  for(int jentry=0;jentry<inputTree->GetEntries();jentry++) {
    inputTree->GetEntry(jentry);

    //
    // Fill histograms
    //
    h1_dB1->Fill(dB1);
    h1_dB2->Fill(dB2);
    if(isTightMuon1 == true)
    {
      h1_dB1_isTight->Fill(dB1);
    }
    if(isTightMuon2 == true)
    {
      h1_dB2_isTight->Fill(dB2);
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
  leg->AddEntry(h1_dB1,"dB1","lp");
  leg->AddEntry(h1_dB1_isTight,"dB1_isTight","lp");

  TLegend *leg2 = new TLegend(0.5373563,0.7097458,0.8577586,0.8474576,NULL,"brNDC");
  leg2->SetBorderSize(1);
  leg2->SetTextFont(62);
  leg2->SetTextSize(0.05084746);
  leg2->SetLineColor(0);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry(h1_dB2,"dB2","lp");
  leg2->AddEntry(h1_dB2_isTight,"dB2_isTight","lp");

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
  h1_dB1->Draw();
  h1_dB1_isTight->Draw("same");
  leg->Draw("same");
  c->Print("dB1.png");

  TCanvas *c2 = new TCanvas("c2","c2",18,40,700,500);
  c2->Range(18.78366,0.5855855,56.2702,0.9308669);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(2);
  c2->SetFrameBorderMode(0);
  c2->SetLogy();
  c2->cd();
  h1_dB2->Draw();
  h1_dB2_isTight->Draw("same");
  leg2->Draw("same");
  c2->Print("dB2.png");
}
