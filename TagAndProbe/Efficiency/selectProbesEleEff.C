//================================================================================================
//
// Select probes for electron efficiencies with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <vector>                         // STL vector class
#include <iostream>                       // standard I/O
#include <iomanip>                        // functions to format standard I/O
#include "Math/LorentzVector.h"           // 4-vector class

// structure for output ntuple
#include "EffData.hh" 
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void selectProbesEleEff(const TString infilename,           // input ntuple
                        const TString outputDir,            // output directory
			const Int_t   effType,              // type of efficiency to compute
		        const Bool_t  doGenMatch = kFALSE,  // match to generator leptons
			const Bool_t  doWeighted = kFALSE   // store events with weights
) {
  gBenchmark->Start("selectProbesEleEff");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
 
  const Double_t TAG_PT_CUT = 25;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  enum { eIDISO, eTrigger };  // event category enum
  if(effType > eTrigger) {
    cout << "Invalid effType option! Exiting..." << endl;
    return;
  }
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:pass:evtNum");

  //
  // Declare output ntuple variables
  //
  Int_t  matchGen;
  Int_t  npv;
  Float_t  genVPt, genVPhi, genVy, genVMass;
  Float_t  rawpfmet, rawpfmetPhi;
  Float_t  type1pfmet, type1pfmetPhi;
  Float_t genmet, genmetPhi, u1, u2;
  Int_t   q1, q2;
  Float_t   pfChIso1, pfChIso2;
  LorentzVector *dilep=0, *lep1=0, *lep2=0;
  LorentzVector *sc1=0, *sc2=0;
  Int_t isVetoEle1, isLooseEle1, isMediumEle1, isTightEle1;
  Int_t isVetoEle2, isLooseEle2, isMediumEle2, isTightEle2;
  Int_t passSingleEleTrigger;
  Int_t matchTrigObj1, matchTrigObj2;
  
  // Read input file and get the TTrees
  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename);	 assert(infile);
  TTree *intree = (TTree*)infile->Get("Events"); assert(intree);
  
  intree->SetBranchAddress("matchGen",   &matchGen);    // event has both leptons matched to MC Z->ll
  intree->SetBranchAddress("npv",        &npv);	        // number of primary vertices
  intree->SetBranchAddress("genVPt",     &genVPt);	        
  intree->SetBranchAddress("genVPhi",    &genVPhi);	        
  intree->SetBranchAddress("genVy",    	 &genVy);	        
  intree->SetBranchAddress("genVMass",   &genVMass);	        
  intree->SetBranchAddress("rawpfmet",   &rawpfmet);	        
  intree->SetBranchAddress("rawpfmetPhi",   &rawpfmetPhi);	        
  intree->SetBranchAddress("type1pfmet",   &type1pfmet);	        
  intree->SetBranchAddress("type1pfmetPhi",   &type1pfmetPhi);	        
  intree->SetBranchAddress("genmet",        &genmet);	        // MET
  intree->SetBranchAddress("genmetPhi",     &genmetPhi);      // phi(MET)
  intree->SetBranchAddress("u1",         &u1);	        // parallel component of recoil
  intree->SetBranchAddress("u2",         &u2);	        // perpendicular component of recoil
  intree->SetBranchAddress("q1",         &q1);	        // charge of tag lepton
  intree->SetBranchAddress("q2",         &q2);	        // charge of probe lepton
  intree->SetBranchAddress("pfChIso1",         &pfChIso1);
  intree->SetBranchAddress("pfChIso2",         &pfChIso2);
  intree->SetBranchAddress("dilep",      &dilep);       // dilepton 4-vector
  intree->SetBranchAddress("lep1",       &lep1);        // tag lepton 4-vector
  intree->SetBranchAddress("lep2",       &lep2);        // probe lepton 4-vector
  intree->SetBranchAddress("sc1",       &sc1);        // sc probe lepton 4-vector
  intree->SetBranchAddress("sc2",       &sc2);        // sc probe lepton 4-vector
  intree->SetBranchAddress("isVetoEle1",       &isVetoEle1);
  intree->SetBranchAddress("isLooseEle1",       &isLooseEle1);
  intree->SetBranchAddress("isMediumEle1",       &isMediumEle1);
  intree->SetBranchAddress("isTightEle1",       &isTightEle1);
  intree->SetBranchAddress("isVetoEle2",       &isVetoEle2);
  intree->SetBranchAddress("isLooseEle2",       &isLooseEle2);
  intree->SetBranchAddress("isMediumEle2",       &isMediumEle2);
  intree->SetBranchAddress("isTightEle2",       &isTightEle2);
  intree->SetBranchAddress("passSingleEleTrigger",       &passSingleEleTrigger);
  intree->SetBranchAddress("matchTrigObj1",       &matchTrigObj1);
  intree->SetBranchAddress("matchTrigObj2",       &matchTrigObj2);

  //
  // loop over events
  //
  double nProbeTotal = 0;
  double nProbePass = 0;
  for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
   
    if (sc1->Pt() < TAG_PT_CUT) continue;
    
    // check GEN match if necessary
    //if(doGenMatch && !matchGen) continue;
    
    Bool_t  pass=kFALSE;
    Float_t mass=0;
    if(effType==eIDISO)
    {
      // probe = slimmed electron
      // pass = MidiumID and Charged Isolation Cut < 0.15 * probe pT
      if(isMediumEle2 && pfChIso2 < 15*sc2->Pt()) {pass=kTRUE;}
      else {pass=kFALSE;}

      mass = dilep->M();
    }

    nProbes+=1; // count probe total
    
    // Fill tree
    data.mass	= mass;
    data.pt 	= sc2->Pt();
    data.eta	= sc2->Eta();
    data.phi	= sc2->Phi();
    data.weight	= 1;
    data.q	= q2;
    data.npv	= npv;
    data.pass	= (pass) ? 1 : 0;

    outTree->Fill();
    
  }  
  delete infile;
  infile=0, intree=0;	   

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectProbesEleEff"); 
}
