// -*- C++ -*-
//
// Package:    Test/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Catherine Aiko Medlock
//         Created:  Tue, 22 Jul 2014 11:26:17 GMT
// This code is just an adaptation of Kevin Sung's original code used for the 8 TeV analysis:
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectZmm.C

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TVector2.h>               // 2D vector class

//
// class declaration
//

class selectZmm : public edm::EDAnalyzer {
   public:
      explicit selectZmm(const edm::ParameterSet&);
      ~selectZmm();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;
enum { eMuMu2HLT=1, eMuMu1HLT1L1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };  // event category enum
//
// static data member definitions
//

//--------------------------------------------------------------------------------------------------------------
// Settings 
//============================================================================================================== 

const Double_t MASS_LOW  = 40;
const Double_t MASS_HIGH = 200;
const Double_t PT_CUT    = 25;
const Double_t ETA_CUT   = 2.4;
const Double_t MUON_MASS = 0.105658369;

const Int_t BOSON_ID  = 23;
const Int_t LEPTON_ID = 13;

Double_t nsel=0;

//TString outFilename = TString("selectZmm.root");
//TFile *outFile_Zmm = new TFile();
TTree *outTree = new TTree();
edm::Service<TFileService> fs_Zmm;

//
// Declare output ntuple variables
//
UInt_t 	runNum=0, lumiSec=0, evtNum=0;
UInt_t  matchGen=0;
UInt_t 	category=0;
UInt_t  npv=0, npu=0;
UInt_t  id_1, id_2;
Double_t x_1, x_2, xPDF_1, xPDF_2;	// Todo 
Double_t scalePDF, weightPDF;		// Todo	
LorentzVector *genV=0;			// Todo
Float_t genVPdgID=0, genVPt=0, genVPhi=0, genVy=0, genVMass=0;
Float_t scale1fb; 			// Todo
Float_t rawpfmet=0, rawpfmetPhi=0;
Float_t type1pfmet=0, type1pfmetPhi=0;
Float_t genmet=0, genmetPhi=0;
Float_t u1=0, u2=0;
Int_t q1=0, q2=0;
LorentzVector *dilep=0, *lep1=0, *lep2=0;
///// muon specific /////
Float_t trkIso1=0, trkIso2=0;		// Todo
Float_t emIso1=0, emIso2=0;		// Todo
Float_t hadIso1=0, hadIso2=0;		// Todo
Float_t pfChIso1=0, pfChIso2=0;
Float_t pfNhIso1=0, pfNhIso2=0;
Float_t pfPhotonIso1=0, pfPhotonIso2=0;
Float_t pfCombIso1=0, pfCombIso2=0;	// Todo
Float_t d01, dz1, d02, dz2; 		// Todo
Float_t muNchi21, muNchi22;		// Todo
UInt_t nPixHits1, nTkLayers1, nPixHits2, nTkLayers2;	// Todo
UInt_t nValidHits1, nMatch1, nValidHits2, nMatch2;	// Todo
UInt_t typeBits1, typeBits2; 		// Todo
Int_t isLooseMuon1=0, isSoftMuon1=0, isTightMuon1=0, isMediumMuon1=0;
Int_t isLooseMuon2=0, isSoftMuon2=0, isTightMuon2=0, isMediumMuon2=0;
Int_t passSingleMuTrigger=0, matchTrigObj1=0, matchTrigObj2=0;

//
// constructors and destructor
//
selectZmm::selectZmm(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
{
   //now do what ever initialization is needed
}


selectZmm::~selectZmm()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectZmm::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // Good vertex requirement
   if (vertices->empty()) return; // Skip the event if no PV found
   const reco::Vertex &PV = vertices->front();
   npv = vertices->size();

   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   if(muons->size()<2) return; // Skip the event if there is no possibility of both tag and probe leptons

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   //const TString singleMu("HLT_IsoMu20_eta2p1_IterTrk02_v1");
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
     //if((const TString)names.triggerName(i)!=singleMu) continue;
     //
     //==================================	
     //Trigger Wildcard selection
     //==================================	
     //
     //==================================	
     //Select Triggers starting from "HLT_IsoMu20_v*"
     //==================================	
     std::string trig = names.triggerName(i);
     bool acceptHLT=false;
     if(trig.find("HLT_IsoMu20_v")!=std::string::npos) acceptHLT=true;
     if(!acceptHLT) continue;
     //
     //==================================	
     //Check selected triggers
     //==================================	
     std::cout<<"Selected Trigger: " <<trig<<std::endl;
     passSingleMuTrigger = triggerBits->accept(i) ? 1 : 0;
     break;
   }

   //
   // SELECTION PROCEDURE:
   //  (1) Find a good muon that passes the soft ID -> this will be the "tag"
   //  (2) Pair the tag with a probe muon so that that tag+probe mass is
   //      inside the Z window
   //
   Bool_t foundTag = kFALSE;
   
   for (unsigned int i1=0;i1 < muons->size();i1++) {
     const pat::Muon &tag = (*muons)[i1];

     Int_t isLooseMuon = tag.isLooseMuon() ? 1 : 0;
     if(!isLooseMuon) continue; // lepton selection

     foundTag = kTRUE;

     // Tag lepton information
     LorentzVector vTag(tag.pt(),tag.eta(),tag.phi(),MUON_MASS);

     lep1 = &vTag;

     q1 = tag.charge();
     pfChIso1 = tag.chargedHadronIso();
     pfNhIso1 = tag.neutralHadronIso();
     pfPhotonIso1 = tag.photonIso();
     isLooseMuon1 = tag.isLooseMuon() ? 1 : 0;
     isSoftMuon1  = tag.isSoftMuon(PV) ? 1 : 0;
     isTightMuon1 = tag.isTightMuon(PV) ? 1 : 0;
     isMediumMuon1 = tag.isMediumMuon() ? 1 : 0;

     // Match to a trigger object
     matchTrigObj1 = 0;
     Int_t nObj = 0, tagTrigObjIdx = triggerObjects->size();
     if(passSingleMuTrigger) {
       for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
           matchTrigObj1 = (sqrt((tag.eta()-obj.eta())*(tag.eta()-obj.eta())+(tag.phi()-obj.phi())*(tag.phi()-obj.phi())) <= 0.2) ? 1 : 0;
           if(matchTrigObj1) tagTrigObjIdx = nObj;
           nObj++;
       }
     }

     for (unsigned int i2=0;i2 < muons->size();i2++) {
       if(i2==i1) continue;
       const pat::Muon &probe = (*muons)[i2];

       // Probe lepton information
       q2 = probe.charge();
       pfChIso2 = probe.chargedHadronIso();
       pfNhIso2 = probe.neutralHadronIso();
       pfPhotonIso2 = probe.photonIso();
       isLooseMuon2 = probe.isLooseMuon() ? 1 : 0;
       isSoftMuon2  = probe.isSoftMuon(PV) ? 1 : 0;
       isTightMuon2 = probe.isTightMuon(PV) ? 1 : 0;
       isMediumMuon2 = probe.isMediumMuon() ? 1 : 0;

       LorentzVector vProbe(probe.pt(),probe.eta(),probe.phi(),MUON_MASS);

       lep2 = &vProbe;

       // Match to a trigger object
       matchTrigObj2 = 0;
       if(passSingleMuTrigger) {
         Int_t nProbeObj = 0;
         for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
             if(matchTrigObj1 && nProbeObj==tagTrigObjIdx) continue; // Do not want to match the tag and the probe to the same object
             matchTrigObj2 = (sqrt((probe.eta()-obj.eta())*(probe.eta()-obj.eta())+(probe.phi()-obj.phi())*(probe.phi()-obj.phi())) <= 0.2) ? 1 : 0;
             nProbeObj++;
         }
       }

       // Dilepton 4-vector
       LorentzVector vDilep = vTag + vProbe;
       dilep = &vDilep;

       // Perform matching of leptons to generator level leptons
       const reco::GenParticle* gen1 = tag.genParticle();
       const reco::GenParticle* gen2 = probe.genParticle();
       Int_t id1, id2;
       Float_t eta1, eta2;
       Float_t phi1, phi2;
       Bool_t match1, match2;
       Bool_t hasGenMatch = kFALSE;
       if(gen1!=NULL && gen2!=NULL) {
         id1 = gen1->pdgId(), id2 = gen2->pdgId();
         eta1 = gen1->eta(), eta2 = gen2->eta();
         phi1 = gen1->phi(), phi2 = gen2->phi();
         match1 = ( fabs(id1)==13 && sqrt((tag.eta()-eta1)*(tag.eta()-eta1)+(tag.phi()-phi1)*(tag.phi()-phi1)) < 0.5 );
         match2 = ( fabs(id2)==13 && sqrt((probe.eta()-eta2)*(probe.eta()-eta2)+(probe.phi()-phi2)*(probe.phi()-phi2)) < 0.5 );
         if(match1 && match2) hasGenMatch = kTRUE;
       }
       matchGen  = hasGenMatch ? 1 : 0;
       id_1  = hasGenMatch ? id1 : -999;
       id_2  = hasGenMatch ? id2 : -999;

       // Type-1 corrected PF MET
       edm::Handle<pat::METCollection> mets;
       iEvent.getByToken(metToken_, mets);

       const pat::MET &met = mets->front();
       TVector2 vtype1pfMET = TVector2(met.px(),met.py());
       type1pfmet    = vtype1pfMET.Mod();
       type1pfmetPhi = vtype1pfMET.Phi();

       // Generator level MET
       TVector2 vgenMET = TVector2(met.genMET()->px(),met.genMET()->py());
       genmet    = vgenMET.Mod();
       genmetPhi = vgenMET.Phi();

       // Raw PF MET
       edm::Handle<pat::PackedCandidateCollection> pfs;
       iEvent.getByToken(pfToken_, pfs);
       Float_t rawpfMETpx=0, rawpfMETpy=0;

       for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
         const pat::PackedCandidate &pf = (*pfs)[jcand];
         rawpfMETpx -= pf.px();
         rawpfMETpy -= pf.py();
       }

       TVector2 vrawpfMET = TVector2(rawpfMETpx,rawpfMETpy);
       rawpfmet    = vrawpfMET.Mod();
        rawpfmetPhi = vrawpfMET.Phi();

       // Generator level lepton information and hadronic recoil
       genV = new LorentzVector(0,0,0,0);
       genVPdgID = 0;
       genVPt    = 0;
       genVPhi   = 0;
       genVy     = 0;
       genVMass  = 0;
       u1        = 0;
       u2        = 0;
       if(gen1!=NULL) {
         const reco::Candidate* mother = gen1->mother(0);
         genVPdgID = mother->pdgId();
         genVPt    = mother->pt();
         genVPhi   = mother->phi();
         genVy     = mother->y();
         genVMass  = mother->mass();
         TVector2 vZPt((vDilep.Pt())*cos(vDilep.Phi()),(vDilep.Pt())*sin(vDilep.Phi()));
         TVector2 vMet(type1pfmet*cos(type1pfmetPhi), type1pfmet*sin(type1pfmetPhi));
         TVector2 vU = -1.0*(vMet+vZPt);
         u1 = ((vDilep.Px())*(vU.Px()) + (vDilep.Py())*(vU.Py()))/(vDilep.Pt()); // u1 = (pT . u)/|pT|
         u2 = ((vDilep.Px())*(vU.Py()) - (vDilep.Py())*(vU.Px()))/(vDilep.Pt()); // u2 = (pT x u)/|pT|
       }
       // Count total number of events selected
       nsel++;
       // Fill tree
       outTree->Fill();
       return;
     } // End of i2 loop
   } // End of i1 loop
   if(!foundTag) return; // Skip event if there is no tag lepton

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
selectZmm::beginJob()
{

  //
  // Set up output ntuple
  //

  //outFile = new TFile(outFilename,"RECREATE");
  //outTree = new TTree("Events","Events");
  outTree = fs_Zmm->make<TTree>("Events", "Events");

  outTree->Branch("runNum",        &runNum,   	   "runNum/i");        // event run number
  outTree->Branch("lumiSec",       &lumiSec,   	   "lumiSec/i");       // event lumi section
  outTree->Branch("evtNum",        &evtNum,   	   "evtNum/i");        // event number
  outTree->Branch("matchGen",      &matchGen,      "matchGen/I");      // event has both leptons matched to MC Z->ll
  outTree->Branch("category",      &category,      "category/I");      // dilepton category
  outTree->Branch("id_1",	   &id_1, 	   "id_1/I");          // PDF info -- parton ID for parton 1
  outTree->Branch("id_2",	   &id_2, 	   "id_2/I");          // PDF info -- parton ID for parton 2
  outTree->Branch("x_1",	   &x_1, 	   "x_1/d");           // PDF info -- x for parton 1
  outTree->Branch("x_2",	   &x_2, 	   "x_2/d");           // PDF info -- x for parton 2
  outTree->Branch("xPDF_1",	   &xPDF_1, 	   "xPDF_1/d");        // PDF info -- x*F for parton 1
  outTree->Branch("xPDF_2",	   &xPDF_2, 	   "xPDF_2/d");        // PDF info -- x*F for parton 2
  outTree->Branch("scalePDF",	   &scalePDF, 	   "scalePDF/d");      // PDF info -- energy scale of parton interaction
  outTree->Branch("weightPDF",	   &weightPDF, 	   "weightPDF/d");     // PDF info -- PDF weight
  outTree->Branch("npv",           &npv,           "npv/I");           // number of vertices
  outTree->Branch("npu",           &npu,           "npu/I");           // number of in-time PU events(MC)
  outTree->Branch("genV", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &genV); // GEN boson 4-vector (signal MC)
  outTree->Branch("genVPt",        &genVPt,        "genVPt/F");        // GEN boson pT (signal MC)
  outTree->Branch("genVPhi",       &genVPhi,       "genVPhi/F");       // GEN boson phi (signal MC)
  outTree->Branch("genVy",         &genVy,         "genVy/F");         // GEN boson rapidity (signal MC)
  outTree->Branch("genVMass",      &genVMass,      "genVMass/F");      // GEN boson mass (signal MC)
  outTree->Branch("scale1fb",      &scale1fb,      "scale1fb/F");      // event weight per 1/fb (MC)
  outTree->Branch("rawpfmet",      &rawpfmet,      "rawpfmet/F");      // Raw PF MET
  outTree->Branch("rawpfmetPhi",   &rawpfmetPhi,   "rawpfmetPhi/F");   // Raw PF MET phi
  outTree->Branch("type1pfmet",    &type1pfmet,    "type1pfmet/F");    // Type-1 corrected PF MET
  outTree->Branch("type1pfmetPhi", &type1pfmetPhi, "type1pfmetPhi/F"); // Type-1 corrected PF MET phi
  outTree->Branch("genmet",        &genmet,        "genmet/F");        // Generator level MET
  outTree->Branch("genmetPhi",     &genmetPhi,     "genmetPhi/F");     // Generator level MET phi
  outTree->Branch("u1",            &u1,            "u1/F");            // parallel component of recoil
  outTree->Branch("u2",            &u2,            "u2/F");            // perpendicular component of recoil
  outTree->Branch("q1",            &q1,            "q1/I");            // tag lepton charge
  outTree->Branch("q2",            &q2,            "q2/I");            // probe lepton charge
  outTree->Branch("dilep","ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &dilep); // dilepton 4-vector
  outTree->Branch("lep1", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep1);  // lepton 4-vector (tag lepton)
  outTree->Branch("lep2", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep2);  // lepton 4-vector (probe lepton)
  outTree->Branch("passSingleMuTrigger", &passSingleMuTrigger, "passSingleMuTrigger/I"); // single muon trigger
  outTree->Branch("matchTrigObj1", &matchTrigObj1,  "matchTrigObj1/I"); // tag lepton match to trigger object
  outTree->Branch("matchTrigObj2", &matchTrigObj2,  "matchTrigObj2/I"); // probe lepton match to trigger object
  ///// muon specific ///// 
  outTree->Branch("trkIso1",       &trkIso1,      "trkIso1/F");       // track isolation of tag lepton
  outTree->Branch("trkIso2",       &trkIso2,      "trkIso2/F");       // track isolation of probe lepton
  outTree->Branch("emIso1",        &emIso1,       "emIso1/F");        // ECAL isolation of tag lepton
  outTree->Branch("emIso2",        &emIso2,       "emIso2/F");        // ECAL isolation of probe lepton
  outTree->Branch("hadIso1",       &hadIso1,      "hadIso1/F");       // HCAL isolation of tag lepton
  outTree->Branch("hadIso2",       &hadIso2,      "hadIso2/F");       // HCAL isolation of probe lepton
  outTree->Branch("pfChIso1",      &pfChIso1,     "pfChIso1/F");      // tag lepton charged hadron isolation
  outTree->Branch("pfChIso2",      &pfChIso2,     "pfChIso2/F");      // probe lepton charged hadron isolation
  outTree->Branch("pfNhIso1",      &pfNhIso1,     "pfNhIso1/F");      // tag lepton neutral hadron isolation
  outTree->Branch("pfNhIso2",      &pfNhIso2,     "pfNhIso2/F");      // probe lepton neutral hadron isolation
  outTree->Branch("pfPhotonIso1",  &pfPhotonIso1, "pfPhotonIso1/F");  // tag lepton photon isolation
  outTree->Branch("pfPhotonIso2",  &pfPhotonIso2, "pfPhotonIso2/F");  // probe lepton photon isolation
  outTree->Branch("pfCombIso1",    &pfCombIso1,   "pfCombIso1/F");    // PF combined isolation of tag lepton
  outTree->Branch("pfCombIso2",    &pfCombIso2,   "pfCombIso2/F");    // PF combined isolation of probe lepton
  outTree->Branch("isLooseMuon1",  &isLooseMuon1, "isLooseMuon1/I");  // tag lepton loose muon ID
  outTree->Branch("isSoftMuon1",   &isSoftMuon1,  "isSoftMuon1/I");   // tag lepton soft muon ID
  outTree->Branch("isTightMuon1",  &isTightMuon1, "isTightMuon1/I");  // tag lepton tight muon ID
  outTree->Branch("isMediumMuon1", &isMediumMuon1,"isMediumMuon1/I"); // tag lepton medium muon ID
  outTree->Branch("isLooseMuon2",  &isLooseMuon2, "isLooseMuon2/I");  // probe lepton loose muon ID
  outTree->Branch("isSoftMuon2",   &isSoftMuon2,  "isSoftMuon2/I");   // probe lepton soft muon ID
  outTree->Branch("isTightMuon2",  &isTightMuon2, "isTightMuon2/I");  // probe lepton tight muon ID
  outTree->Branch("isMediumMuon2", &isMediumMuon2,"isMediumMuon2/I"); // probe lepton medium muon ID
  outTree->Branch("d01",	   &d01,	  "d01/F");	      // transverse impact parameter of tag lepton
  outTree->Branch("d02",	   &d02,	  "d02/F");	      // transverse impact parameter of probe lepton
  outTree->Branch("dz1",	   &dz1,	  "dz1/F");	      // longitudinal impact parameter of tag lepton
  outTree->Branch("dz2",	   &dz2,	  "dz2/F");	      // longitudinal impact parameter of probe lepton
  outTree->Branch("muNchi21",	   &muNchi21,	  "muNchi21/F");      // muon fit normalized chi^2 of tag lepton
  outTree->Branch("muNchi22",	   &muNchi22,	  "muNchi22/F");      // muon fit normalized chi^2 of probe lepton
  outTree->Branch("nPixHits1",     &nPixHits1,    "nPixHits1/i");     // number of pixel hits of tag muon
  outTree->Branch("nPixHits2",     &nPixHits2,    "nPixHits2/i");     // number of pixel hits of probe muon
  outTree->Branch("nTkLayers1",    &nTkLayers1,   "nTkLayers1/i");    // number of tracker layers of tag muon
  outTree->Branch("nTkLayers2",    &nTkLayers2,   "nTkLayers2/i");    // number of tracker layers of probe muon
  outTree->Branch("nMatch1",       &nMatch1,      "nMatch1/i");       // number of matched segments of tag muon
  outTree->Branch("nMatch2",       &nMatch2,      "nMatch2/i");       // number of matched segments of probe muon
  outTree->Branch("nValidHits1",   &nValidHits1,  "nValidHits1/i");   // number of valid muon hits of tag muon
  outTree->Branch("nValidHits2",   &nValidHits2,  "nValidHits2/i");   // number of valid muon hits of probe muon
  outTree->Branch("typeBits1",     &typeBits1,    "typeBits1/i");     // muon type of tag muon
  outTree->Branch("typeBits2",     &typeBits2,    "typeBits2/i");     // muon type of probe muon
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectZmm::endJob() 
{
   // Save tree in output file
   //outFile->Write();
   //outFile->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "Z -> m m" << std::endl;
  std::cout << nsel << " +/- " << sqrt(nsel) << " per 1/fb" << std::endl;
  std::cout << std::endl;

  //std::cout << std::endl;
  //std::cout << "  <> Output saved in " << outFilename << "/" << std::endl;
  //std::cout << std::endl;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectZmm::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectZmm::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectZmm::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectZmm::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectZmm::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectZmm);
