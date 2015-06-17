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
//         Created:  Tue, 22 Jul 2014 11:26:17 Gmt_We
// This code is just an adaptation of Kevin Sung's original code used for the 8 TeV analysis:
// https://github.com/jaylawhorn/mitewk/blob/master/Selection/selectWe.C

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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

class selectWe : public edm::EDAnalyzer {
   public:
      explicit selectWe(const edm::ParameterSet&);
      ~selectWe();

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
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
};

//
// constants, enums and typedefs
//

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;

//
// static data member definitions
//

//--------------------------------------------------------------------------------------------------------------
// Settings 
//============================================================================================================== 

const Double_t PT_CUT    = 20;
const Double_t ETA_CUT   = 2.5;
const Double_t ELE_MASS  = 0.000511;

const Double_t ECAL_GAP_LOW  = 1.4442;
const Double_t ECAL_GAP_HIGH = 1.566;

// Count total number of events selected
Double_t nsel_We=0;

TString outFilename = TString("selectWe.root");
TFile *outFile = new TFile();
TTree *outTree_We = new TTree();

//
// Declare output ntuple variables
//
Int_t  npv_We=0;
Float_t genVPdgID_We=0, genVPt_We=0, genVPhi_We=0, genVy_We=0, genVMass_We=0;
Float_t genLepPdgID_We=0, genLepPt_We=0, genLepPhi_We=0;
Float_t rawpfmet_We=0, rawpfmetPhi_We=0;
Float_t type1pfmet_We=0, type1pfmetPhi_We=0;
Float_t genmet_We=0, genmetPhi_We=0;
Float_t mt_We=0, u1_We=0, u2_We=0;
Int_t q_We=0;
LorentzVector *lep_We=0;
///// electron specific /////
Float_t pfChIso_We=0, pfGamIso_We=0, pfNeuIso_We=0;
Int_t isVetoEle_We=0, isLooseEle_We=0, isMediumEle_We=0, isTightEle_We=0;
LorentzVector *sc_We=0;
Int_t passSingleEleTrigger_We=0, matchTrigObj_We=0;

//
// constructors and destructor
//
selectWe::selectWe(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
   metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
   pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
{
   //now do what ever initialization is needed
}


selectWe::~selectWe()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
selectWe::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   // good vertex requirement
   if (vertices->empty()) return; // skip the event if no PV found
   npv_We = vertices->size();
//   const reco::Vertex &PV = vertices->front();

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   if(electrons->size()==0) return; // skip the event if there are no electrons

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   const TString singleEle("HLT_Ele27_eta2p1_WP85_Gsf_v1");
   for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
       if((const TString)names.triggerName(i)!=singleEle) continue;
       passSingleEleTrigger_We = triggerBits->accept(i) ? 1 : 0;
       break;
   }

   //
   // SELECTION PROCEDURE:
   //  (1) Look for 1 good electron matched to trigger
   //  (2) Reject event if another electron is present passing looser cuts
   //
   Int_t  nLooseLep=0;
   Bool_t passSel=kFALSE;
   Int_t  goodEleIdx=0;
   for (unsigned int jElectron=0;jElectron < electrons->size();jElectron++) {
     const pat::Electron &ele = (*electrons)[jElectron];

     // check ECAL gap
     if(fabs(ele.superCluster()->eta())>=ECAL_GAP_LOW && fabs(ele.superCluster()->eta())<=ECAL_GAP_HIGH) continue;

     isLooseEle_We = ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
     isTightEle_We = ele.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

     if(fabs(ele.superCluster()->eta())> 2.5)  continue;    // lepton |eta| cut
     if(ele.superCluster()->energy()   < 20 )  continue;    // lepton pT cut
     if(isLooseEle_We && ele.chargedHadronIso() <= 0.15*ele.pt()) nLooseLep++; // loose lepton selection
     if(nLooseLep>1) { // extra lepton veto
       passSel=kFALSE;
       break;
     }

     if(fabs(ele.superCluster()->eta()) > ETA_CUT) continue; // lepton |eta| cut
     if(ele.superCluster()->energy()    < PT_CUT ) continue; // lepton pT cut
     if(!(isTightEle_We && ele.chargedHadronIso() <= 0.15*ele.pt()))  continue; // lepton selection

     passSel = kTRUE;
     goodEleIdx = jElectron;
   }
   if(passSel==kFALSE) return;
   // Count total number of events selected
   nsel_We++;

   const pat::Electron &goodEle = (*electrons)[goodEleIdx];
   // Lepton information
   q_We           = goodEle.charge();
   pfChIso_We     = goodEle.chargedHadronIso();
   pfGamIso_We    = goodEle.photonIso();
   pfNeuIso_We    = goodEle.neutralHadronIso();
   isVetoEle_We   = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-veto");
   isLooseEle_We  = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-loose");
   isMediumEle_We = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-medium");
   isTightEle_We  = goodEle.electronID("cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");

   // Match to a trigger object
   matchTrigObj_We = 0;
   if(passSingleEleTrigger_We) {
     for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
         matchTrigObj_We = (sqrt((goodEle.eta()-obj.eta())*(goodEle.eta()-obj.eta())+(goodEle.phi()-obj.phi())*(goodEle.phi()-obj.phi())) <= 0.2) ? 1 : 0;
     }
   }

   // Lepton and supercluster 4-vectors
   LorentzVector vLep(goodEle.pt(),goodEle.eta(),goodEle.phi(),ELE_MASS);
   LorentzVector vSC(goodEle.superCluster()->energy()*(goodEle.pt()/goodEle.p()),goodEle.superCluster()->eta(),goodEle.superCluster()->phi(),ELE_MASS);

   lep_We = &vLep;
   sc_We = &vSC;

   // Type-1 corrected PF MET (default)
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);

   const pat::MET &met = mets->front();
   TVector2 vtype1pfMET_We = TVector2(met.px(),met.py());
   type1pfmet_We    = vtype1pfMET_We.Mod();
   type1pfmetPhi_We = vtype1pfMET_We.Phi();

   // Generator level MET
   TVector2 vgenMET_We = TVector2(met.genMET()->px(),met.genMET()->py());
   genmet_We    = vgenMET_We.Mod();
   genmetPhi_We = vgenMET_We.Phi();

   // Raw PF MET
   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   Float_t rawpfMETpx=0, rawpfMETpy=0;

   for(unsigned int jcand=0; jcand < pfs->size(); jcand++) {
     const pat::PackedCandidate &pf = (*pfs)[jcand];
     rawpfMETpx -= pf.px();
     rawpfMETpy -= pf.py();
   }

   TVector2 vrawpfMET_We = TVector2(rawpfMETpx,rawpfMETpy);
   rawpfmet_We    = vrawpfMET_We.Mod();
   rawpfmetPhi_We = vrawpfMET_We.Phi();

   // Lepton transverse mass
   mt_We = sqrt(2.0*vLep.Pt()*vtype1pfMET_We.Mod()*(1.0-cos(vLep.Phi()-vtype1pfMET_We.Phi())));

   // Generator level lepton information and hadronic recoil
   const reco::GenParticle* genLep = goodEle.genLepton();
   if(genLep!=NULL) {
     genLepPdgID_We = genLep->pdgId();
     genLepPt_We    = genLep->pt();
     genLepPhi_We   = genLep->phi();
     const reco::Candidate* mother = genLep->mother(0);
     genVPdgID_We = mother->pdgId();
     genVPt_We    = mother->pt();
     genVPhi_We   = mother->phi();
     genVy_We     = mother->y();
     genVMass_We  = mother->mass();
     TVector2 vWPt(genVPt_We*cos(genVPhi_We),genVPt_We*sin(genVPhi_We));
     TVector2 vLepPt(vLep.Px(),vLep.Py());
     TVector2 vU = -1.0*(vtype1pfMET_We+vLepPt);
     u1_We = (vWPt.Px()*vU.Px()+vWPt.Py()*vU.Py())/genVPt_We; // u1_We = (pT . u)/|pT|
     u2_We = (vWPt.Px()*vU.Px()-vWPt.Py()*vU.Py())/genVPt_We; // u1_We = (pT x u)/|pT|
   }

   // Fill tree
   outTree_We->Fill();

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
selectWe::beginJob()
{

  //
  // Set up output ntuple
  //

  outFile = new TFile(outFilename,"RECREATE");
  outTree_We = new TTree("Events","Events");

  outTree_We->Branch("npv",           &npv_We,           "npv/I");          // number of primary vertices
  outTree_We->Branch("genVPt",        &genVPt_We,        "genVPt/F");       // GEN boson pT (signal MC)
  outTree_We->Branch("genVPhi",       &genVPhi_We,       "genVPhi/F");      // GEN boson phi (signal MC)
  outTree_We->Branch("genVy",         &genVy_We,         "genVy/F");        // GEN boson rapidity (signal MC)
  outTree_We->Branch("genVMass",      &genVMass_We,      "genVMass/F");     // GEN boson mass (signal MC)
  outTree_We->Branch("genLepPt",      &genLepPt_We,      "genLepPt/F");     // GEN lepton pT (signal MC)
  outTree_We->Branch("genLepPhi",     &genLepPhi_We,     "genLepPhi/F");    // GEN lepton phi (signal MC)
  outTree_We->Branch("rawpfmet",      &rawpfmet_We,      "rawpfmet/F");     // Raw PF MET
  outTree_We->Branch("rawpfmetPhi",   &rawpfmetPhi_We,   "rawpfmetPhi/F");  // Raw PF MET phi
  outTree_We->Branch("type1pfmet",    &type1pfmet_We,    "type1pfmet/F");   // Type-1 corrected PF MET
  outTree_We->Branch("type1pfmetPhi", &type1pfmetPhi_We, "type1pfmetPhi/F");// Type-1 corrected PF MET phi
  outTree_We->Branch("genmet",        &genmet_We,        "genmet/F");       // Generator level MET
  outTree_We->Branch("genmetPhi",     &genmetPhi_We,     "genmetPhi/F");    // Generator level MET phi
  outTree_We->Branch("mt",            &mt_We,            "mt/F");           // transverse mass
  outTree_We->Branch("u1",            &u1_We,            "u1/F");           // parallel component of recoil
  outTree_We->Branch("u2",            &u2_We,            "u2/F");           // perpendicular component of recoil 
  outTree_We->Branch("q",             &q_We,             "q/I");            // lepton charge
  outTree_We->Branch("lep", "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &lep_We);   // lepton 4-vector
  outTree_We->Branch("sc",  "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >", &sc_We);    // supercluster 4-vector
  ///// electron specific /////
  outTree_We->Branch("pfChIso",       &pfChIso_We,       "pfChIso/F");      // PF charged hadron isolation of electron
  outTree_We->Branch("pfGamIso",      &pfGamIso_We,      "pfGamIso/F");     // PF photon isolation of electron
  outTree_We->Branch("pfNeuIso",      &pfNeuIso_We,      "pfNeuIso/F");     // PF neutral hadron isolation of electron
  outTree_We->Branch("isVetoEle",     &isVetoEle_We,     "isVetoEle/I");    // lepton veto electron ID
  outTree_We->Branch("isLooseEle",    &isLooseEle_We,    "isLooseEle/I");   // lepton loose electron ID
  outTree_We->Branch("isMediumEle",   &isMediumEle_We,   "isMediumEle/I");  // lepton medium electron ID
  outTree_We->Branch("isTightEle",    &isTightEle_We,    "isTightEle/I");   // lepton tight electron ID
  outTree_We->Branch("passSingleEleTrigger", &passSingleEleTrigger_We, "passSingleEleTrigger/I"); // single electron trigger
  outTree_We->Branch("matchTrigObj",  &matchTrigObj_We,   "matchTrigObj/I"); // lepton match to trigger object
}

// ------------ method called once each job just after ending the event loop  ------------
void 
selectWe::endJob() 
{
   // Save tree to output file
   outFile->Write();
   outFile->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================

  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << "W -> e nu" << std::endl;
  std::cout << " pT > " << PT_CUT << std::endl;
  std::cout << " |eta| < " << ETA_CUT << std::endl;
  std::cout << nsel_We << " +/- " << sqrt(nsel_We) << " per 1/fb" << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "  <> Output saved in " << outFilename << "/" << std::endl;
  std::cout << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
selectWe::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
selectWe::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
selectWe::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
selectWe::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
selectWe::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(selectWe);
