
//
// Mode
//

enum Mode {
    MuonTagAndProbe         =1,
    MuonFakeRate            =2,
    ElectronTagAndProbe     =3,
    ElectronFakeRate        =4,
    MuonTagAndProbeMC       =5,
    ElectronTagAndProbeMC   =6
};

//
// Egamma cut based
//

    enum CutType {
        DETAIN          = (1<<0),
        DPHIIN          = (1<<1),
        SIGMAIETAIETA   = (1<<2),
        HOE             = (1<<3),
        OOEMOOP         = (1<<4),
        D0VTX           = (1<<5),
        DZVTX           = (1<<6),
        ISO             = (1<<7),
        VTXFIT          = (1<<8),
        MHITS           = (1<<9)
    };

// all possible cuts pass
static const unsigned int PassAll         = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | OOEMOOP | D0VTX | DZVTX | ISO | VTXFIT | MHITS;
static const unsigned int PassNoIso       = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | OOEMOOP | D0VTX | DZVTX | VTXFIT | MHITS;

//
// tag and probe / FR options
//

enum Option {

    // Muon selection
    MuonIsoDenominator=10,
    MuonIDDenominator=11,
    MuonIDIsoDenominator=12,
    MuonIDIsoRndTagDenominator=13,
    MuonIDIsoRndTagDblNoDzDenominator=14,

    MuonIsoNumerator=20,
    MuonIDNumerator=21,
    MuonTrigSglNumerator=22,
    MuonTrigDblLeadNumerator=23,
    MuonTrigDblTrailNumerator=24,
    MuonTrigDblNumerator=25,

    // Electron selection
    ElectronIsoDenominator=30,
    ElectronIDDenominator=31,
    ElectronIDIsoDenominator=32,
    ElectronIDIsoRndTagDenominator=33,
    ElectronIDIsoRndTagDblNoDzDenominator=34,

    ElectronIsoNumerator=40,
    ElectronIDNumerator=41,
    ElectronTrigSglNumerator=42,
    ElectronTrigDblLeadNumerator=43,
    ElectronTrigDblTrailNumerator=44,
    ElectronTrigDblNumerator=45

};

//
// triggers
//

struct TriggerResults {

    // singles
    UInt_t HLT_IsoMu24_eta2p1_tag_;
    UInt_t HLT_IsoMu24_eta2p1_probe_;
    UInt_t HLT_Ele27_WP80_tag_;
    UInt_t HLT_Ele27_WP80_probe_;

    // doubles
    UInt_t HLT_Ele17_Ele8_LeadingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_TrailingLeg_probe_;
    UInt_t HLT_Ele17_Ele8_probe_;
    UInt_t HLT_Mu17_Mu8_LeadingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_TrailingLeg_probe_;
    UInt_t HLT_Mu17_Mu8_probe_;

};

//
// map branches
//

void mapExtraBranches(LeptonTree *leptonTree, TriggerResults &triggerResults) 
{
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_tag", &triggerResults.HLT_IsoMu24_eta2p1_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_IsoMu24_eta2p1_probe", &triggerResults.HLT_IsoMu24_eta2p1_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_tag", &triggerResults.HLT_Ele27_WP80_tag_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele27_WP80_probe", &triggerResults.HLT_Ele27_WP80_probe_);

    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_LeadingLeg_probe", &triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_TrailingLeg_probe", &triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Ele17_Ele8_probe", &triggerResults.HLT_Ele17_Ele8_probe_);

    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_LeadingLeg_probe", &triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_TrailingLeg_probe", &triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_);
    leptonTree->tree_->SetBranchAddress("HLT_Mu17_Mu8_probe", &triggerResults.HLT_Mu17_Mu8_probe_);

}

//
// probes
//

bool isProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode) 
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    // muons
    if (option == MuonIsoDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                   return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsHPASS) != LeptonTree::PassMuIsHPASS) return false;
    }
    if (option == MuonIDDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                   return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ 
                + leptonTree.pfnhiso03_ - 0.5 * leptonTree.dbeta03_))/leptonTree.probe_.Pt() > 0.10)     return false;
    }

    if (option == MuonIDIsoDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                   return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_
                + leptonTree.pfnhiso03_ - 0.5 * leptonTree.dbeta03_))/leptonTree.probe_.Pt() > 0.10)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsHPASS) != LeptonTree::PassMuIsHPASS) return false;
    }

    if (option == MuonIDIsoRndTagDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                   return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_
                + leptonTree.pfnhiso03_ - 0.5 * leptonTree.dbeta03_))/leptonTree.probe_.Pt() > 0.10)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsHPASS) != LeptonTree::PassMuIsHPASS) return false;
    }

    if (option == MuonIDIsoRndTagDblNoDzDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_IsoMu24_eta2p1_tag_ == 0 && !isMC)                                   return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if ((triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) && !isMC)         return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_
                + leptonTree.pfnhiso03_ - 0.5 * leptonTree.dbeta03_))/leptonTree.probe_.Pt() > 0.10)     return false;
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsHPASS) != LeptonTree::PassMuIsHPASS) return false;
    }


    // electrons
    if (option == ElectronIsoDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                       return false;
        if ((leptonTree.mediumId_ & PassNoIso) != PassNoIso)                                        return false;
        if (leptonTree.mhit_ != 0)                                                                  return false;
        if (!leptonTree.chargesAgree_)                                                              return false;     
        if (fabs(leptonTree.sceta_) > 1.4442 && fabs(leptonTree.sceta_) < 1.566)                    return false;        
        if (fabs(leptonTree.sceta_) < 1.4442 && leptonTree.hoe_ > 0.1)                              return false;
        if (fabs(leptonTree.sceta_) < 1.556  && leptonTree.hoe_ > 0.075)                            return false;
    }
    if (option == ElectronIDDenominator) {
        if (leptonTree.tag_.Pt() < 30.0)                                                            return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                                       return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ + leptonTree.pfnhiso03_ 
                - leptonTree.ea03_ * TMath::Max(double(0.0), leptonTree.rhoIsoAllCentral_)))/leptonTree.probe_.Pt() > 0.09) return false;
        if (fabs(leptonTree.sceta_) > 1.4442 && fabs(leptonTree.sceta_) < 1.566)                    return false;
    }

    if (option == ElectronIDIsoDenominator) {        
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if ((leptonTree.mediumId_ & PassNoIso) != PassNoIso)                                        return false;
        if (leptonTree.mhit_ != 0)                                                                  return false;
        if (!leptonTree.chargesAgree_)                                                              return false;
        if (fabs(leptonTree.sceta_) > 1.4442 && fabs(leptonTree.sceta_) < 1.566)                    return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ + leptonTree.pfnhiso03_
                - leptonTree.ea03_ * TMath::Max(double(0.0), leptonTree.rhoIsoAllCentral_)))/leptonTree.probe_.Pt() > 0.09) return false;
    }

    if (option == ElectronIDIsoRndTagDenominator) 
    {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;    
        if ((leptonTree.mediumId_ & PassNoIso) != PassNoIso)                                        return false;
        if (leptonTree.mhit_ != 0)                                                                  return false;
        if (!leptonTree.chargesAgree_)                                                              return false;
        if (fabs(leptonTree.sceta_) < 1.4442 && leptonTree.hoe_ > 0.1)                              return false;
        if (fabs(leptonTree.sceta_) < 1.556  && leptonTree.hoe_ > 0.075)                            return false;
        if (fabs(leptonTree.sceta_) > 1.4442 && fabs(leptonTree.sceta_) < 1.566)                    return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ + leptonTree.pfnhiso03_
                - leptonTree.ea03_ * TMath::Max(double(0.0), leptonTree.rhoIsoAllCentral_)))/leptonTree.probe_.Pt() > 0.09) return false;
    }

    if (option == ElectronIDIsoRndTagDblNoDzDenominator) {
        if (leptonTree.tag_.Pt() < 32.0)                                           return false;
        if (triggerResults.HLT_Ele27_WP80_tag_ == 0 && !isMC)                      return false;
        if ((leptonTree.qTag_ < 0 && leptonTree.rnd_ < 0.5) ||
            (leptonTree.qTag_ > 0 && leptonTree.rnd_ >= 0.5))                      return false;
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0 && !isMC)        return false;
        if ((leptonTree.mediumId_ & PassNoIso) != PassNoIso)                                        return false;
        if (leptonTree.mhit_ != 0)                                                                  return false;
        if (!leptonTree.chargesAgree_)                                                              return false;
        if (fabs(leptonTree.sceta_) < 1.4442 && leptonTree.hoe_ > 0.1)                              return false;
        if (fabs(leptonTree.sceta_) < 1.556  && leptonTree.hoe_ > 0.075)                            return false;
        if (fabs(leptonTree.sceta_) > 1.4442 && fabs(leptonTree.sceta_) < 1.566)                    return false;
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ + leptonTree.pfnhiso03_
                - leptonTree.ea03_ * TMath::Max(double(0.0), leptonTree.rhoIsoAllCentral_)))/leptonTree.probe_.Pt() > 0.09) return false;
    }
 
    return true;
}

bool isPassProbe(unsigned int option, const LeptonTree &leptonTree, const TriggerResults &triggerResults, const Mode mode) 
{

    bool isMC = false;
    if (mode == ElectronTagAndProbeMC || mode == MuonTagAndProbeMC) isMC = true;

    // muons
    if (option == MuonIsoNumerator) {
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_
                + leptonTree.pfnhiso03_ - 0.5 * leptonTree.dbeta03_))/leptonTree.probe_.Pt() > 0.10)     return false;
    }
    if (option == MuonIDNumerator) {
        if ((leptonTree.leptonSelection_ & LeptonTree::PassMuIsHPASS) != LeptonTree::PassMuIsHPASS) return false;
    }

    if (option == MuonTrigSglNumerator) {
        if (triggerResults.HLT_IsoMu24_eta2p1_probe_ == 0 && !isMC)   return false;
    }

    if (option == MuonTrigDblLeadNumerator) {
        if ((triggerResults.HLT_Mu17_Mu8_LeadingLeg_probe_ == 0) && !isMC) return false;
    }

    if (option == MuonTrigDblTrailNumerator) {
        if ((triggerResults.HLT_Mu17_Mu8_TrailingLeg_probe_ == 0) && !isMC) return false;
    }

    if (option == MuonTrigDblNumerator) {
        if ((triggerResults.HLT_Mu17_Mu8_probe_ == 0) && !isMC) return false;
    }

    // electrons
    if (option == ElectronIsoNumerator) {
        if ((leptonTree.pfchiso03_ + TMath::Max(double(0.0), leptonTree.pfemiso03_ + leptonTree.pfnhiso03_
                - leptonTree.ea03_ * TMath::Max(double(0.0), leptonTree.rhoIsoAllCentral_)))/leptonTree.probe_.Pt() > 0.09) return false;
    }
    if (option == ElectronIDNumerator) {
        if ((leptonTree.mediumId_ & PassNoIso) != PassNoIso)                                        return false;
        if (leptonTree.mhit_ != 0)                                                                  return false;
        if (!leptonTree.chargesAgree_)                                                              return false;
    }

    if (option == ElectronTrigSglNumerator) {
        if (triggerResults.HLT_Ele27_WP80_probe_ == 0 && !isMC)   return false;
    }

    if (option == ElectronTrigDblLeadNumerator) {
        if (triggerResults.HLT_Ele17_Ele8_LeadingLeg_probe_ == 0 && !isMC) return false;
    }

    if (option == ElectronTrigDblTrailNumerator) {
        if (triggerResults.HLT_Ele17_Ele8_TrailingLeg_probe_ == 0 && !isMC) return false;
    }

    if (option == ElectronTrigDblNumerator) {
        if (triggerResults.HLT_Ele17_Ele8_probe_ == 0 && !isMC) return false;
    }

    return true;
}

