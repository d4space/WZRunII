import FWCore.ParameterSet.Config as cms

process = cms.Process("table")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      'file:/u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_2_0/src/WZ/CSA14/miniAODsample/WZJetsTo3LNu_Tune4C_PU20bx25/484D51C6-2673-E411-8AB0-001E67398412.root')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
    maxEventsToPrint = cms.untracked.int32(3),
    printVertex = cms.untracked.bool(False),
    printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
    src = cms.InputTag("prunedGenParticles"),
    )

process.p = cms.Path(process.printTree)

