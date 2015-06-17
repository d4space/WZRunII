import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        '/store/user/arapyan/WplusToENu_CT10_13TeV-powheg-pythia8/Spring14dr-PU-S14-POSTLS170_V6-v1/140723_170241/0000/miniAOD-prod_PAT_19.root'
'file:/u/user/sangilpark/miniAODSample/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PU20bx25.root'
#    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.selectZmm = cms.EDAnalyzer('selectZmm',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.p = cms.Path(process.selectZmm)
