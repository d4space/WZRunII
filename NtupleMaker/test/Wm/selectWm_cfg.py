import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/a/arapyan/public/forCatherine/miniAOD-prod_PAT.root'
	'file:/u/user/sangilpark/WorkDir/WZ_13TeV/CMSSW_7_4_1_patch1/src/WZ/CSA14/test/datasample/DYJetsToLNu_M50_13TeV_bx25sample_MINIAODSIM.root'
#    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('selectWm.root')
	)

process.demo = cms.EDAnalyzer('selectWm',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.p = cms.Path(process.demo)
