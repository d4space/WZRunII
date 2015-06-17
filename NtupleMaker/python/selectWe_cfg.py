import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/mc/Spring14miniaod/WminusToENu_CT10_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/020492A6-161A-E411-AF9E-00266CFE79A4.root',
#      'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/mc/Spring14miniaod/WminusToENu_CT10_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/66EE9628-D918-E411-B779-1CC1DE1D023A.root',
 
#'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/mc/Spring14miniaod/WplusToENu_CT10_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1E86DC0F-612D-E411-B6A7-001E6739753A.root',

'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02215B44-2D70-E411-90A3-0025905A60B8.root',
#        '/store/user/arapyan/WplusToENu_CT10_13TeV-powheg-pythia8/Spring14dr-PU-S14-POSTLS170_V6-v1/140723_170241/0000/miniAOD-prod_PAT_19.root'
#    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.demo = cms.EDAnalyzer('selectWe',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electrons = cms.InputTag("slimmedElectrons"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.p = cms.Path(process.demo)
