import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/a/arapyan/public/forCatherine/miniAOD-prod_PAT.root'
       #'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/38D1C54C-0F02-E511-A54E-AC853D9F5256.root', # DYJetsToLL sample
       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0EB45900-F9FE-E411-8F52-00266CFFA704.root', # WJetsToLNu sample
    )
)

process.demo = cms.EDAnalyzer("MiniAODTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.p = cms.Path(process.demo)
