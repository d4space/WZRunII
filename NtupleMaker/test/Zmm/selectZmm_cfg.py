import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Demo")

options = VarParsing.VarParsing("analysis")
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        '/store/user/arapyan/WplusToENu_CT10_13TeV-powheg-pythia8/Spring14dr-PU-S14-POSTLS170_V6-v1/140723_170241/0000/miniAOD-prod_PAT_19.root'
	'file:/d3/scratch/khakim/RunII/miniAODsamples/Spring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_amcatnloFXFX_Asympt25ns_MCRUN2_74_V9_v3/5A752801-7514-E511-9EB2-00259073E4C8.root'

## DYJetsToLL 13TeV bx50 dataset-----------------------------------------
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/38D1C54C-0F02-E511-A54E-AC853D9F5256.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/3E8D4650-0F02-E511-A60B-0026189438E0.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/B80DA399-0F02-E511-8DA8-0025901D4B06.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/BED3C659-0F02-E511-A55D-002354EF3BE1.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/C6A67F57-0F02-E511-9545-00261894387A.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/D2F31955-0F02-E511-8C03-0025905A60F4.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/DEAB7C6F-0F02-E511-BF61-0025B3E05D8C.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/F07D4B4E-0F02-E511-B7F3-0025905A6118.root',
#       'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/FE850B5A-0F02-E511-9E51-0025905A613C.root'

# DYJetsToLL 13TeV bx50 dataset----------------------------------------
#	'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/38D1C54C-0F02-E511-A54E-AC853D9F5256.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/3E8D4650-0F02-E511-A60B-0026189438E0.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/B80DA399-0F02-E511-8DA8-0025901D4B06.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/BED3C659-0F02-E511-A55D-002354EF3BE1.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/C6A67F57-0F02-E511-9545-00261894387A.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/D2F31955-0F02-E511-8C03-0025905A60F4.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/DEAB7C6F-0F02-E511-BF61-0025B3E05D8C.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/F07D4B4E-0F02-E511-B7F3-0025905A6118.root',
#        'root://xrootd-cms.infn.it//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/10000/FE850B5A-0F02-E511-9E51-0025905A613C.root'

#    fileNames = cms.untracked.vstring(options.inputFiles
    )
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('selectZmm.root')
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
