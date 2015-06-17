# CRAB3 Tutorial : https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#Prerequisites_to_run_the_tutoria
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DYJetsToLL_M50_13TeV_bx50'
config.General.workArea = 'crab_projects'
#config.General.transferOutputs = True
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'selectZmm_cfg.py'

config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-StartupFlat10to50bx50Raw_MCRUN2_74_V8-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased' 	# used for MC
#config.Data.splitting = 'LumiBased'	# used for RD
config.Data.unitsPerJob = 10 		# ex)unitPerJob = 10, 1 job will take 10 file if FileBased

#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'	# used for RD
#config.Data.runRange = '193093-193999' # used for RD

#config.Data.outLFNDirBase = '/store/user/spak/MINIAOD/DYJetsToLL_13TeV/' # or '/store/group/<subdir>'
config.Data.outLFNDirBase = '/store/user/spak/test/DYJetsToLL_13TeV/' # or '/store/group/<subdir>'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

#config.Site.storageSite = <site where the user has permission to write>
#config.Site.whitelist = ['T2_KR_KNU']
config.Site.storageSite = 'T2_KR_KNU'
