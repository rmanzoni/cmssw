datasets = { 
  "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" : '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
  'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : '/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
  'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : '/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
  'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' : '/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM',
  'DYJetsToLL_M-50_HT-600toinf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' :'/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM' 
}
sample=4
import sys

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = ''

## MC 
config.General.workArea = list(datasets)[int(sample)]
                                                
#config.General.workArea = 'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Asympt25ns_MCRUN2_74_V9'

## DATA

config.section_('JobType')
config.JobType.psetName    = 'runFrameworkMC.py'
config.JobType.pluginName  = 'Analysis'

config.JobType.inputFiles = ['../data/weightfile.root']
## MC
config.JobType.pyCfgParams = ['globalTag=74X_mcRun2_asymptotic_v2']

#config.JobType.inputFiles  = ['PY8_RunIISpring15DR74_bx50_MC_PFCHS.db','PY8_RunIISpring15DR74_bx50_MC_Puppi.db']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

## MC
config.Data.inputDataset = datasets[list(datasets)[int(sample)]]

config.Data.unitsPerJob = 1
config.Data.inputDBS  = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000
config.Data.publication = False


config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_DE_DESY'
config.Data.outLFNDirBase = '/store/user/rfriese/mvamet/skimming/2016-01-13/'
config.General.workArea = '/nfs/dust/cms/user/rfriese/crab_mvamet_skim-2016-01-13'


#  LocalWords:  MINIAODSIM
