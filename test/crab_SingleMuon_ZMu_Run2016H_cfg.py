from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'SingleMuonRun2016H_ZMu_PromptReco_v10'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'RunTree_collisions_cfg.py'
config.JobType.outputFiles = ['DTNtuple.root']

config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2016H-ZMu-PromptReco-v2/RAW-RECO'

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = '282408-999990' 

config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 150  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/user/battilan/DTNtuples/ZMu/POS/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

