from CRABClient.UserUtilities import config
config = config()


config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/B0ToKstarJpsi'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "test_cfg.py" 
config.JobType.allowUndistributedCMSSW = True


config.section_('Data')
config.Data.inputDataset = '/BdtoJpsiKstar_Jpsito2Mu_KstartoKPi_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv3-124X_mcRun3_2022_realistic_postEE_v1-v2/MINIAODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/group/phys_bphys/cpv/B0ToKstarJpsi_9_9_24'
config.Data.outputDatasetTag = 'B0ToKstarJpsi_9_9_24'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'


