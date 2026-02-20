import sys


# if using a list
# source /cvmfs/cms.cern.ch/common/crab-setup.sh
# python3  xx.py xx.txt
#from CRABClient.UserUtilities import config
#config = config()

# if using this file directly
from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.requestName = 'NanoPost_2022_MC_2026Feb17_BsToJPsiPhi'
config.General.workArea = '/afs/cern.ch/work/y/yilai/gamma/crab_projects_MC_'+config.General.requestName
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "test_mc_2022.py"
#config.JobType.psetName = "test_mc_test.py"
#config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.scriptArgs = ['isMC=1','era=UL2018','dataRun=X','isVjets=0']
config.JobType.inputFiles = ['BDh_postproc.py', 'BDh_Producer.py', 'test_mc_2022.py', 'test_mc_test.py']
config.JobType.outputFiles = ['test_mc_Skim.root']
#config.JobType.sendPythonFolder = True
config.section_("Data")
#config.Data.inputDataset = '/BuToD0K_D0ToKs2Pi_Run3/yilai-Run3Summer22_MiniAODv4-5443ed9e0a49f9c5d5f0b2fff4804347/USER'
#config.Data.inputDataset = '/LambdaBToJpsiLambda_JpsiFilter_MuFilter_LambdaFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.Data.inputDataset = '/LambdaBToJpsiLambda_Unbiased_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
#config.Data.inputDataset = '/lambdab_lambda2K_2/yilai-Run3Summer22_MiniAODv4-9542347e91aacc6744c44ee907c98ff2/USER'
#config.Data.inputDataset = '/lambdab_lambda2Pi/yilai-Run3Summer22_MiniAODv4-9542347e91aacc6744c44ee907c98ff2/USER'
#config.Data.inputDataset = '/JPsiMuMu_JPsiNoFilter_2MuPtEtaFilter_TuneCP5_13p6TeV-pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v4/MINIAODSIM'
#config.Data.inputDataset = '/BuToJpsiK_JpsiToMuMu_MuFilter_Pt-2_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
#config.Data.inputDataset = '/Jpsito2Mu_JpsiPT8_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-MUO_POG_130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
#config.Data.inputDataset = '/Upsilonto2Mu_UpsilonFilter_2MuFilter_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
#config.Data.inputDataset = '/QCDB-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v3/MINIAODSIM'
#config.Data.inputDataset = '/BuToJpsiKstar_Unbiased_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.Data.inputDataset = '/ButoJpsiK_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.Data.inputDataset = '/LambdabtoJpsiKp_Jpsito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.Data.inputDataset = '/BsToJPsiPhi_JPsiToMuMu_PhiToKK_EtaPtFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM'
#config.Data.inputDataset = '/QCDB-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.General.requestName = 'NanoPost_2022_MC_2026Feb18_QCDB_4Jets_HT_100to200'

#config.Data.inputDataset = '/QCDB-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.General.requestName = 'NanoPost_2022_MC_2026Feb18_QCDB_4Jets_HT_200to400'
#
#config.Data.inputDataset = '/QCDB-4Jets_HT-40to100_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'
#config.General.requestName = 'NanoPost_2022_MC_2026Feb18_QCDB_4Jets_HT_40to100'

#config.Data.inputDataset = '/QCDB-4Jets_HT-40to100_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'
#config.Data.inputDataset = '/QCDB-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v4/MINIAODSIM'
#config.Data.inputDataset = '/QCDB-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM'

#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100000
#config.Data.totalUnits = 5
config.JobType.maxMemoryMB = 2000  ## 2500*4
#config.JobType.maxJobRuntimeMin = 1315  ## 21.9 hours
config.JobType.numCores = 2
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 10000

#config.Data.outLFNDirBase = '/store/user/yilai/NanoPost'
config.Data.outLFNDirBase = '/store/group/phys_bphys/yilai/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.storageSite = "T3_US_FNALLPC"

# config.section_("User")
#config.User.voGroup = 'dcms'

#if __name__ == '__main__':
#    f=open(sys.argv[1]) 
#    content = f.readlines()
#    content = [x.strip() for x in content] 
#    from CRABAPI.RawCommand import crabCommand
#    n=100
#    for dataset in content :
#        config.Data.inputDataset = dataset
#        n+=1
#        nnn="%s"%n
#        config.General.requestName = "BDh_NanoPost_2022_v1_"+dataset.split('/')[1][:30]+dataset.split('/')[2][:30]+nnn
#        config.Data.outputDatasetTag = dataset.split('/')[2][:30]+nnn
#        #print('submit', config.Data.inputDataset, config.General.requestName, config.Data.outputDatasetTag)
#        crabCommand('submit', config = config)
#

