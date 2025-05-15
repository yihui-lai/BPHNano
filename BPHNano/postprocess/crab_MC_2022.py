import sys
#from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
# source /cvmfs/cms.cern.ch/common/crab-setup.sh

#config = Configuration()
config = config()

config.section_("General")
config.General.requestName = 'BDh_NanoPost_2022_MC_may14'
config.General.workArea = '/afs/cern.ch/work/y/yilai/gamma/crab_projects_MC_may14'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "test_mc_2022.py"
#config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.scriptArgs = ['isMC=1','era=UL2018','dataRun=X','isVjets=0']
config.JobType.inputFiles = ['BDh_postproc.py', 'BDh_Producer.py', 'test_mc_2022.py']
config.JobType.outputFiles = ['test_mc_Skim.root']
#config.JobType.sendPythonFolder = True
config.section_("Data")
config.Data.inputDataset = '/BuToD0K_D0ToKs2Pi_Run3/yilai-Run3Summer22_MiniAODv4-5443ed9e0a49f9c5d5f0b2fff4804347/USER'
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 5
config.JobType.maxMemoryMB = 2000  ## 2500*4
#config.JobType.maxJobRuntimeMin = 1315  ## 21.9 hours
config.JobType.numCores = 2
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 10000

config.Data.outLFNDirBase = '/store/user/yilai/NanoPost'
#config.Data.outLFNDirBase = '/store/group/phys_bph/gamma/ntuples/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.section_("Site")
#config.Site.storageSite = "T2_CH_CERN"

config.Site.storageSite = "T3_US_FNALLPC"

# config.section_("User")
#config.User.voGroup = 'dcms'

if __name__ == '__main__':
    f=open(sys.argv[1]) 
    content = f.readlines()
    content = [x.strip() for x in content] 
    from CRABAPI.RawCommand import crabCommand
    n=100
    for dataset in content :
        config.Data.inputDataset = dataset
        n+=1
        nnn="%s"%n
        config.General.requestName = "BDh_NanoPost_2022_v1_"+dataset.split('/')[1][:30]+dataset.split('/')[2][:30]+nnn
        config.Data.outputDatasetTag = dataset.split('/')[2][:30]+nnn
        #print('submit', config.Data.inputDataset, config.General.requestName, config.Data.outputDatasetTag)
        crabCommand('submit', config = config)


