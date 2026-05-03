import sys

# if using a list
# source /cvmfs/cms.cern.ch/common/crab-setup.sh
# python3  xx.py xx.txt
from CRABClient.UserUtilities import config
config = config()

# if using this file directly
#from WMCore.Configuration import config
#config = config()


config.section_("General")
requestName = 'BNMC_Feb25'
config.General.requestName = requestName
config.General.workArea = '/afs/cern.ch/work/y/yilai/gamma/'+config.General.requestName
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "run_mc_Run3Summer22.py"
#config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script_mc.sh'
config.JobType.scriptArgs = []
config.JobType.inputFiles = ['postproc_mc.py', "run_mc_Run3Summer22.py", "crab_script_mc.sh"]
config.JobType.outputFiles = ['out_step1_Skim.root']
#config.JobType.sendPythonFolder = True
config.section_("Data")
config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-22Sep2023-v1/MINIAOD'
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 400000 # events
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB = 2000  ## 2500*4
#config.JobType.maxJobRuntimeMin = 1315  ## 21.9 hours
config.JobType.numCores = 1

config.Data.outLFNDirBase = '/store/group/phys_bphys/yilai/eta_2mu2pi/postprocess_v2/2022pre/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.storageSite = "T3_US_FNALLPC"

# config.section_("User")
#config.User.voGroup = 'dcms'

if __name__ == '__main__':
    f=open(sys.argv[1]) 
    content = f.readlines()
    content = [x.strip() for x in content] 
    from CRABAPI.RawCommand import crabCommand
    n=10
    for dataset in content :
        config.Data.inputDataset = dataset
        n+=1
        nnn="%s"%n
        config.General.requestName = requestName +dataset.split('/')[1][:30]+dataset.split('/')[2][:30]+nnn
        config.Data.outputDatasetTag = dataset.split('/')[2][:30]+nnn
        print(config.General.requestName, config.Data.outputDatasetTag)
        crabCommand('submit', config = config)
        #crabCommand('submit', config = config, dryrun = True)

