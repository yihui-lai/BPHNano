import sys
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
# source /cvmfs/cms.cern.ch/common/crab-setup.sh

config = Configuration()
#config = config()

config.section_("General")
config.General.requestName = 'BDh_NanoPost_2022_Data_may14_1'
config.General.workArea = '/afs/cern.ch/work/y/yilai/gamma/crab_projects_data_may14'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "test_data_2022.py"
#config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script_data.sh'
config.JobType.scriptArgs = ['isMC=1','era=UL2018','dataRun=X','isVjets=0']
config.JobType.inputFiles = ['BDh_postproc.py', 'BDh_postproc_data.py', 'BDh_Producer.py', 'test_mc_2022.py', 'test_data_2022.py']
config.JobType.outputFiles = ['test_data_Skim.root']
#config.JobType.sendPythonFolder = True
config.section_("Data")
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-22Sep2023-v1/MINIAOD'
config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022G-22Sep2023-v1/MINIAOD'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100000
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 3
#config.Data.totalUnits = 5
config.JobType.maxMemoryMB = 2000  ## 2500*4
#config.JobType.maxJobRuntimeMin = 1315  ## 21.9 hours
config.JobType.numCores = 2
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Muon.json'

config.Data.outLFNDirBase = '/store/user/yilai/NanoPost'
#config.Data.outLFNDirBase = '/store/group/phys_bph/gamma/ntuples/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.section_("Site")
#config.Site.storageSite = "T2_CH_CERN"

config.Site.storageSite = "T3_US_FNALLPC"

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
#        #config.Data.unitsPerJob = 2000000
#        n+=1
#        nnn="%s"%n
#        config.General.requestName = "BDh_NanoPost_2022_v1_"+dataset.split('/')[1][:30]+dataset.split('/')[2][:30]+nnn
#        config.Data.outputDatasetTag = dataset.split('/')[2][:30]+nnn
#        #crabCommand('submit', config = config)


