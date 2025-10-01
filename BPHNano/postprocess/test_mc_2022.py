# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO --conditions auto:phase1_2022_realistic --datatier NANOAOD --era Run3 --eventcontent NANOAODSIM --filein /store/mc/Run3Summer22MiniAODv4/DstarToD0Pi_D0To2Mu_MuFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v1/30000/71ec4425-d76b-446d-9d89-a7b250c56568.root --fileout file:test_mc.root --nThreads 8 -n 1000 --no_exec --python_filename test_mc.py --scenario pp --step NANO --mc --customise_commands=process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3


from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
#options.register('globalTag', '130X_mcRun3_2022_realistic_v5',
#    #'130X_dataRun3_Prompt_v3',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Global tag"
#)
#
#options.register('outputFiles', 'BPHnano.root',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "outputFiles"
#)
options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Adds gen info/matching"
)
#options.register('skip', 0,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.int,
#    "Skip first N events"
#)
#options.setDefault('maxEvents', 10)
#options.setDefault('tag', 'test')
options.parseArguments()
#
#print("options.inputFiles = ", options.inputFiles)
#print("options.outputFiles = ", options.outputFiles)
#print("options.maxEvents = ", options.maxEvents)
#print("options.skip = ", options.skip)
print("options.isMC = ", options.isMC)

process = cms.Process('NANO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('PhysicsTools.BPHNano.nanoBPH_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/y/yilai/gamma/LambdaBToJpsiLambda_JpsiFilter_MuFilter_LambdaFilter.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECO nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    #SelectEvents = cms.untracked.PSet(
    #    SelectEvents = cms.vstring('filtering_step')
    #),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:test_mc.root'),
    #outputCommands = process.NANOAODSIMEventContent.outputCommands
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
      "keep edmTriggerResults_*_*_*",
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# BPH 
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.BPHNano.nanoBPH_cff import *
process = nanoAOD_customizeMC(process)

#BDK
#process = nanoAOD_customizeBDh_MC(process)
##process.nanoAOD_BPH_step = cms.Path(process.nanoSequence + cms.Sequence(cms.Task(lhcInfoTable)) + cms.Sequence(genWeightsTableTask))


# Lambdab0 -> lambda0 + J/psi
#process = nanoAOD_customizeMuonBPH(process,  True)
#process = nanoAOD_customizeDiMuonBPH(process,True)
process = nanoAOD_customizeTrackBPH(process, True)
#process = nanoAOD_customizeLambda(process,   True)

# Lambdab0 -> lambda0 + hh
process = nanoAOD_customizeLambdahh(process, True)

process.nanoAOD_BPH_step = cms.Path(process.nanoSequence + cms.Sequence(genWeightsTableTask))

# No filter at this point, move to postprocess
#process.BDhFilter = cms.EDFilter("BDhFilter",
#    minNumber = cms.uint32(1),
#    src = cms.InputTag("BDh","B")
#)
#process.BDhFilterSequence = cms.Sequence(process.BDhFilter)
#process.filtering_step = cms.Path(process.BDhFilterSequence)


# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
process.schedule = cms.Schedule(process.nanoAOD_BPH_step,process.endjob_step,process.NANOAODSIMoutput_step)
#process.schedule = cms.Schedule(process.nanoAOD_BPH_step,process.filtering_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 2
process.options.numberOfStreams = 0
# filter all path with the production filter sequence
#for path in process.paths:
#        if not path in ['endjob_step']: continue
#        getattr(process,path).insert(0, process.BDhFilterSequence)


# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon 

#call to customisation function nanoAOD_customizeCommon imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeCommon(process)

# End of customisation functions


# Customisation from command line

process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
