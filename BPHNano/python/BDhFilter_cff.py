import FWCore.ParameterSet.Config as cms
BDhFilter = cms.EDFilter("BDhFilter",
    src = cms.InputTag("BDh","B"),
    minNumber = cms.uint32( 1 )
)

BDhFilterSequence   = cms.Sequence(BDhFilter)
