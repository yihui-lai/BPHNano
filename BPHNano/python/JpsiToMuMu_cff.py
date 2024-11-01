import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

JpsiToMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonBPH', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonBPH', 'SelectedTransientMuons'),
    lep1Selection = cms.string(''),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1.'
                                 ' && 2.8 < mass() && mass() <3.4 && charge() == 0'
                                 ' && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection =  cms.string('2.9<userFloat("fitted_mass")'
                                   ' && userFloat("fitted_mass") < 3.3'
                                   ' && userFloat("sv_prob")>0.001')
)

CountDiMuonBPH = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("JpsiToMuMu:SelectedDiLeptons")
)  


JpsiToMuMuTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("JpsiToMuMu:SelectedDiLeptons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("JpsiToMuMu"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
          fitted_mass = Var("userFloat('fitted_mass')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
          svprob = Var("userFloat('sv_prob')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
          vtx_x =Var("userFloat('vtx_x')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
          vtx_y = Var("userFloat('vtx_y')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
          vtx_z = Var("userFloat('vtx_z')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),

    )
)

JpsiToMuMuSequence = cms.Sequence(JpsiToMuMu)
JpsiToMuMuTables = cms.Sequence(JpsiToMuMuTable)

