import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

########################### Selections ###########################

MuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonBPH', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonBPH', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.0 && abs(eta) < 2.4 && isLooseMuon && isGlobalMuon'),
    lep2Selection = cms.string('pt > 1.0 && abs(eta) < 2.4 && isLooseMuon && isGlobalMuon'),
    preVtxSelection  = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1.'
                                  '&& 2.9 < mass() && mass() < 3.3 '
                                  '&& charge() == 0'),
    postVtxSelection = cms.string('2.9 < userFloat("fitted_mass") && userFloat("fitted_mass") < 3.3'
                                  '&& userFloat("sv_prob") > 0.001')
)

CountDiMuonBPH = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("MuMu:SelectedDiLeptons")
)  

########################### Tables ###########################

MuMuTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("MuMu:SelectedDiLeptons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("MuMu"),
    doc  = cms.string("Dilepton collections"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
          fitted_mass = Var("userFloat('fitted_mass')", float, doc="Fitted dilepton mass"),
          fitted_massErr = Var("userFloat('fitted_massErr')", float, doc="Fitted dilepton massErr"),
          svprob = Var("userFloat('sv_prob')", float, doc="Vtx fit probability"),
          vtx_x =Var("userFloat('vtx_x')", float, doc="Vtx position in x"),
          vtx_y = Var("userFloat('vtx_y')", float, doc="Vtx position in y"),
          vtx_z = Var("userFloat('vtx_z')", float, doc="Vtx position in z"),
          lep_deltaR = Var("userFloat('lep_deltaR')", float, doc="lep_deltaR"),
          l1_idx = Var("userInt('l1_idx')", int, doc="l1_idx"),
          l2_idx = Var("userInt('l2_idx')", int, doc="l2_idx"),
    )
)

MuMuSequence = cms.Sequence(MuMu)
MuMuTables = cms.Sequence(MuMuTable)
