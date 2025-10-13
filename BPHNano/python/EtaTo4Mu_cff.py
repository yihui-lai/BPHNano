import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

########################### Selections ###########################
EtaTo4Mu = cms.EDProducer(
    'EtaTo4MuBuilder',
    src = cms.InputTag('muonBPH', 'AllMuons'),
    transientTracksSrc = cms.InputTag('muonBPH', 'AllTransientMuons'),
    #src = cms.InputTag('muonBPH', 'SelectedMuons'),
    #transientTracksSrc = cms.InputTag('muonBPH', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 4 && abs(eta) < 2.4 && isMediumMuon && isGlobalMuon'),
    lep2Selection = cms.string('pt > 3 && abs(eta) < 2.4 && isLooseMuon && isGlobalMuon'),
    lep3Selection = cms.string('pt > 2 && abs(eta) < 2.4 && isLooseMuon && isGlobalMuon'),
    lep4Selection = cms.string('pt > 2 && abs(eta) < 2.4 && isLooseMuon && isGlobalMuon'),
    preVtxSelection  = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1.'
        '&& abs(userCand("l1").vz - userCand("l3").vz) <= 1.'
        '&& abs(userCand("l1").vz - userCand("l4").vz) <= 1.'
        '&& pt > 1. && ((mass > 0.45 && mass < 0.6)||(mass > 0.9 && mass < 1.0))'
        '&& charge() == 0'
    ),
    postVtxSelection = cms.string('userFloat("sv_prob") > 0.0 && userFloat("fitted_mass") > 0.45 && userFloat("fitted_mass") < 1.2'),
)

CountEtaTo4MuonBPH = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("EtaTo4Mu:Selected4Leptons")
)

########################### Tables ###########################
EtaTo4MuTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("EtaTo4Mu:Selected4Leptons"),
    cut = cms.string(""),
    name = cms.string("EtaTo4Mu"),
    doc  = cms.string("Four-muon collections"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        fitted_mass = Var("userFloat('fitted_mass')", float, doc="Fitted four-lepton mass"),
        fitted_massErr = Var("userFloat('fitted_massErr')", float, doc="Fitted four-lepton mass error"),
        svprob = Var("userFloat('sv_prob')", float, doc="Vertex fit probability"),
        vtx_x = Var("userFloat('vtx_x')", float, doc="Vertex position x"),
        vtx_y = Var("userFloat('vtx_y')", float, doc="Vertex position y"),
        vtx_z = Var("userFloat('vtx_z')", float, doc="Vertex position z"),
        dca_avg = Var("userFloat('dca_avg')", float, doc="Average DCA between lepton tracks"),

        lep_min_deltaR = Var("userFloat('lep_min_deltaR')", float, doc="Minimum ΔR between the 4 muons"),
        lep_max_deltaR = Var("userFloat('lep_max_deltaR')", float, doc="Maximum ΔR between the 4 muons"),
        lep_avg_deltaR = Var("userFloat('lep_avg_deltaR')", float, doc="Average ΔR between the 4 muons"),

        l1_idx = Var("userInt('l1_idx')", int, doc="Index of lepton 1"),
        l2_idx = Var("userInt('l2_idx')", int, doc="Index of lepton 2"),
        l3_idx = Var("userInt('l3_idx')", int, doc="Index of lepton 3"),
        l4_idx = Var("userInt('l4_idx')", int, doc="Index of lepton 4"),
    )
)

########################### MC Matching ###########################

EtaTo4MuBPHMCMatch = cms.EDProducer(
    "MCMatcher",
    src         = EtaTo4MuTable.src,
    matched     = cms.InputTag("finalGenParticlesBPH"),
    mcPdgId     = cms.vint32(221, 331),
    checkCharge = cms.bool(False),
    mcStatus    = cms.vint32(2),
    maxDeltaR   = cms.double(0.05),
    maxDPtRel   = cms.double(0.5),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True),
)

EtaTo4MuBPHMCTable = cms.EDProducer(
    "CandMCMatchTableProducerBPH",
    recoObjects = EtaTo4MuTable.src,
    genParts    = cms.InputTag("finalGenParticlesBPH"),
    mcMap       = cms.InputTag("EtaTo4MuBPHMCMatch"),
    objName     = EtaTo4MuTable.name,
    objType     = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("EtaTo4Mu"),
    docString   = cms.string("MC matching to status==2 4-muon resonance"),
)

########################### Sequences ###########################

EtaTo4MuSequence = cms.Sequence(EtaTo4Mu)
EtaTo4MuTables = cms.Sequence(EtaTo4MuTable)
EtaTo4MuMCSequence = cms.Sequence(EtaTo4Mu + EtaTo4MuBPHMCMatch)
EtaTo4MuMCTables = cms.Sequence(EtaTo4MuTable + EtaTo4MuBPHMCTable)

