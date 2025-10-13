import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *


# Used to study lambdaB -> lambda0(->\p\pi) + 2h
BtoD0Kstar = cms.EDProducer("BtoD0KstarProducer",    
   # which beamSpot to reference
   beamSpot        = cms.InputTag('offlineBeamSpot'),
   vertices        = cms.InputTag('offlineSlimmedPrimaryVertices'),
   tracks          = cms.InputTag("packedPFCandidates"),
   lostTracks      = cms.InputTag("lostTracks"),
   # all Tracks
   tkNHitsCut      = cms.int32(3),    # Number of valid hits on track
   minTrackPt      = cms.double(1.0),
   maxTrackEta     = cms.double(2.4), # Eta of track
   tkChi2Cut       = cms.double(30.), # Track normalized Chi2
   tkIPSigXYCut    = cms.double(0.5), # Track IP significance
   # Ks0
   TrkSigXYCut       = cms.double(4),
   vtxChi2Cut        = cms.double(6.63),
   vtxDecaySigXYCut  = cms.double(-1),                # (Not used) KLM XY decay distance significance
   vtxDecaySigXYZCut = cms.double(-1),                # (Not used) KLM XYZ decay distance significance
   cosThetaXYCut     = cms.double(0.995),
   cosThetaXYZCut    = cms.double(-1),                # cos(angleXYZ) between x and p of V0 candidate
   Ks0_l_xyzSigCut   = cms.double(1),                 # Ks flight distancesignificance from D0
   # D0
   DtkPtCut            = cms.double(1.0),   # Pt cut of track 3, 4
   diTrack2_dca        = cms.double(0.2),    # 0.2cm is a decent cut to remove comb.bkg
   Trk34SigXYCut       = cms.double(0.5),    # track DCA significance
   D0_PtCut            = cms.double(3), # D0 Pt cut, GeV
   D0vtxDecaySigXYCut  = cms.double(2), # D0 XY distance significance from PV
   # B
   BtkPtCut        = cms.double(0.9), # Pt cut of track 5
   B_PtCut         = cms.double(4), # B Pt cut, GeV
   Btrk_dcaSigCut  = cms.double(1), # B trk dca Sig
   verbose         = cms.int32(0)
)

D0Table = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BtoD0Kstar", "D0"),
    cut       = cms.string(""),
    name      = cms.string("D0"),
    doc       = cms.string("D0 Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        )
)
D0MCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = D0Table.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(421),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)
D0MCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = D0Table.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("D0MCMatch"),
    objName = D0Table.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("D0"),
    docString = cms.string("MC matching to status==2 D0"),
)

KstarTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BtoD0Kstar", "Kstar"),
    cut       = cms.string(""),
    name      = cms.string("Kstar"),
    doc       = cms.string("Kstar Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        )
)
KstarMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = KstarTable.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(323),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)
KstarMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = KstarTable.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("KstarMCMatch"),
    objName = KstarTable.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("Kstar"),
    docString = cms.string("MC matching to status==2 Kstar"),
)

BTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BtoD0Kstar", "B"),
    cut       = cms.string(""),
    name      = cms.string("B"),
    doc       = cms.string("B Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        rawmass = ufloat('rawmass'),
        massErr = ufloat('massErr'),
        chi2 = ufloat('chi2'),
        prob = ufloat('prob'),
        lxyz = ufloat('lxyz'),
        lxyzSig = ufloat('lxyzSig'),
        MC_pt = ufloat('MC_pt'),
        MC_eta = ufloat('MC_eta'),
        MC_phi = ufloat('MC_phi'),
        MC_mass = ufloat('MC_mass'),
        MC_massErr = ufloat('MC_massErr'),
        MC_chi2 = ufloat('MC_chi2'),
        MC_prob = ufloat('MC_prob'),
        MC_lxyz = ufloat('MC_lxyz'),
        MC_lxyzSig = ufloat('MC_lxyzSig'),
        d0_idx = uint('d0_idx'),
        kstar_idx = uint('kstar_idx'),
        D0_mass = ufloat('D0_mass'),
        Kstar_mass = ufloat('Kstar_mass'),
        )
)
BMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = BTable.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(521),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)
BMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = BTable.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("BMCMatch"),
    objName = BTable.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("B"),
    docString = cms.string("MC matching to status==2 B"),
)

# Gen match
BDKstarGen = cms.EDProducer("BDKstarGen",
   genParticle = cms.InputTag('finalGenParticlesBPH'),
)
BDKstarGenmatchTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BDKstarGen", "BDKstarGenmatch"),
    cut       = cms.string(""),
    name      = cms.string("BGenmatch"),
    doc       = cms.string("genpart Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        idx_B    = uint('idx_B'),
        B_charge = uint('B_charge'),
        idx_D0 = uint("idx_D0"),
        idx_D0_dau1 = uint("idx_D0_dau1"),
        idx_D0_dau2 = uint("idx_D0_dau2"),
        idx_Kstar = uint("idx_Kstar"),
        idx_Kstar_pi = uint("idx_Kstar_pi"),
        idx_Ks  = uint("idx_Ks"),
        idx_Ks_pi1 = uint("idx_Ks_pi1"),
        idx_Ks_pi2 = uint("idx_Ks_pi2")
        )
)

Countgenpart = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDKstarGen", "BDKstarGenmatch")
)


BDKstarSequence = cms.Sequence(BtoD0Kstar)
BDKstarSequenceTable = cms.Sequence(D0Table + KstarTable + BTable)

BDKstarSequenceMC = cms.Sequence(BDKstarGen + BtoD0Kstar + D0MCMatch + KstarMCMatch + BMCMatch)
BDKstarSequenceMCTable = cms.Sequence(BDKstarGenmatchTable + D0Table + KstarTable + BTable + D0MCTable + KstarMCTable + BMCTable)



