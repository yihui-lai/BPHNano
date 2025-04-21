import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *

BDh = cms.EDProducer("BDhFitter",    
   # which beamSpot to reference
   beamSpot = cms.InputTag('offlineBeamSpot'),
   vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
   tracks          = cms.InputTag("packedPFCandidates"),
   lostTracks      = cms.InputTag("lostTracks"),
   genParticle = cms.InputTag('finalGenParticlesBPH'),

   tkNHitsCut = cms.int32(3), # Number of valid hits on track
   tkPtCut = cms.double(0.3), # Pt of track
   tkChi2Cut = cms.double(30.), # Track normalized Chi2
   mPiPiCut = cms.double(0.7), # invariant mass of track pair, assuming both tracks are charged pions
   vtxChi2Cut = cms.double(20), # Vertex KLM chi2
   vtxDecaySigXYCut = cms.double(3), # KLM XY decay distance significance
   vtxDecaySigXYZCut = cms.double(0.1), # KLM XYZ decay distance significance
   cosThetaXYCut = cms.double(0.9),  # cos(angleXY) between x and p of V0 candidate
   cosThetaXYZCut = cms.double(-1), # cos(angleXYZ) between x and p of V0 candidate
   kShortMassCut = cms.double(0.05), # V0 mass window +- pdg value

   # Not used for now
   tkIPSigXYCut = cms.double(-1.),
   tkIPSigZCut = cms.double(-1.),
   vtxDecayXYCut = cms.double(-1.),
   ssVtxDecayXYCut = cms.double(-1.),
   allowSS = cms.bool(False),
   innerOuterTkDCAThreshold = cms.double(5.),
   innerTkDCACut = cms.double(1.),
   outerTkDCACut = cms.double(1.),
   allowWideAngleVtx = cms.bool(False),
   innerHitPosCut = cms.double(4.),
)


# Ks0
Ks0CandidatesTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('BDh','SelectedV0Collection'),
    cut = cms.string(""),
    name = cms.string("Kshort"),
    doc = cms.string("Kshort Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = cms.PSet(        
        CandVars,
        # track match
        Trk1_idx = uint('Trk1_idx'),
        Trk1_pt  = ufloat('Trk1_pt'),
        Trk1_eta  = ufloat('Trk1_eta'),
        Trk1_phi  = ufloat('Trk1_phi'),
        Trk1_dca_beamspot = ufloat("Trk1_dca_beamspot"),
        Trk1_dcaErr_beamspot = ufloat("Trk1_dcaErr_beamspot"),
        Trk1_PixelHits  = ufloat('Trk1_PixelHits'),
        Trk1_bestTrackHits  = ufloat('Trk1_bestTrackHits'),
        #Trk1_pseudoTrackHits  = ufloat('Trk1_pseudoTrackHits'),
        Trk1_ipsigXY  = ufloat('Trk1_ipsigXY'),
        Trk1_ipsigZ  = ufloat('Trk1_ipsigZ'),
        Trk1_normalizedChi2  = ufloat('Trk1_normalizedChi2'),

        Trk2_idx = uint('Trk2_idx'),
        Trk2_pt  = ufloat('Trk2_pt'),
        Trk2_eta  = ufloat('Trk2_eta'),
        Trk2_phi  = ufloat('Trk2_phi'),
        Trk2_dca_beamspot = ufloat("Trk2_dca_beamspot"),
        Trk2_dcaErr_beamspot = ufloat("Trk2_dcaErr_beamspot"),
        Trk2_PixelHits  = ufloat('Trk2_PixelHits'),
        Trk2_bestTrackHits  = ufloat('Trk2_bestTrackHits'),
        #Trk2_pseudoTrackHits  = ufloat('Trk2_pseudoTrackHits'),
        Trk2_ipsigXY  = ufloat('Trk2_ipsigXY'),
        Trk2_ipsigZ  = ufloat('Trk2_ipsigZ'),
        Trk2_normalizedChi2  = ufloat('Trk2_normalizedChi2'),
        # closest xing point
        trk1_dot_trk2 = ufloat('trk1_dot_trk2'),
        cxPtx = ufloat('cxPtx'),
        cxPty = ufloat('cxPty'),
        cxPtz = ufloat('cxPtz'),
        cxPtR2 = ufloat('cxPtR2'),
        dca = ufloat('dca'),
        massSquared = ufloat('massSquared'),
        # fit
        KLM_vtx_x  = ufloat('KLM_vtx_x'),
        KLM_vtx_y  = ufloat('KLM_vtx_y'),
        KLM_vtx_z  = ufloat('KLM_vtx_z'),
        KLM_chi2  = ufloat('KLM_chi2'),
        KLM_ndof = ufloat('KLM_ndof'),
        KLM_normalizedChi2 = ufloat('KLM_normalizedChi2'),
        KLM_cos_theta_XY = ufloat('KLM_cos_theta_XY'),
        KLM_cos_theta_XYZ = ufloat('KLM_cos_theta_XYZ'),
        #KLM_Trk1_pt  = ufloat('KLM_Trk1_pt'),
        #KLM_Trk1_eta  = ufloat('KLM_Trk1_eta'),
        #KLM_Trk1_phi  = ufloat('KLM_Trk1_phi'),
        #KLM_Trk2_pt  = ufloat('KLM_Trk2_pt'),
        #KLM_Trk2_eta  = ufloat('KLM_Trk2_eta'),
        #KLM_Trk2_phi  = ufloat('KLM_Trk2_phi'),
        KLM_Ks0_pt = ufloat('KLM_Ks0_pt'),
        KLM_Ks0_eta = ufloat('KLM_Ks0_eta'),
        KLM_Ks0_phi = ufloat('KLM_Ks0_phi'),
        KLM_Ks0_mass   = ufloat('KLM_Ks0_mass'),
        KLM_distMagXY       = ufloat("KLM_distMagXY"),
        KLM_sigmaDistMagXY  = ufloat("KLM_sigmaDistMagXY"),
        KLM_distMagXYZ      = ufloat("KLM_distMagXYZ"),
        KLM_sigmaDistMagXYZ = ufloat("KLM_sigmaDistMagXYZ"),
        # fit
        Kin_vtx_x = ufloat('Kin_vtx_x'),
        Kin_vtx_y = ufloat('Kin_vtx_y'),
        Kin_vtx_z = ufloat('Kin_vtx_z'),
        Kin_chi2 = ufloat('Kin_chi2'),
        Kin_dof = ufloat('Kin_dof'),
        Kin_prob = ufloat('Kin_prob'),
        Kin_pt  = ufloat('Kin_pt'),
        Kin_eta = ufloat('Kin_eta'),
        Kin_phi = ufloat('Kin_phi'),
        Kin_mass = ufloat('Kin_mass'),
        Kin_massErr = ufloat('Kin_massErr'),
        Kin_cos_theta_2D = ufloat('Kin_cos_theta_2D'),
        Kin_fitted_cos_theta_2D = ufloat('Kin_fitted_cos_theta_2D'),
        Kin_l_xy = ufloat('Kin_l_xy'),
        Kin_l_xy_unc = ufloat('Kin_l_xy_unc'),
        Kin_vtx_cxx = ufloat('Kin_vtx_cxx'),
        Kin_vtx_cyy = ufloat('Kin_vtx_cyy'),
        Kin_vtx_czz = ufloat('Kin_vtx_czz'),
        Kin_vtx_cyx = ufloat('Kin_vtx_cyx'),
        Kin_vtx_czx = ufloat('Kin_vtx_czx'),
        Kin_vtx_czy = ufloat('Kin_vtx_czy'),
        Kin_trk1_pt =  ufloat('Kin_trk1_pt'),
        Kin_trk1_eta = ufloat('Kin_trk1_eta'),
        Kin_trk1_phi = ufloat('Kin_trk1_phi'),
        Kin_trk2_pt =  ufloat('Kin_trk2_pt'),
        Kin_trk2_eta = ufloat('Kin_trk2_eta'),
        Kin_trk2_phi = ufloat('Kin_trk2_phi'),
        )
)

CountKs0Candidates = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag('BDh','SelectedV0Collection')
)

Ks0CandidatesMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = Ks0CandidatesTable.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(310),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(1.0),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

Ks0CandidatesMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = Ks0CandidatesTable.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("Ks0CandidatesMCMatch"),
    objName = Ks0CandidatesTable.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("kshort"),
    docString = cms.string("MC matching to status==2 kshort"),
)



# Gen match
genpartTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BDh", "Genmatch"),
    cut       = cms.string(""),
    name      = cms.string("Genmatch"),
    doc       = cms.string("genpart Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        channelFlag= uint('channelFlag'),
        b_numberOfDaughters= uint('b_numberOfDaughters'),
        d0_numberOfDaughters= uint('d0_numberOfDaughters'),
        ks0_numberOfDaughters= uint('ks0_numberOfDaughters'),
        Ks_flight_distance= ufloat('Ks_flight_distance'),
        Ks_flight_distance_2D= ufloat('Ks_flight_distance_2D'),
        idx_bpion= uint('idx_bpion'),
        idx_dstar0= uint('idx_dstar0'),
        idx_dstar_decay= uint('idx_dstar_decay'),
        idx_d0= uint('idx_d0'),
        idx_ks= uint('idx_ks'),
        idx_kspip= uint('idx_kspip'),
        idx_kspim= uint('idx_kspim'),
        idx_pip= uint('idx_pip'),
        idx_pim= uint('idx_pim'),
        idx_kspip_iso03= ufloat('idx_kspip_iso03'),
        idx_kspim_iso03 = ufloat('idx_kspim_iso03'),
        idx_pip_iso03 = ufloat('idx_pip_iso03'),
        idx_pim_iso03 = ufloat('idx_pim_iso03'),
        idx_kspip_pt= ufloat('idx_kspip_pt'),
        idx_kspim_pt = ufloat('idx_kspim_pt'),
        idx_pip_pt = ufloat('idx_pip_pt'),
        idx_pim_pt = ufloat('idx_pim_pt'),
        idx_kspip_dr= ufloat('idx_kspip_dr'),
        idx_kspim_dr = ufloat('idx_kspim_dr'),
        idx_pip_dr = ufloat('idx_pip_dr'),
        idx_pim_dr = ufloat('idx_pim_dr'),
        idx_pi0_1= uint('idx_pi0_1'),
        idx_pi0_2= uint('idx_pi0_2'),
        )
)
Countgenpart = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDh", "Genmatch")
)


# D0
D0Table = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BDh", "D0"),
    cut       = cms.string(""),
    name      = cms.string("D0"),
    doc       = cms.string("D0 Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        # track match
        Ks0_idx = uint('Ks0_idx'),
        Trk3_idx = uint('Trk3_idx'),
        Trk3_pt  = ufloat('Trk3_pt'),
        Trk3_eta  = ufloat('Trk3_eta'),
        Trk3_phi  = ufloat('Trk3_phi'),
        Trk3_dca_beamspot = ufloat("Trk3_dca_beamspot"),
        Trk3_dcaErr_beamspot = ufloat("Trk3_dcaErr_beamspot"),
        Trk3_PixelHits  = ufloat('Trk3_PixelHits'),
        Trk3_bestTrackHits  = ufloat('Trk3_bestTrackHits'),
        #Trk3_pseudoTrackHits  = ufloat('Trk3_pseudoTrackHits'),
        Trk3_ipsigXY  = ufloat('Trk3_ipsigXY'),
        Trk3_ipsigZ  = ufloat('Trk3_ipsigZ'),
        Trk3_normalizedChi2  = ufloat('Trk3_normalizedChi2'),
        Trk4_idx = uint('Trk4_idx'),
        Trk4_pt  = ufloat('Trk4_pt'),
        Trk4_eta  = ufloat('Trk4_eta'),
        Trk4_phi  = ufloat('Trk4_phi'),
        Trk4_dca_beamspot = ufloat("Trk4_dca_beamspot"),
        Trk4_dcaErr_beamspot = ufloat("Trk4_dcaErr_beamspot"),
        Trk4_PixelHits  = ufloat('Trk4_PixelHits'),
        Trk4_bestTrackHits  = ufloat('Trk4_bestTrackHits'),
        #Trk4_pseudoTrackHits  = ufloat('Trk4_pseudoTrackHits'),
        Trk4_ipsigXY  = ufloat('Trk4_ipsigXY'),
        Trk4_ipsigZ  = ufloat('Trk4_ipsigZ'),
        Trk4_normalizedChi2  = ufloat('Trk4_normalizedChi2'),
        # closest xing point
        trk3_dot_trk4 = ufloat('trk3_dot_trk4'),
        cxPtx = ufloat('cxPtx'),
        cxPty = ufloat('cxPty'),
        cxPtz = ufloat('cxPtz'),
        cxPtR2 = ufloat('cxPtR2'),
        dca = ufloat('dca'),
        # fit
        Kin_vtx_x = ufloat('Kin_vtx_x'),
        Kin_vtx_y = ufloat('Kin_vtx_y'),
        Kin_vtx_z = ufloat('Kin_vtx_z'),
        Kin_chi2 = ufloat('Kin_chi2'),
        Kin_dof = ufloat('Kin_dof'),
        Kin_prob = ufloat('Kin_prob'),
        Kin_cos_theta_2D = ufloat('Kin_cos_theta_2D'),
        Kin_l_xy = ufloat('Kin_l_xy'),
        Kin_l_xy_unc = ufloat('Kin_l_xy_unc'),
        Kin_trk3_pt =  ufloat('Kin_trk3_pt'),
        Kin_trk3_eta = ufloat('Kin_trk3_eta'),
        Kin_trk3_phi = ufloat('Kin_trk3_phi'),
        Kin_trk4_pt =  ufloat('Kin_trk4_pt'),
        Kin_trk4_eta = ufloat('Kin_trk4_eta'),
        Kin_trk4_phi = ufloat('Kin_trk4_phi'),
        Kin_Ks0_pt = ufloat('Kin_Ks0_pt'),
        Kin_Ks0_eta = ufloat('Kin_Ks0_eta'),
        Kin_Ks0_phi = ufloat('Kin_Ks0_phi'),
        Kin_D0_pt = ufloat('Kin_D0_pt'),
        Kin_D0_eta = ufloat('Kin_D0_eta'),
        Kin_D0_phi = ufloat('Kin_D0_phi'),
        Kin_D0_mass   = ufloat('Kin_D0_mass'),
        Ks_flight_distance = ufloat('Ks_flight_distance'),
        Ks_flight_distance_2D = ufloat('Ks_flight_distance_2D'),
        )
)

CountD0 = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDh", "D0")
)

BDhMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = D0Table.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(421),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(1.0),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)


BDhMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = D0Table.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("BDhMCMatch"),
    objName = D0Table.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("D0"),
    docString = cms.string("MC matching to status==2 D0"),
)


# Pion Tracks 
PionTrackTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src  = cms.InputTag("BDh:SelectedTracks"),
    cut  = cms.string(""),
    name = cms.string("Track"),
    doc  = cms.string("track collection"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        vx              = Var("vx()", float, doc="x coordinate of vtx position [cm]"),
        vy              = Var("vy()", float, doc="y coordinate of vtx position [cm]"),
        vz              = Var("vz()", float, doc="z coordinate of vtx position [cm]"),
        bt_pt           = Var("userFloat('bt_pt')", float, doc="pt of best track [GeV]"),
        trackHighPurity = Var("userFloat('trackHighPurity')", float, doc="trackHighPurity"),
        ipsigXY_bs         = Var("userFloat('ipsigXY_bs')", float, doc="ipsigXY"),
        ipsigZ_bs          = Var("userFloat('ipsigZ_bs')", float, doc="ipsigZ"),
        ipsigXY_pv         = Var("userFloat('ipsigXY_pv')", float, doc="ipsigXY"),
        ipsigZ_pv          = Var("userFloat('ipsigZ_pv')", float, doc="ipsigZ"),
        ipXY_bs         = Var("userFloat('ipXY_bs')", float, doc="ipXY"),
        ipZ_bs          = Var("userFloat('ipZ_bs')", float,  doc="ipZ"),
        ipXY_pv         = Var("userFloat('ipXY_pv')", float, doc="ipXY"),
        ipZ_pv          = Var("userFloat('ipZ_pv')", float,  doc="ipZ"),
        # User variables defined in plugins/TrackMerger.cc
        dz      = Var("userFloat('dz')", float, doc="dz signed wrt first PV [cm]"),
        dxy     = Var("userFloat('dxy')", float, doc="dxy (with sign) wrt first PV [cm]"),
        dzS     = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV [cm]"),
        dxyS    = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV [cm]"),
        nValidHits      = Var("userInt('nValidHits')", int, doc="Number of valid hits"),
        ptErr      = Var("userFloat('ptErr')", float, doc="Pt uncertainty"),
        normChi2   = Var("userFloat('normChi2')", float, doc="Track fit chi-squared divided by n.d.o.f."),
        nValidPixelHits = Var("userInt('nValidPixelHits')", float, doc="Number of pixel hits"),
        ),
)

PionTrackMCMatch = cms.EDProducer("MCMatcher",              # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = PionTrackTable.src,                        # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(211),                     # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
    checkCharge = cms.bool(False),                          # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                            # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),                         # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                          # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                 # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
)


PionTrackMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects   = PionTrackMCMatch.src,
    genParts      = cms.InputTag("finalGenParticlesBPH"),
    mcMap         = cms.InputTag("PionTrackMCMatch"),
    objName       = PionTrackTable.name,
    objType       = PionTrackTable.name,
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("track"),
    docString     = cms.string("MC matching to status==1 pions"),
)



BDhSequence = cms.Sequence(BDh + Ks0CandidatesMCMatch + BDhMCMatch + PionTrackMCMatch)
BDhSequenceTable = cms.Sequence(Ks0CandidatesTable + D0Table + Ks0CandidatesMCTable + PionTrackTable + BDhMCTable + genpartTable + PionTrackMCTable)


