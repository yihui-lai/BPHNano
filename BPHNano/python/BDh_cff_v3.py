import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *


xgboost_models = [
    ('Run2022-20250504-1855-Event0', '2022'),                     # include pT
]


savetrack=False

BDh = cms.EDProducer("BDhFitter_v3",    
   # which beamSpot to reference
   xgboost_models = cms.vstring(),
   xgboost_variable_names = cms.vstring(),
   beamSpot        = cms.InputTag('offlineBeamSpot'),
   vertices        = cms.InputTag('offlineSlimmedPrimaryVertices'),
   tracks          = cms.InputTag("packedPFCandidates"),
   lostTracks      = cms.InputTag("lostTracks"),
   # Tracks
   tkNHitsCut      = cms.int32(3),    # Number of valid hits on track
   tkPtCut         = cms.double(0.5), # Pt cut of track 1, 2
   tkEtaCut        = cms.double(2.4), # Eta of track
   tkChi2Cut       = cms.double(30.), # Track normalized Chi2
   tkIPSigXYCut    = cms.double(0.5), # Track IP significance
   # diTrack 1
   TrkSigXYCut       = cms.double(4),                 # track DCA significance
   vtxChi2Cut        = cms.double(6.63),              # Vertex KLM chi2
   vtxDecaySigXYCut  = cms.double(-1),                # (Not used) KLM XY decay distance significance
   vtxDecaySigXYZCut = cms.double(-1),                # (Not used) KLM XYZ decay distance significance
   cosThetaXYCut     = cms.double(0.995),             # cos(theta3D)  between x(D0->Ks0) and p(Ks0) candidate, a loose cut
   cosThetaXYZCut    = cms.double(-1),                # cos(angleXYZ) between x and p of V0 candidate
   # diTrack 2
   DtkPtCut          = cms.double(0.6),   # Pt cut of track 3, 4
   diTrack2_dca      = cms.double(0.2),    # 0.2cm is a decent cut to remove comb.bkg
   Trk34SigXYCut     = cms.double(0.5),    # track DCA significance
   # reco ks0
   Ks0_l_xyzSigCut = cms.double(3), # Ks flight distancesignificance from D0
   # D0
   D0_PtCut            = cms.double(3), # D0 Pt cut, GeV 
   D0vtxDecaySigXYCut  = cms.double(2), # D0 XY distance significance from PV
   # B
   BtkPtCut    = cms.double(0.9), # Pt cut of track 5
   B_PtCut             = cms.double(4), # B Pt cut, GeV
   Btrk_dcaSigCut      = cms.double(1), # B trk dca Sig
   # mass
   kShortMassCut   = cms.double(0.03), # Ks mass window +- pdg value
   D0MassCut       = cms.double(0.06), # D0 mass window +- pdg value
   BMassCut        = cms.double(0.12), # Bu mass window +- pdg value
   savetrack       = cms.bool(savetrack),
   verbose         = cms.int32(0)
)

for entry in xgboost_models:
    BDh.xgboost_models.append(entry[0]),
    BDh.xgboost_variable_names.append(entry[1])

# Gen match

BDhGen = cms.EDProducer("BDhGen",
   genParticle = cms.InputTag('finalGenParticlesBPH'),
)
GenmatchTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BDhGen", "Genmatch"),
    cut       = cms.string(""),
    name      = cms.string("Genmatch"),
    doc       = cms.string("genpart Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        channelFlag            = uint('channelFlag'),
        b_numberOfDaughters    = uint('b_numberOfDaughters'),
        d0_numberOfDaughters   = uint('d0_numberOfDaughters'),
        ks0_numberOfDaughters  = uint('ks0_numberOfDaughters'),
        ks0_flight_distance    = ufloat('ks0_flight_distance'),
        ks0_flight_distance_2D = ufloat('ks0_flight_distance_2D'),
        idx_b                  = uint('idx_b'),
        b_charge               = uint('b_charge'),
        idx_b_pion             = uint('idx_b_pion'),
        idx_b_kaon             = uint('idx_b_kaon'),
        idx_dstar0             = uint('idx_dstar0'),
        idx_dstar_decay        = uint('idx_dstar_decay'),
        idx_d0                 = uint('idx_d0'),
        idx_ks                 = uint('idx_ks'),
        idx_kspip              = uint('idx_kspip'),
        idx_kspim              = uint('idx_kspim'),
        idx_pip                = uint('idx_pip'),
        idx_pim                = uint('idx_pim'),
        idx_pi0_1              = uint('idx_pi0_1'),
        idx_pi0_2              = uint('idx_pi0_2'),
        )
)
Countgenpart = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDhGen", "Genmatch")
)

if savetrack:
    # Tracks 
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
            dxy                    =    ufloat('dxy'),
            dz                     =    ufloat('dz'),
            dxySig                 =    ufloat('dxySig'),
            dzSig                  =    ufloat('dzSig'),
            bt_pt                  =    ufloat('bt_pt'),
            bt_ptErr               =    ufloat('ptErr'),
            trackHighPurity        =    ufloat('trackHighPurity'),
            dxy_bs                 =    ufloat('dxy_bs'),
            dz_bs                  =    ufloat('dz_bs'),
            dxySig_bs              =    ufloat('dxySig_bs'),
            dzSig_bs               =    ufloat('dzSig_bs'),
            dxy_pv                 =    ufloat('dxy_pv'),
            dz_pv                  =    ufloat('dz_pv'),
            dxySig_pv              =    ufloat('dxySig_pv'),
            dzSig_pv               =    ufloat('dzSig_pv'),
            normChi2               =    ufloat('normChi2'),
            nValidPixelHits        =    uint('nValidPixelHits'),
            nValidHits             =    uint('nValidHits')
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
    # Tracks
    DiTrackTable = cms.EDProducer(
        "SimpleCompositeCandidateFlatTableProducer",
        src  = cms.InputTag("BDh:SelectedDiTracks"),
        cut  = cms.string(""),
        name = cms.string("DiTrack"),
        doc  = cms.string("ditrack collection"),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(
            CandVars,
            leg1_idx           = Var("userInt('leg1_idx')", int, doc="leg1_idx"),
            leg2_idx           = Var("userInt('leg2_idx')", int, doc="leg2_idx"),
            cxPtR2             = ufloat('cxPtR2'),
            cxPtx              = ufloat('cxPtx'),
            cxPty              = ufloat('cxPty'),
            cxPtz              = ufloat('cxPtz'),
            dot                = ufloat('dot'),
            dca                = ufloat('dca'),
            massSquared        = ufloat('massSquared'),
            trk1_bs_dca        = ufloat('trk1_bs_dca'),
            trk2_bs_dca        = ufloat('trk2_bs_dca'),
            trk1_pv_dca        = ufloat('trk1_pv_dca'),
            trk2_pv_dca        = ufloat('trk2_pv_dca'),
            trk1_bs_dcaErr     = ufloat('trk1_bs_dcaErr'),
            trk2_bs_dcaErr     = ufloat('trk2_bs_dcaErr'),
            trk1_pv_dcaErr     = ufloat('trk1_pv_dcaErr'),
            trk2_pv_dcaErr     = ufloat('trk2_pv_dcaErr'),
            ),
    )


# B
BTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BDh", "B"),
    cut       = cms.string(""),
    name      = cms.string("B"),
    doc       = cms.string("B Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
        B_fitted_mass             = ufloat('B_fitted_mass'),
        DiTrack_idx1              = uint('DiTrack_idx1'),
        DiTrack_idx2              = uint('DiTrack_idx2'),
        Track_idx1                = uint('Track_idx1'),
        Track_idx2                = uint('Track_idx2'),
        Track_idx3                = uint('Track_idx3'),
        Track_idx4                = uint('Track_idx4'),
        DiTrk1_cxPtR2             = ufloat('DiTrk1_cxPtR2'),
        DiTrk1_cxPtz              = ufloat('DiTrk1_cxPtz'),
        DiTrk1_dot                = ufloat('DiTrk1_dot'),
        DiTrk1_dca                = ufloat('DiTrk1_dca'),
        DiTrk1_massSquared        = ufloat('DiTrk1_massSquared'),
        DiTrk1_trk1_bs_dca     = ufloat('DiTrk1_trk1_bs_dca'),
        DiTrk1_trk2_bs_dca     = ufloat('DiTrk1_trk2_bs_dca'),
        DiTrk1_trk1_pv_dca     = ufloat('DiTrk1_trk1_pv_dca'),
        DiTrk1_trk2_pv_dca     = ufloat('DiTrk1_trk2_pv_dca'),
        DiTrk1_trk1_bs_dcaSig     = ufloat('DiTrk1_trk1_bs_dcaSig'),
        DiTrk1_trk2_bs_dcaSig     = ufloat('DiTrk1_trk2_bs_dcaSig'),
        DiTrk1_trk1_pv_dcaSig     = ufloat('DiTrk1_trk1_pv_dcaSig'),
        DiTrk1_trk2_pv_dcaSig     = ufloat('DiTrk1_trk2_pv_dcaSig'),

        DiTrk2_cxPtR2             = ufloat('DiTrk2_cxPtR2'),
        DiTrk2_cxPtz              = ufloat('DiTrk2_cxPtz'),
        DiTrk2_dot                = ufloat('DiTrk2_dot'),
        DiTrk2_dca                = ufloat('DiTrk2_dca'),
        DiTrk2_massSquared        = ufloat('DiTrk2_massSquared'),
        DiTrk2_trk1_bs_dca     = ufloat('DiTrk2_trk1_bs_dca'),
        DiTrk2_trk2_bs_dca     = ufloat('DiTrk2_trk2_bs_dca'),
        DiTrk2_trk1_pv_dca     = ufloat('DiTrk2_trk1_pv_dca'),
        DiTrk2_trk2_pv_dca     = ufloat('DiTrk2_trk2_pv_dca'),
        DiTrk2_trk1_bs_dcaSig     = ufloat('DiTrk2_trk1_bs_dcaSig'),
        DiTrk2_trk2_bs_dcaSig     = ufloat('DiTrk2_trk2_bs_dcaSig'),
        DiTrk2_trk1_pv_dcaSig     = ufloat('DiTrk2_trk1_pv_dcaSig'),
        DiTrk2_trk2_pv_dcaSig     = ufloat('DiTrk2_trk2_pv_dcaSig'),

        Ks0_Kin_vtx_x       =    ufloat("Ks0_Kin_vtx_x"),
        Ks0_Kin_vtx_y       =    ufloat("Ks0_Kin_vtx_y"),
        Ks0_Kin_vtx_r       =    ufloat("Ks0_Kin_vtx_r"),
        Ks0_Kin_vtx_z       =    ufloat("Ks0_Kin_vtx_z"),
        Ks0_Kin_chi2        =    ufloat("Ks0_Kin_chi2"),
        Ks0_Kin_dof         =    ufloat("Ks0_Kin_dof"),
        Ks0_Kin_prob        =    ufloat("Ks0_Kin_prob"),
        Ks0_Kin_pt          =    ufloat("Ks0_Kin_pt"),
        Ks0_Kin_eta         =    ufloat("Ks0_Kin_eta"),
        Ks0_Kin_phi         =    ufloat("Ks0_Kin_phi"),
        Ks0_Kin_mass        =    ufloat("Ks0_Kin_mass"),
        Ks0_Kin_massErr     =    ufloat("Ks0_Kin_massErr"),
        Ks0_Kin_trk1_pt     =    ufloat("Ks0_Kin_trk1_pt"), 
        Ks0_Kin_trk1_eta    =    ufloat("Ks0_Kin_trk1_eta"),
        Ks0_Kin_trk1_phi    =    ufloat("Ks0_Kin_trk1_phi"),
        Ks0_Kin_trk2_pt     =    ufloat("Ks0_Kin_trk2_pt"),
        Ks0_Kin_trk2_eta    =    ufloat("Ks0_Kin_trk2_eta"),
        Ks0_Kin_trk2_phi    =    ufloat("Ks0_Kin_trk2_phi"),
        Ks0_Kin_bs_alpha_2D =  ufloat("Ks0_Kin_bs_alpha_2D"),
        Ks0_Kin_pv_alpha_2D = ufloat("Ks0_Kin_pv_alpha_2D"),
        Ks0_Kin_pv_alpha_3D = ufloat("Ks0_Kin_pv_alpha_3D"),
        Ks0_Kin_d0_alpha_2D = ufloat("Ks0_Kin_d0_alpha_2D"),
        Ks0_Kin_d0_alpha_3D = ufloat("Ks0_Kin_d0_alpha_3D"),
        Ks0_Kin_bs_l_xy     = ufloat("Ks0_Kin_bs_l_xy"),
        Ks0_Kin_bs_l_xySig  = ufloat("Ks0_Kin_bs_l_xySig"),
        Ks0_Kin_pv_l_xy     = ufloat("Ks0_Kin_pv_l_xy"),
        Ks0_Kin_pv_l_xySig  = ufloat("Ks0_Kin_pv_l_xySig"),
        Ks0_Kin_pv_l_xyz    = ufloat("Ks0_Kin_pv_l_xyz"),
        Ks0_Kin_pv_l_xyzSig = ufloat("Ks0_Kin_pv_l_xyzSig"),
        Ks0_Kin_d0_l_xy     = ufloat("Ks0_Kin_d0_l_xy"),
        Ks0_Kin_d0_l_xySig  = ufloat("Ks0_Kin_d0_l_xySig"),
        Ks0_Kin_d0_l_xyz    = ufloat("Ks0_Kin_d0_l_xyz"),
        Ks0_Kin_d0_l_xyzSig = ufloat("Ks0_Kin_d0_l_xyzSig"),
        Ks0_Kin_d0_dca      = ufloat("Ks0_Kin_d0_dca"),
        Ks0_Kin_d0_dcaSig   = ufloat("Ks0_Kin_d0_dcaSig"),

        D0_premass    = ufloat('D0_premass'),
        D0_Kin_vtx_x    =    ufloat("D0_Kin_vtx_x"),
        D0_Kin_vtx_y    =    ufloat("D0_Kin_vtx_y"),
        D0_Kin_vtx_r    =    ufloat("D0_Kin_vtx_r"),
        D0_Kin_vtx_z    =    ufloat("D0_Kin_vtx_z"),
        D0_Kin_chi2     =    ufloat("D0_Kin_chi2"),
        D0_Kin_dof      =    ufloat("D0_Kin_dof"),
        D0_Kin_prob     =    ufloat("D0_Kin_prob"),
        D0_Kin_pt       =    ufloat("D0_Kin_pt"),
        D0_Kin_eta      =    ufloat("D0_Kin_eta"),
        D0_Kin_phi      =    ufloat("D0_Kin_phi"),
        D0_Kin_mass     =    ufloat("D0_Kin_mass"),
        D0_Kin_massErr  =    ufloat("D0_Kin_massErr"),
        D0_Kin_trk3_pt   =    ufloat("D0_Kin_trk3_pt"),
        D0_Kin_trk3_eta  =    ufloat("D0_Kin_trk3_eta"),
        D0_Kin_trk3_phi  =    ufloat("D0_Kin_trk3_phi"),
        D0_Kin_trk4_pt   =    ufloat("D0_Kin_trk4_pt"),
        D0_Kin_trk4_eta  =    ufloat("D0_Kin_trk4_eta"),
        D0_Kin_trk4_phi  =    ufloat("D0_Kin_trk4_phi"),
        D0_Kin_ks0_pt   =    ufloat("D0_Kin_ks0_pt"),
        D0_Kin_ks0_eta  =    ufloat("D0_Kin_ks0_eta"),
        D0_Kin_ks0_phi  =    ufloat("D0_Kin_ks0_phi"),
        D0_Kin_bs_alpha_2D =  ufloat("D0_Kin_bs_alpha_2D"),
        D0_Kin_pv_alpha_2D = ufloat("D0_Kin_pv_alpha_2D"),
        D0_Kin_pv_alpha_3D = ufloat("D0_Kin_pv_alpha_3D"),
        D0_Kin_b_alpha_2D = ufloat("D0_Kin_b_alpha_2D"),
        D0_Kin_b_alpha_3D = ufloat("D0_Kin_b_alpha_3D"),
        D0_Kin_bs_l_xy     = ufloat("D0_Kin_bs_l_xy"),
        D0_Kin_bs_l_xySig  = ufloat("D0_Kin_bs_l_xySig"),
        D0_Kin_pv_l_xy     = ufloat("D0_Kin_pv_l_xy"),
        D0_Kin_pv_l_xySig  = ufloat("D0_Kin_pv_l_xySig"),
        D0_Kin_pv_l_xyz    = ufloat("D0_Kin_pv_l_xyz"),
        D0_Kin_pv_l_xyzSig = ufloat("D0_Kin_pv_l_xyzSig"),
        D0_Kin_b_l_xy     = ufloat("D0_Kin_b_l_xy"),
        D0_Kin_b_l_xySig  = ufloat("D0_Kin_b_l_xySig"),
        D0_Kin_b_l_xyz    = ufloat("D0_Kin_b_l_xyz"),
        D0_Kin_b_l_xyzSig = ufloat("D0_Kin_b_l_xyzSig"),
        D0_Kin_b_dca      = ufloat("D0_Kin_b_dca"),
        D0_Kin_b_dcaSig   = ufloat("D0_Kin_b_dcaSig"),
        
        BTrack_idx  = uint('BTrack_idx'),
        B_premass    = ufloat('B_premass'),
        B_Kin_vtx_x    =    ufloat("B_Kin_vtx_x"),
        B_Kin_vtx_y    =    ufloat("B_Kin_vtx_y"),
        B_Kin_vtx_r    =    ufloat("B_Kin_vtx_r"),
        B_Kin_vtx_z    =    ufloat("B_Kin_vtx_z"),
        B_Kin_chi2     =    ufloat("B_Kin_chi2"),
        B_Kin_dof      =    ufloat("B_Kin_dof"),
        B_Kin_prob     =    ufloat("B_Kin_prob"),
        B_Kin_pt       =    ufloat("B_Kin_pt"),
        B_Kin_eta      =    ufloat("B_Kin_eta"),
        B_Kin_phi      =    ufloat("B_Kin_phi"),
        B_Kin_mass     =    ufloat("B_Kin_mass"),
        B_Kin_massErr  =    ufloat("B_Kin_massErr"),
        B_Kin_trk_pt   =    ufloat("B_Kin_trk_pt"),
        B_Kin_trk_charge   =    ufloat("B_Kin_trk_charge"),
        B_Kin_trk_eta  =    ufloat("B_Kin_trk_eta"),
        B_Kin_trk_phi  =    ufloat("B_Kin_trk_phi"),
        B_Kin_D0_pt   =    ufloat("B_Kin_D0_pt"),
        B_Kin_D0_eta  =    ufloat("B_Kin_D0_eta"),
        B_Kin_D0_phi  =    ufloat("B_Kin_D0_phi"),
        B_Kin_bs_alpha_2D =  ufloat("B_Kin_bs_alpha_2D"),
        B_Kin_pv_alpha_2D = ufloat("B_Kin_pv_alpha_2D"),
        B_Kin_pv_alpha_3D = ufloat("B_Kin_pv_alpha_3D"),
        B_Kin_bs_l_xy     = ufloat("B_Kin_bs_l_xy"),
        B_Kin_bs_l_xySig  = ufloat("B_Kin_bs_l_xySig"),
        B_Kin_pv_l_xy     = ufloat("B_Kin_pv_l_xy"),
        B_Kin_pv_l_xySig  = ufloat("B_Kin_pv_l_xySig"),
        B_Kin_pv_l_xyz    = ufloat("B_Kin_pv_l_xyz"),
        B_Kin_pv_l_xyzSig = ufloat("B_Kin_pv_l_xyzSig"),
        B_Kin_trk_bs_dca      = ufloat("B_Kin_trk_bs_dca"),
        B_Kin_trk_bs_dcaSig   = ufloat("B_Kin_trk_bs_dcaSig"),
        B_Kin_trk_pv_dca      = ufloat("B_Kin_trk_pv_dca"),
        B_Kin_trk_pv_dcaSig   = ufloat("B_Kin_trk_pv_dcaSig"),
        B_Kin_trk_b_dca      = ufloat("B_Kin_trk_b_dca"),
        B_Kin_trk_b_dcaSig   = ufloat("B_Kin_trk_b_dcaSig"),

        B_kaon_premass          =    ufloat('B_kaon_premass'),
        B_kaon_Kin_vtx_x        =    ufloat("B_kaon_Kin_vtx_x"),
        B_kaon_Kin_vtx_y        =    ufloat("B_kaon_Kin_vtx_y"),
        B_kaon_Kin_vtx_r        =    ufloat("B_kaon_Kin_vtx_r"),
        B_kaon_Kin_vtx_z        =    ufloat("B_kaon_Kin_vtx_z"),
        B_kaon_Kin_chi2         =    ufloat("B_kaon_Kin_chi2"),
        B_kaon_Kin_dof          =    ufloat("B_kaon_Kin_dof"),
        B_kaon_Kin_prob         =    ufloat("B_kaon_Kin_prob"),
        B_kaon_Kin_pt           =    ufloat("B_kaon_Kin_pt"),
        B_kaon_Kin_eta          =    ufloat("B_kaon_Kin_eta"),
        B_kaon_Kin_phi          =    ufloat("B_kaon_Kin_phi"),
        B_kaon_Kin_mass         =    ufloat("B_kaon_Kin_mass"),
        B_kaon_Kin_massErr      =    ufloat("B_kaon_Kin_massErr"),
        B_kaon_Kin_trk_pt       =    ufloat("B_kaon_Kin_trk_pt"),
        B_kaon_Kin_trk_charge   =    ufloat("B_kaon_Kin_trk_charge"),
        B_kaon_Kin_trk_eta      =    ufloat("B_kaon_Kin_trk_eta"),
        B_kaon_Kin_trk_phi      =    ufloat("B_kaon_Kin_trk_phi"),
        B_kaon_Kin_D0_pt        =    ufloat("B_kaon_Kin_D0_pt"),
        B_kaon_Kin_D0_eta       =    ufloat("B_kaon_Kin_D0_eta"),
        B_kaon_Kin_D0_phi       =    ufloat("B_kaon_Kin_D0_phi"),
        B_kaon_Kin_bs_alpha_2D  =    ufloat("B_kaon_Kin_bs_alpha_2D"),
        B_kaon_Kin_pv_alpha_2D  =    ufloat("B_kaon_Kin_pv_alpha_2D"),
        B_kaon_Kin_pv_alpha_3D  =    ufloat("B_kaon_Kin_pv_alpha_3D"),
        B_kaon_Kin_bs_l_xy      =    ufloat("B_kaon_Kin_bs_l_xy"),
        B_kaon_Kin_bs_l_xySig   =    ufloat("B_kaon_Kin_bs_l_xySig"),
        B_kaon_Kin_pv_l_xy      =    ufloat("B_kaon_Kin_pv_l_xy"),
        B_kaon_Kin_pv_l_xySig   =    ufloat("B_kaon_Kin_pv_l_xySig"),
        B_kaon_Kin_pv_l_xyz     =    ufloat("B_kaon_Kin_pv_l_xyz"),
        B_kaon_Kin_pv_l_xyzSig  =    ufloat("B_kaon_Kin_pv_l_xyzSig"),
        )
)

for entry in xgboost_models:
    setattr(BTable.variables,
            "xgb_%s" % entry[1],
            Var("userFloat('xgb_%s')" % entry[1], float, doc = "New XGBoost MVA id")
            )

CountB = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDh", "B")
)

BMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = BTable.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(521),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),                            # Minimum deltaR for the match
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

if savetrack:
    BDhSequence = cms.Sequence(BDh)
    BDhSequenceTable = cms.Sequence(PionTrackTable + DiTrackTable + BTable )
    BDhSequenceMC = cms.Sequence(BDh + PionTrackMCMatch + BMCMatch + BDhGen)
    BDhSequenceMCTable = cms.Sequence(PionTrackTable + DiTrackTable + PionTrackMCTable + BTable + BMCTable + GenmatchTable)
else:
    BDhSequence = cms.Sequence(BDh )
    BDhSequenceTable = cms.Sequence(BTable )
    BDhSequenceMC = cms.Sequence(BDh + BMCMatch + BDhGen)
    BDhSequenceMCTable = cms.Sequence(BTable + BMCTable + GenmatchTable)



