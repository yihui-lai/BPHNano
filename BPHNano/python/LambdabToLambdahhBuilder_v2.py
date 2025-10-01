import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *


xgboost_models = [
]

savetrack=True

# Used to study lambdaB -> lambda0(->\p\pi) + 2h
LambdabToLambdahhv2 = cms.EDProducer("LambdabToLambdahhBuilder_v2",    
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
   tkIPSigXYCut    = cms.double(-1), # Track IP significance for all tracks
   # diTrack 1
   TrkSigXYCut       = cms.double(3),                 # track DCA significance
   vtxChi2Cut        = cms.double(6.63),              # Vertex KLM chi2
   vtxDecaySigXYCut  = cms.double(-1),                # (Not used) KLM XY decay distance significance
   vtxDecaySigXYZCut = cms.double(-1),                # (Not used) KLM XYZ decay distance significance
   cosThetaXYCut     = cms.double(0.9),             # cos(theta3D)  between x(ldb0->ld0) and p(ld0) candidate, a loose cut
   cosThetaXYZCut    = cms.double(-1),                # cos(angleXYZ) between x and p of V0 candidate
   # diTrack 2
   Lambdab_tkPtCut          = cms.double(0.6),   # Pt cut of track 3, 4
   diTrack2_dca      = cms.double(1),    # 0.2cm is a decent cut to remove comb.bkg
   Trk34SigXYCut     = cms.double(-1),    # track DCA significance
   # reco ld0
   ld0_l_xyzSigCut = cms.double(3), # ld0 flight distancesignificance from ldb0
   # ldb0
   ldb0_PtCut            = cms.double(5), # ldb0 Pt cut, GeV 
   ldb0_vtxDecaySigXYCut  = cms.double(-1), # ldb0 XY distance significance from PV
   # mass
   ld0_MassCut   = cms.double(0.01), # ld0 mass window +- pdg value
   ldb0_MassCut  = cms.double(0.15), # ldb0 mass window +- pdg value
   savetrack       = cms.bool(savetrack),
   verbose         = cms.int32(0)
)

for entry in xgboost_models:
    LambdabToLambdahhv2.xgboost_models.append(entry[0]),
    LambdabToLambdahhv2.xgboost_variable_names.append(entry[1])

if savetrack:
    # Tracks 
    PionTrackTable = cms.EDProducer(
        "SimpleCompositeCandidateFlatTableProducer",
        src  = cms.InputTag("LambdabToLambdahhv2:SelectedTracks"),
        cut  = cms.string(""),
        name = cms.string("Track"),
        doc  = cms.string("track collection"),
        singleton = cms.bool(False),
        extension = cms.bool(False),
        variables = cms.PSet(
            CandVars,
            isPacked  = Var("userInt('isPacked')", bool, doc="track from packedCandidate collection"),
            isLostTrk = Var("userInt('isLostTrk')", bool, doc="track from lostTrack collection"),
            DCASig  = Var("userFloat('DCASig')", float, doc="significance of xy-distance of closest approach wrt beamspot"),
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
        mcPdgId     = cms.vint32(321),                     # one or more PDG ID (321 = charged kaon, 211 = charged pion); absolute values (see below)
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
        src  = cms.InputTag("LambdabToLambdahhv2:SelectedDiTracks"),
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


# Ldb0
LambdabTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("LambdabToLambdahhv2", "Ldb0"),
    cut       = cms.string(""),
    name      = cms.string("Ldb0"),
    doc       = cms.string("Ldb0 Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        CandVars,
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
        DiTrk1_trk1_mass          = ufloat('DiTrk1_trk1_mass'),
        DiTrk1_trk2_mass          = ufloat('DiTrk1_trk2_mass'),
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

        Ld0_Kin_vtx_x       =    ufloat("Ld0_Kin_vtx_x"),
        Ld0_Kin_vtx_y       =    ufloat("Ld0_Kin_vtx_y"),
        Ld0_Kin_vtx_r       =    ufloat("Ld0_Kin_vtx_r"),
        Ld0_Kin_vtx_z       =    ufloat("Ld0_Kin_vtx_z"),
        Ld0_Kin_chi2        =    ufloat("Ld0_Kin_chi2"),
        Ld0_Kin_dof         =    ufloat("Ld0_Kin_dof"),
        Ld0_Kin_prob        =    ufloat("Ld0_Kin_prob"),
        Ld0_Kin_pt          =    ufloat("Ld0_Kin_pt"),
        Ld0_Kin_eta         =    ufloat("Ld0_Kin_eta"),
        Ld0_Kin_phi         =    ufloat("Ld0_Kin_phi"),
        Ld0_Kin_mass        =    ufloat("Ld0_Kin_mass"),
        Ld0_Kin_mass_hyp    =    uint("Ld0_Kin_mass_hyp"),
        Ld0_Kin_massErr     =    ufloat("Ld0_Kin_massErr"),
        Ld0_Kin_trk1_pt     =    ufloat("Ld0_Kin_trk1_pt"), 
        Ld0_Kin_trk1_eta    =    ufloat("Ld0_Kin_trk1_eta"),
        Ld0_Kin_trk1_phi    =    ufloat("Ld0_Kin_trk1_phi"),
        Ld0_Kin_trk2_pt     =    ufloat("Ld0_Kin_trk2_pt"),
        Ld0_Kin_trk2_eta    =    ufloat("Ld0_Kin_trk2_eta"),
        Ld0_Kin_trk2_phi    =    ufloat("Ld0_Kin_trk2_phi"),
        Ld0_Kin_bs_alpha_2D =  ufloat("Ld0_Kin_bs_alpha_2D"),
        Ld0_Kin_pv_alpha_2D = ufloat("Ld0_Kin_pv_alpha_2D"),
        Ld0_Kin_pv_alpha_3D = ufloat("Ld0_Kin_pv_alpha_3D"),
        Ld0_Kin_ldb0_alpha_2D = ufloat("Ld0_Kin_ldb0_alpha_2D"),
        Ld0_Kin_ldb0_alpha_3D = ufloat("Ld0_Kin_ldb0_alpha_3D"),
        Ld0_Kin_bs_l_xy     = ufloat("Ld0_Kin_bs_l_xy"),
        Ld0_Kin_bs_l_xySig  = ufloat("Ld0_Kin_bs_l_xySig"),
        Ld0_Kin_pv_l_xy     = ufloat("Ld0_Kin_pv_l_xy"),
        Ld0_Kin_pv_l_xySig  = ufloat("Ld0_Kin_pv_l_xySig"),
        Ld0_Kin_pv_l_xyz    = ufloat("Ld0_Kin_pv_l_xyz"),
        Ld0_Kin_pv_l_xyzSig = ufloat("Ld0_Kin_pv_l_xyzSig"),
        Ld0_Kin_ldb0_l_xy     = ufloat("Ld0_Kin_ldb0_l_xy"),
        Ld0_Kin_ldb0_l_xySig  = ufloat("Ld0_Kin_ldb0_l_xySig"),
        Ld0_Kin_ldb0_l_xyz    = ufloat("Ld0_Kin_ldb0_l_xyz"),
        Ld0_Kin_ldb0_l_xyzSig = ufloat("Ld0_Kin_ldb0_l_xyzSig"),
        Ld0_Kin_ldb0_dca      = ufloat("Ld0_Kin_ldb0_dca"),
        Ld0_Kin_ldb0_dcaSig   = ufloat("Ld0_Kin_ldb0_dcaSig"),

        Ldb0_premass      = ufloat('Ldb0_premass'),
        Ldb0_Kin_vtx_x    =    ufloat("Ldb0_Kin_vtx_x"),
        Ldb0_Kin_vtx_y    =    ufloat("Ldb0_Kin_vtx_y"),
        Ldb0_Kin_vtx_r    =    ufloat("Ldb0_Kin_vtx_r"),
        Ldb0_Kin_vtx_z    =    ufloat("Ldb0_Kin_vtx_z"),
        Ldb0_Kin_chi2     =    ufloat("Ldb0_Kin_chi2"),
        Ldb0_Kin_dof      =    ufloat("Ldb0_Kin_dof"),
        Ldb0_Kin_prob     =    ufloat("Ldb0_Kin_prob"),
        Ldb0_Kin_pt       =    ufloat("Ldb0_Kin_pt"),
        Ldb0_Kin_eta      =    ufloat("Ldb0_Kin_eta"),
        Ldb0_Kin_phi      =    ufloat("Ldb0_Kin_phi"),
        Ldb0_Kin_mass     =    ufloat("Ldb0_Kin_mass"),
        Ldb0_Kin_massErr  =    ufloat("Ldb0_Kin_massErr"),
        Ldb0_Kin_mass_fixlb0  =    ufloat("Ldb0_Kin_mass_fixlb0"),
        Ldb0_Kin_massErr_fixlb0  =    ufloat("Ldb0_Kin_massErr_fixlb0"),
        Ldb0_Kin_trk3_pt   =    ufloat("Ldb0_Kin_trk3_pt"),
        Ldb0_Kin_trk3_eta  =    ufloat("Ldb0_Kin_trk3_eta"),
        Ldb0_Kin_trk3_phi  =    ufloat("Ldb0_Kin_trk3_phi"),
        Ldb0_Kin_trk4_pt   =    ufloat("Ldb0_Kin_trk4_pt"),
        Ldb0_Kin_trk4_eta  =    ufloat("Ldb0_Kin_trk4_eta"),
        Ldb0_Kin_trk4_phi  =    ufloat("Ldb0_Kin_trk4_phi"),
        Ldb0_Kin_ld0_pt   =    ufloat("Ldb0_Kin_ld0_pt"),
        Ldb0_Kin_ld0_eta  =    ufloat("Ldb0_Kin_ld0_eta"),
        Ldb0_Kin_ld0_phi  =    ufloat("Ldb0_Kin_ld0_phi"),

        Ldb0_Kin_trk3_iso04  =    ufloat("Ldb0_Kin_trk3_iso04"),
        Ldb0_Kin_trk4_iso04  =    ufloat("Ldb0_Kin_trk4_iso04"),
        Ldb0_Kin_ld0_iso04   =    ufloat("Ldb0_Kin_ld0_iso04"),
        Ldb0_Kin_bs_alpha_2D =  ufloat("Ldb0_Kin_bs_alpha_2D"),
        Ldb0_Kin_pv_alpha_2D = ufloat("Ldb0_Kin_pv_alpha_2D"),
        Ldb0_Kin_pv_alpha_3D = ufloat("Ldb0_Kin_pv_alpha_3D"),
        Ldb0_Kin_bs_l_xy     = ufloat("Ldb0_Kin_bs_l_xy"),
        Ldb0_Kin_bs_l_xySig  = ufloat("Ldb0_Kin_bs_l_xySig"),
        Ldb0_Kin_pv_l_xy     = ufloat("Ldb0_Kin_pv_l_xy"),
        Ldb0_Kin_pv_l_xySig  = ufloat("Ldb0_Kin_pv_l_xySig"),
        Ldb0_Kin_pv_l_xyz    = ufloat("Ldb0_Kin_pv_l_xyz"),
        Ldb0_Kin_pv_l_xyzSig = ufloat("Ldb0_Kin_pv_l_xyzSig"),
        )
)

for entry in xgboost_models:
    setattr(LambdabTable.variables,
            "xgb_%s" % entry[1],
            Var("userFloat('xgb_%s')" % entry[1], float, doc = "New XGBoost MVA id")
            )

LambdabMCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = LambdabTable.src,                      # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(5122),                            # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                            # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)


LambdabMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = LambdabTable.src,
    genParts = cms.InputTag("finalGenParticlesBPH"),
    mcMap = cms.InputTag("LambdabMCMatch"),
    objName = LambdabTable.name,
    objType = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("Ldb0"),
    docString = cms.string("MC matching to status==2 Ldb0"),
)



# Gen match
Lambdab0Gen = cms.EDProducer("Lambdab0Gen",
   genParticle = cms.InputTag('finalGenParticlesBPH'),
)
LambdaBGenmatchTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("Lambdab0Gen", "LambdaBGenmatch"),
    cut       = cms.string(""),
    name      = cms.string("LambdaBGenmatch"),
    doc       = cms.string("genpart Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        idx_lambdaB    = uint('idx_lambdaB'),
        lambdaB_charge = uint('lambdaB_charge'),
        channelFlag            = uint('channelFlag'),
        lambdaB_numberOfDaughters = uint('lambdaB_numberOfDaughters'),
        lambda0_numberOfDaughters = uint('lambda0_numberOfDaughters'),
        idx_lambda0 = uint("idx_lambda0"),
        idx_b_kaon1 = uint("idx_b_kaon1"),
        idx_b_kaon2 = uint("idx_b_kaon2"),
        idx_b_pion1 = uint("idx_b_pion1"),
        idx_b_pion2 = uint("idx_b_pion2"),
        idx_proton  = uint("idx_proton"),
        idx_pion_from_lambda0 = uint("idx_pion_from_lambda0")
        )
)
Countgenpart = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("Lambdab0Gen", "LambdaBGenmatch")
)

if savetrack:
    LambdabToLambdahhv2Sequence = cms.Sequence(LambdabToLambdahhv2)
    LambdabToLambdahhv2SequenceTable = cms.Sequence(PionTrackTable + DiTrackTable + LambdabTable )
    LambdabToLambdahhv2SequenceMC = cms.Sequence(LambdabToLambdahhv2 + LambdabMCMatch + PionTrackMCMatch + Lambdab0Gen)
    LambdabToLambdahhv2SequenceMCTable = cms.Sequence(PionTrackTable + DiTrackTable + LambdabTable + LambdabMCTable + PionTrackMCTable + LambdaBGenmatchTable)
else:
    LambdabToLambdahhv2Sequence = cms.Sequence(LambdabToLambdahhv2 )
    LambdabToLambdahhv2SequenceTable = cms.Sequence(LambdabTable )
    LambdabToLambdahhv2SequenceMC = cms.Sequence(LambdabToLambdahhv2 + LambdabMCMatch + Lambdab0Gen)
    LambdabToLambdahhv2SequenceMCTable = cms.Sequence(LambdabTable + LambdabMCTable + LambdaBGenmatchTable)



