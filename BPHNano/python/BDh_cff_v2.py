import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *


savetrack=False

BDh = cms.EDProducer("BDhFitter_v2",    
   # which beamSpot to reference
   beamSpot = cms.InputTag('offlineBeamSpot'),
   vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
   tracks          = cms.InputTag("packedPFCandidates"),
   lostTracks      = cms.InputTag("lostTracks"),
   genParticle = cms.InputTag('finalGenParticlesBPH'),
   # Tracks
   tkNHitsCut = cms.int32(3), # Number of valid hits on track
   tkPtCut = cms.double(0.3), # Pt of track
   tkEtaCut = cms.double(2.4), # Eta of track
   tkChi2Cut = cms.double(30.), # Track normalized Chi2
   # diTracks
   vtxChi2Cut = cms.double(20), # Vertex KLM chi2
   vtxDecaySigXYCut = cms.double(3), # KLM XY decay distance significance
   vtxDecaySigXYZCut = cms.double(0.1), # KLM XYZ decay distance significance
   cosThetaXYCut = cms.double(0.9),  # cos(angleXY) between x and p of V0 candidate
   cosThetaXYZCut = cms.double(-1), # cos(angleXYZ) between x and p of V0 candidate
   # reco ks0
   mPiPiCut = cms.double(0.7), # invariant mass of track pair, assuming both tracks are charged pions
   kShortMassCut = cms.double(0.1), # V0 mass window +- pdg value
   D0MassCut = cms.double(0.2), # V0 mass window +- pdg value
   BMassCut = cms.double(0.3), # V0 mass window +- pdg value
   savetrack = cms.bool(savetrack),
   verbose = cms.int32(0)
)

# Gen match
GenmatchTable = cms.EDProducer(
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
        ks0_flight_distance= ufloat('ks0_flight_distance'),
        ks0_flight_distance_2D= ufloat('ks0_flight_distance_2D'),
        idx_b= uint('idx_b'),
        b_charge= uint('b_charge'),
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
        DiTrack_idx1  = uint('DiTrack_idx1'),
        DiTrack_idx2  = uint('DiTrack_idx2'),
        Track_idx1  = uint('Track_idx1'),
        Track_idx2  = uint('Track_idx2'),
        Track_idx3  = uint('Track_idx3'),
        Track_idx4  = uint('Track_idx4'),
        BTrack_idx  = uint('BTrack_idx'),
        DiTrk1_cxPtR2     =    ufloat('DiTrk1_cxPtR2'),
        DiTrk1_cxPtz      =    ufloat('DiTrk1_cxPtz'),
        DiTrk1_dot        =    ufloat('DiTrk1_dot'),
        DiTrk1_trk1_dcaSig=    ufloat('DiTrk1_trk1_dcaSig'),
        DiTrk1_trk2_dcaSig=    ufloat('DiTrk1_trk2_dcaSig'),
        DiTrk1_dca        =    ufloat('DiTrk1_dca'),
        DiTrk1_massSquared=    ufloat('DiTrk1_massSquared'),
        DiTrk1_KLM_vtx_x  =    ufloat('DiTrk1_KLM_vtx_x'),
        DiTrk1_KLM_vtx_y  =    ufloat('DiTrk1_KLM_vtx_y'),
        DiTrk1_KLM_vtx_z  =    ufloat('DiTrk1_KLM_vtx_z'),
        DiTrk1_KLM_chi2           =      ufloat("DiTrk1_KLM_chi2"),
        DiTrk1_KLM_ndof           =      ufloat("DiTrk1_KLM_ndof"),
        DiTrk1_KLM_normalizedChi2 =      ufloat("DiTrk1_KLM_normalizedChi2"),
        DiTrk1_KLM_distMagXY      =      ufloat("DiTrk1_KLM_distMagXY"),
        DiTrk1_KLM_distMagXYErr   =      ufloat("DiTrk1_KLM_distMagXYErr"),
        DiTrk1_KLM_cos_theta_XY   =      ufloat("DiTrk1_KLM_cos_theta_XY"),
        DiTrk2_cxPtR2     =    ufloat('DiTrk2_cxPtR2'),
        DiTrk2_cxPtz      =    ufloat('DiTrk2_cxPtz'),
        DiTrk2_dot        =    ufloat('DiTrk2_dot'),
        DiTrk2_trk1_dcaSig=    ufloat('DiTrk2_trk1_dcaSig'),
        DiTrk2_trk2_dcaSig=    ufloat('DiTrk2_trk2_dcaSig'),
        DiTrk2_dca        =    ufloat('DiTrk2_dca'),
        DiTrk2_massSquared=    ufloat('DiTrk2_massSquared'),
        DiTrk2_KLM_vtx_x  =    ufloat('DiTrk2_KLM_vtx_x'),
        DiTrk2_KLM_vtx_y  =    ufloat('DiTrk2_KLM_vtx_y'),
        DiTrk2_KLM_vtx_z  =    ufloat('DiTrk2_KLM_vtx_z'),
        DiTrk2_KLM_chi2           =      ufloat("DiTrk2_KLM_chi2"),
        DiTrk2_KLM_ndof           =      ufloat("DiTrk2_KLM_ndof"),
        DiTrk2_KLM_normalizedChi2 =      ufloat("DiTrk2_KLM_normalizedChi2"),
        DiTrk2_KLM_distMagXY      =      ufloat("DiTrk2_KLM_distMagXY"),
        DiTrk2_KLM_distMagXYErr   =      ufloat("DiTrk2_KLM_distMagXYErr"),
        DiTrk2_KLM_cos_theta_XY   =      ufloat("DiTrk2_KLM_cos_theta_XY"),
        Ks0_Kin_vtx_x    =    ufloat("Ks0_Kin_vtx_x"),
        Ks0_Kin_vtx_y    =    ufloat("Ks0_Kin_vtx_y"),
        Ks0_Kin_vtx_z    =    ufloat("Ks0_Kin_vtx_z"),
        Ks0_Kin_chi2     =    ufloat("Ks0_Kin_chi2"),
        Ks0_Kin_dof      =    ufloat("Ks0_Kin_dof"),
        Ks0_Kin_prob     =    ufloat("Ks0_Kin_prob"),
        Ks0_Kin_pt       =    ufloat("Ks0_Kin_pt"),
        Ks0_Kin_eta      =    ufloat("Ks0_Kin_eta"),
        Ks0_Kin_phi      =    ufloat("Ks0_Kin_phi"),
        Ks0_Kin_mass     =    ufloat("Ks0_Kin_mass"),
        Ks0_Kin_massErr  =    ufloat("Ks0_Kin_massErr"),
        Ks0_Kin_fitted_cos_theta_2D  =    ufloat("Ks0_Kin_fitted_cos_theta_2D"),
        Ks0_Kin_l_xy      =    ufloat("Ks0_Kin_l_xy"),
        Ks0_Kin_l_xy_unc  =    ufloat("Ks0_Kin_l_xy_unc"),
        Ks0_Kin_trk1_pt   =    ufloat("Ks0_Kin_trk1_pt"), 
        Ks0_Kin_trk1_eta  =    ufloat("Ks0_Kin_trk1_eta"),
        Ks0_Kin_trk1_phi  =    ufloat("Ks0_Kin_trk1_phi"),
        Ks0_Kin_trk2_pt   =    ufloat("Ks0_Kin_trk2_pt"),
        Ks0_Kin_trk2_eta  =    ufloat("Ks0_Kin_trk2_eta"),
        Ks0_Kin_trk2_phi  =    ufloat("Ks0_Kin_trk2_phi"),
        D0_premass    = ufloat('D0_premass'),
        D0_Kin_vtx_x    =    ufloat("D0_Kin_vtx_x"),
        D0_Kin_vtx_y    =    ufloat("D0_Kin_vtx_y"),
        D0_Kin_vtx_z    =    ufloat("D0_Kin_vtx_z"),
        D0_Kin_chi2     =    ufloat("D0_Kin_chi2"),
        D0_Kin_dof      =    ufloat("D0_Kin_dof"),
        D0_Kin_prob     =    ufloat("D0_Kin_prob"),
        D0_Kin_pt       =    ufloat("D0_Kin_pt"),
        D0_Kin_eta      =    ufloat("D0_Kin_eta"),
        D0_Kin_phi      =    ufloat("D0_Kin_phi"),
        D0_Kin_mass     =    ufloat("D0_Kin_mass"),
        D0_Kin_massErr  =    ufloat("D0_Kin_massErr"),
        D0_Kin_fitted_cos_theta_2D  =    ufloat("D0_Kin_fitted_cos_theta_2D"),
        D0_Kin_l_xy      =    ufloat("D0_Kin_l_xy"),
        D0_Kin_l_xy_unc  =    ufloat("D0_Kin_l_xy_unc"),
        D0_Kin_trk3_pt   =    ufloat("D0_Kin_trk3_pt"),
        D0_Kin_trk3_eta  =    ufloat("D0_Kin_trk3_eta"),
        D0_Kin_trk3_phi  =    ufloat("D0_Kin_trk3_phi"),
        D0_Kin_trk4_pt   =    ufloat("D0_Kin_trk4_pt"),
        D0_Kin_trk4_eta  =    ufloat("D0_Kin_trk4_eta"),
        D0_Kin_trk4_phi  =    ufloat("D0_Kin_trk4_phi"),
        D0_Kin_ks0_pt   =    ufloat("D0_Kin_ks0_pt"),
        D0_Kin_ks0_eta  =    ufloat("D0_Kin_ks0_eta"),
        D0_Kin_ks0_phi  =    ufloat("D0_Kin_ks0_phi"),
        ks0_flight_distance = ufloat("ks0_flight_distance"),
        ks0_flight_distance_2D = ufloat("ks0_flight_distance_2D"),
        B_premass    = ufloat('B_premass'),
        B_Kin_vtx_x    =    ufloat("B_Kin_vtx_x"),
        B_Kin_vtx_y    =    ufloat("B_Kin_vtx_y"),
        B_Kin_vtx_z    =    ufloat("B_Kin_vtx_z"),
        B_Kin_chi2     =    ufloat("B_Kin_chi2"),
        B_Kin_dof      =    ufloat("B_Kin_dof"),
        B_Kin_prob     =    ufloat("B_Kin_prob"),
        B_Kin_pt       =    ufloat("B_Kin_pt"),
        B_Kin_eta      =    ufloat("B_Kin_eta"),
        B_Kin_phi      =    ufloat("B_Kin_phi"),
        B_Kin_mass     =    ufloat("B_Kin_mass"),
        B_Kin_massErr  =    ufloat("B_Kin_massErr"),
        B_Kin_fitted_cos_theta_2D  =    ufloat("B_Kin_fitted_cos_theta_2D"),
        B_Kin_l_xy      =    ufloat("B_Kin_l_xy"),
        B_Kin_l_xy_unc  =    ufloat("B_Kin_l_xy_unc"),
        B_Kin_trk_pt   =    ufloat("B_Kin_trk_pt"),
        B_Kin_trk_charge   =    ufloat("B_Kin_trk_charge"),
        B_Kin_trk_eta  =    ufloat("B_Kin_trk_eta"),
        B_Kin_trk_phi  =    ufloat("B_Kin_trk_phi"),
        B_Kin_D0_pt   =    ufloat("B_Kin_D0_pt"),
        B_Kin_D0_eta  =    ufloat("B_Kin_D0_eta"),
        B_Kin_D0_phi  =    ufloat("B_Kin_D0_phi"),
        D0_flight_distance = ufloat("D0_flight_distance"),
        D0_flight_distance_2D = ufloat("D0_flight_distance_2D")
        )
)

CountB = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BDh", "B")
)

#D0MCMatch = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
#    src         = D0Table.src,                      # final reco collection
#    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
#    mcPdgId     = cms.vint32(421),                            # one or more PDG ID (13 = mu); absolute values (see below)
#    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
#    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
#    maxDeltaR   = cms.double(0.3),                            # Minimum deltaR for the match
#    maxDPtRel   = cms.double(1.0),                            # Minimum deltaPt/Pt for the match
#    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
#    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
#)
#
#
#D0MCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
#    recoObjects = D0Table.src,
#    genParts = cms.InputTag("finalGenParticlesBPH"),
#    mcMap = cms.InputTag("D0MCMatch"),
#    objName = D0Table.name,
#    objType = cms.string("Other"),
#    objBranchName = cms.string("genPart"),
#    genBranchName = cms.string("B"),
#    docString = cms.string("MC matching to status==2 D0"),
#)

if savetrack:
    BDhSequence = cms.Sequence(BDh)
    BDhSequenceTable = cms.Sequence(PionTrackTable + DiTrackTable + BTable )
    BDhSequenceMC = cms.Sequence(BDh + PionTrackMCMatch)
    BDhSequenceMCTable = cms.Sequence(PionTrackTable + GenmatchTable + DiTrackTable + PionTrackMCTable + BTable )
else:
    BDhSequence = cms.Sequence(BDh )
    BDhSequenceTable = cms.Sequence(BTable )
    BDhSequenceMC = cms.Sequence(BDh )
    BDhSequenceMCTable = cms.Sequence(GenmatchTable + BTable )




