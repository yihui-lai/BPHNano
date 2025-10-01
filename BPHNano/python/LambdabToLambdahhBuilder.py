import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *

########################### LambdaB0 -> Lambda0 hh ###########################

LambdabToLambdahh = cms.EDProducer(
    'LambdabToLambdahhBuilder',
    lambda0 = cms.InputTag('LambdaToPPi', 'SelectedLambdaCollection'),
    v0TransientTracks = cms.InputTag('LambdaToPPi', 'SelectedLambda'),
    tracks = cms.InputTag('tracksBPH', 'SelectedTracks'),
    transientTracks = cms.InputTag('tracksBPH', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection  = cms.string('pt > 6'
                                  '&& 5.0 < mass && mass < 6.8'),
    postVtxSelection = cms.string('5.3 < userFloat("fitted_mass") && userFloat("fitted_mass") < 6.5'
                                  '&& userFloat("sv_prob") > 0.001'
                                  '&& userFloat("fitted_pt") > 6'
                                  '&& abs(userFloat("fitted_cos_theta_2D")) > 0.7'),  # 94% eff.
)

########################### Tables ###########################

LambdabToLambdahhTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("LambdabToLambdahh"),
    cut       = cms.string(""),
    name      = cms.string("LambdabToLambdahh"),
    doc       = cms.string("LambdabToLambdahh Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        # pre-fit quantities
        CandVars,
        trk1_idx    = uint('trk1_idx'),
        trk2_idx    = uint('trk2_idx'),
        lambda_idx   = uint('lambda_idx'),
        trk3_idx    = uint('trk3_idx'),
        trk4_idx    = uint('trk4_idx'),
        trk1_mass   = ufloat('trk1_mass'),
        trk2_mass   = ufloat('trk2_mass'),
        trk12_dr            = ufloat('trk12_dr'),
        trk34_dr            = ufloat('trk34_dr'),
        trk34_dca           = ufloat('trk34_dca'),
        trk34_cxPtR         = ufloat('trk34_cxPtR'),
        trk1234_min_dr      = ufloat('trk1234_min_dr'),
        trk1234_max_dr      = ufloat('trk1234_max_dr'),
        # vtx info
        chi2      = ufloat('sv_chi2'),
        svprob    = ufloat('sv_prob'),
        cos2D     = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        l_xy      = ufloat('l_xy'),
        l_xy_unc  = ufloat('l_xy_unc'),
        # lambda0 flys
        ld0_ldb_fit_cos2D = ufloat('ld0_ldb_fitted_cos_theta_2D'),
        ld0_ldb_l_xy      = ufloat('ld0_ldb_l_xy'),
        ld0_ldb_l_xy_unc  = ufloat('ld0_ldb_l_xy_unc'),
        # post-fit momentum /masses
        mhh_fullfit    = ufloat('fitted_mhh'),
        mlambda_h1_fullfit = ufloat('fitted_mlambda_h1'),
        mlambda_h2_fullfit = ufloat('fitted_mlambda_h2'),
        fit_mass       = ufloat('fitted_mass'),
        fit_massErr    = ufloat('fitted_massErr'),
        fit_pt         = ufloat('fitted_pt'),
        fit_eta        = ufloat('fitted_eta'),
        fit_phi        = ufloat('fitted_phi'),
        alt_pk_mass    = ufloat('alt_pk_mass'),
        alt_pk_pt      = ufloat('alt_pk_pt'),
        alt_pk_eta     = ufloat('alt_pk_eta'),
        alt_pk_phi     = ufloat('alt_pk_phi'),
        alt_kp_mass    = ufloat('alt_kp_mass'),
        alt_kp_pt      = ufloat('alt_kp_pt'),
        alt_kp_eta     = ufloat('alt_kp_eta'),
        alt_kp_phi     = ufloat('alt_kp_phi'),
        alt_pp_mass    = ufloat('alt_pp_mass'),
        alt_pp_pt      = ufloat('alt_pp_pt'),
        alt_pp_eta     = ufloat('alt_pp_eta'),
        alt_pp_phi     = ufloat('alt_pp_phi'),

        # vertex
        vtx_x   = ufloat('vtx_x'),
        vtx_y   = ufloat('vtx_y'),
        vtx_z   = ufloat('vtx_z'),
        vtx_cxx = ufloat('vtx_cxx'),
        vtx_cyy = ufloat('vtx_cyy'),
        vtx_czz = ufloat('vtx_czz'),
        vtx_cyx = ufloat('vtx_cyx'),
        vtx_czx = ufloat('vtx_czx'),
        vtx_czy = ufloat('vtx_czy'),
        # post-fit tracks/leptons
        #trk3
        fit_trk3_pt  = ufloat('fitted_trk3_pt'),
        fit_trk3_eta = ufloat('fitted_trk3_eta'),
        fit_trk3_phi = ufloat('fitted_trk3_phi'),
        #trk4
        fit_trk4_pt  = ufloat('fitted_trk4_pt'),
        fit_trk4_eta = ufloat('fitted_trk4_eta'),
        fit_trk4_phi = ufloat('fitted_trk4_phi'),
        #lambda
        fit_lambda_pt  = ufloat('fitted_lambda_pt'),
        fit_lambda_eta = ufloat('fitted_lambda_eta'),
        fit_lambda_phi = ufloat('fitted_lambda_phi'),
        # isolation 
        trk3_iso04   = ufloat('trk3_iso04'),
        trk4_iso04   = ufloat('trk4_iso04'),
        lambda_iso04 = ufloat('lambda_iso04'),
        iso04        = ufloat('iso04'),
        iso04_ntrk        = ufloat('iso04_ntrk'),
        trk1_svip2d     = ufloat('trk1_svip2d'),
        trk1_svip2d_err = ufloat('trk1_svip2d_err'),
        trk2_svip2d     = ufloat('trk2_svip2d'),
        trk2_svip2d_err = ufloat('trk2_svip2d_err'),
        trk3_svip2d     = ufloat('trk3_svip2d'),
        trk3_svip2d_err = ufloat('trk3_svip2d_err'),
        trk4_svip2d     = ufloat('trk4_svip2d'),
        trk4_svip2d_err = ufloat('trk4_svip2d_err'),

    )
)

CountLambdabToLambdahh = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("LambdabToLambdahh")
)


LambdabToLambdahhBPHMCMatch = cms.EDProducer("MCMatcher",                  # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = LambdabToLambdahhTable.src,                           # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPH"),       # final mc-truth particle collection
    mcPdgId     = cms.vint32(5122),                             # one or more PDG ID (443 = J/psi); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(2),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

LambdabToLambdahhBPHMCTable = cms.EDProducer("CandMCMatchTableProducerBPH",
    recoObjects = LambdabToLambdahhTable.src,
    genParts    = cms.InputTag("finalGenParticlesBPH"),
    mcMap       = cms.InputTag("LambdabToLambdahhBPHMCMatch"),
    objName     = LambdabToLambdahhTable.name,
    objType     = cms.string("Other"),
    objBranchName = cms.string("genPart"),
    genBranchName = cms.string("LambdabToLambdahh"),
    docString   = cms.string("MC matching to status==2 LambdaB0"),
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

########################### Sequencies  ############################
LambdabToLambdahhSequence = cms.Sequence( LambdabToLambdahh  )
LambdabToLambdahhTables   = cms.Sequence( LambdabToLambdahhTable )
LambdabToLambdahhSequenceMC = cms.Sequence( LambdabToLambdahh + LambdabToLambdahhBPHMCMatch + Lambdab0Gen)
LambdabToLambdahhTablesMC   = cms.Sequence( LambdabToLambdahhTable + LambdabToLambdahhBPHMCTable + LambdaBGenmatchTable)



