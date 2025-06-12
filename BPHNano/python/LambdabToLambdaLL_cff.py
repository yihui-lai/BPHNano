import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *

########################### B-> K* ll ###########################

LambdabToLambdaMuMu = cms.EDProducer(
    'LambdabToLambdaLLBuilder',
    dileptons = cms.InputTag("MuMu:SelectedDiLeptons"),
    leptonTransientTracks = cms.InputTag('muonBPH', 'SelectedTransientMuons'),
    ditracks = cms.InputTag('LambdaToPPi', 'SelectedLambdaCollection'),
    transientTracks = cms.InputTag('tracksBPH', 'SelectedTransientTracks'),
    v0TransientTracks = cms.InputTag('LambdaToPPi', 'SelectedLambda'),
    PUtracks = cms.InputTag('tracksBPH', 'SelectedTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection  = cms.string('pt > 1.'
                                  '&& 5.2 < mass && mass < 6.0'),
    postVtxSelection = cms.string('5.2 < userFloat("fitted_mass") && userFloat("fitted_mass") < 6.0'
                                  '&& userFloat("sv_prob") > 0.001'
                                  '&& userFloat("fitted_cos_theta_2D") > 0.'),
    dileptonMassContraint = cms.double(3.0969)
)

########################### Tables ###########################

LambdabToLambdaMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("LambdabToLambdaMuMu"),
    cut       = cms.string(""),
    name      = cms.string("LambdabToLambdaMuMu"),
    doc       = cms.string("LambdabToLambdaMuMu Variables"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        # pre-fit quantities
        CandVars,
        l1_idx      = uint('l1_idx'),
        l2_idx      = uint('l2_idx'),
        trk1_idx    = uint('trk1_idx'),
        trk2_idx    = uint('trk2_idx'),
        lambda_idx   = uint('lambda_idx'),
        trk1_mass   = ufloat('trk1_mass'),
        trk2_mass   = ufloat('trk2_mass'),
        min_dr      = ufloat('min_dr'),
        max_dr      = ufloat('max_dr'),
        # vtx info
        chi2      = ufloat('sv_chi2'),
        svprob    = ufloat('sv_prob'),
        cos2D     = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        l_xy      = ufloat('l_xy'),
        l_xy_unc  = ufloat('l_xy_unc'),
        # post-fit momentum /masses
        mll_fullfit    = ufloat('fitted_mll'),
        mlambda_fullfit = ufloat('fitted_ditrack_mass'),
        fit_mass       = ufloat('fitted_mass'),
        fit_massErr    = ufloat('fitted_massErr'),
        fit_pt         = ufloat('fitted_pt'),
        fit_eta        = ufloat('fitted_eta'),
        fit_phi        = ufloat('fitted_phi'),
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
        #l1
        fit_l1_pt  = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        #l2
        fit_l2_pt  = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        #lambda
        fit_lambda_pt  = ufloat('fitted_lambda_pt'),
        fit_lambda_eta = ufloat('fitted_lambda_eta'),
        fit_lambda_phi = ufloat('fitted_lambda_phi'),
        # isolation 
        l1_iso04   = ufloat('l1_iso04'),
        l2_iso04   = ufloat('l2_iso04'),
        lambda_iso04 = ufloat('lambda_iso04'),

        trk1_svip2d     = ufloat('trk1_svip2d'),
        trk1_svip2d_err = ufloat('trk1_svip2d_err'),
        trk2_svip2d     = ufloat('trk2_svip2d'),
        trk2_svip2d_err = ufloat('trk2_svip2d_err'),
    )
)

CountLambdabToLambdaMuMu = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("LambdabToLambdaMuMu")
)

########################### Sequencies  ############################
LambdabToLambdaMuMuSequence = cms.Sequence( LambdabToLambdaMuMu  )
LambdabToLambdaMuMuTables   = cms.Sequence( LambdabToLambdaMuMuTable )
