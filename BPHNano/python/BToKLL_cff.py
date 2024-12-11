import FWCore.ParameterSet.Config as cms
from PhysicsTools.BPHNano.common_cff import *

########################### B-> K ll ###########################

BToKmumu = cms.EDProducer(
    'BToTrkLLBuilder',
    dileptons = cms.InputTag("JpsiToMuMu:SelectedDiLeptons"),
    leptonTransientTracks = cms.InputTag('muonBPH', 'SelectedTransientMuons'),
    dileptonKinVtxs = cms.InputTag('JpsiToMuMu:SelectedDiLeptonKinVtxs'),
    kaons = cms.InputTag('tracksBPH', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPH', 'SelectedTransientTracks'),
    trackMass = cms.double(493.677),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    PUtracks = cms.InputTag('tracksBPH', 'SelectedTracks'),
    preVtxSelection  = cms.string('pt > 5.0'
                                  '&& 4.5 < mass && mass < 6.'
                                  '&& userFloat("min_dr") > 0.03'),
    postVtxSelection = cms.string('5. < userFloat("fitted_mass") && userFloat("fitted_mass") < 6.' 
                                  '&& userFloat("sv_prob") > 0.001'
                                  '&& userFloat("fitted_cos_theta_2D") > 0.9'),
    dileptonMassContraint = cms.double(3.0969)
)

########################### Tables ###########################

BToKmumuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("BToKmumu"),
    cut       = cms.string(""),
    name      = cms.string("BToKMuMu"),
    doc       = cms.string("BToKMuMu Variable"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        # pre-fit quantities
        CandVars,
        l1_idx = uint('l1_idx'),
        l2_idx = uint('l2_idx'),
        k_idx  = uint('trk_idx'),
        minDR  = ufloat('min_dr'),
        maxDR  = ufloat('max_dr'),
        # fit and vtx info
        svprob   = ufloat('sv_prob'),
        l_xy     = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),

        vtx_x   = ufloat('vtx_x'),
        vtx_y   = ufloat('vtx_y'),
        vtx_z   = ufloat('vtx_z'),
        vtx_cxx = ufloat('vtx_cxx'),
        vtx_cyy = ufloat('vtx_cyy'),
        vtx_czz = ufloat('vtx_czz'),
        vtx_cyx = ufloat('vtx_cyx'),
        vtx_czx = ufloat('vtx_czx'),
        vtx_czy = ufloat('vtx_czy'),
        beamspot_x  = ufloat('beamspot_x'),
        beamspot_y  = ufloat('beamspot_y'),
        beamspot_z  = ufloat('beamspot_z'),
        # Mll
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D     = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass    = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt      = ufloat('fitted_pt'),
        fit_eta     = ufloat('fitted_eta'),
        fit_phi     = ufloat('fitted_phi'),
        fit_l1_pt   = ufloat('fitted_l1_pt'),
        fit_l1_eta  = ufloat('fitted_l1_eta'),
        fit_l1_phi  = ufloat('fitted_l1_phi'),
        fit_l2_pt   = ufloat('fitted_l2_pt'),
        fit_l2_eta  = ufloat('fitted_l2_eta'),
        fit_l2_phi  = ufloat('fitted_l2_phi'),
        fit_k_pt    = ufloat('fitted_trk_pt'),
        fit_k_eta   = ufloat('fitted_trk_eta'),
        fit_k_phi   = ufloat('fitted_trk_phi'),
        k_svip2d    = ufloat('k_svip2d'),
        k_svip2d_err = ufloat('k_svip2d_err'),
        l1_iso04    = ufloat('l1_iso04'),
        l2_iso04    = ufloat('l2_iso04'),
        k_iso04     = ufloat('trk_iso04'),
        # constraint properties
        constraint_sv_prob  = ufloat('constraint_sv_prob'),
        constraint_mass     = ufloat('constraint_mass'),
        constraint_massErr  = ufloat('constraint_massErr'),
        constraint_mll      = ufloat('constraint_mll')
    )
)

CountBToKmumu = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src       = cms.InputTag("BToKmumu")
)    

########################### Sequencies  ############################
BToKMuMuSequence = cms.Sequence(BToKmumu)
BToKMuMuTables   = cms.Sequence(BToKmumuTable)