/// original authors: Yihui Lai (Princeton)
// takes 5 tracks and make a B
// # B+ -> D0 K*(892)+, with D0->K- pi+, K*+->K_S0 pi+, K_S0->pi+ pi-


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <memory>
#include <typeinfo>
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "helper.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include <memory>
#include "KinVtxFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include <tuple>  // For std::tuple

#include "PhysicsTools/BPHNano/plugins/XGBooster.h"

#include <iostream>
#include <chrono>

// pdg mass constants
namespace {

  const double KS0_MASS = 0.497611;
  const double KS0_MASS_window = 0.02;
  const double D0_MASS = 1.86484;
  const double D0_MASS_window = 0.02;
  const double KSTAR_MASS = 0.89166;
  const double KSTAR_MASS_window = 0.075;
  const double BPLUS_MASS = 5.27934;
  const double BPLUS_MASS_window = 0.3;

}  // namespace

class BtoD0KstarProducer : public edm::stream::EDProducer<> {
public:
  explicit BtoD0KstarProducer(const edm::ParameterSet &theParameters):
      bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      token_beamSpot(consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"))),
      token_vertices(consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("tracks"))),
      lostTracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("lostTracks"))),
      tkNHitsCut_( theParameters.getParameter<int>("tkNHitsCut")),
      minTrackPt_( theParameters.getParameter<double>("minTrackPt")),
      maxTrackEta_( theParameters.getParameter<double>("maxTrackEta")),
      tkChi2Cut_( theParameters.getParameter<double>("tkChi2Cut")),
      tkIPSigXYCut_( theParameters.getParameter<double>("tkIPSigXYCut")),
      vtxChi2Cut_( theParameters.getParameter<double>("vtxChi2Cut")),
      vtxDecaySigXYCut_( theParameters.getParameter<double>("vtxDecaySigXYCut")),
      TrkSigXYCut_( theParameters.getParameter<double>("TrkSigXYCut")),
      vtxDecaySigXYZCut_( theParameters.getParameter<double>("vtxDecaySigXYZCut")),
      cosThetaXYCut_( theParameters.getParameter<double>("cosThetaXYCut")),
      cosThetaXYZCut_( theParameters.getParameter<double>("cosThetaXYZCut")),
      DtkPtCut_( theParameters.getParameter<double>("DtkPtCut")),
      diTrack2_dca_( theParameters.getParameter<double>("diTrack2_dca")),
      Trk34SigXYCut_( theParameters.getParameter<double>("Trk34SigXYCut")),
      Ks0_l_xyzSigCut_( theParameters.getParameter<double>("Ks0_l_xyzSigCut")),
      D0_PtCut_( theParameters.getParameter<double>("D0_PtCut")),
      D0vtxDecaySigXYCut_( theParameters.getParameter<double>("D0vtxDecaySigXYCut")),
      B_PtCut_( theParameters.getParameter<double>("B_PtCut")),
      Btrk_dcaSigCut_( theParameters.getParameter<double>("Btrk_dcaSigCut")),
      verbose( theParameters.getParameter<int>("verbose"))
        {
          produces<pat::CompositeCandidateCollection>("D0");
          produces<pat::CompositeCandidateCollection>("Kstar");
          produces<pat::CompositeCandidateCollection>("B");
  }

  ~BtoD0KstarProducer() override {}

  virtual void produce(edm::Event&, const edm::EventSetup&);
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  const edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTracksToken_;
  const int tkNHitsCut_;
  const double minTrackPt_;
  const double maxTrackEta_;
  const double tkChi2Cut_;
  const double tkIPSigXYCut_;
  const double vtxChi2Cut_;
  const double vtxDecaySigXYCut_;
  const double TrkSigXYCut_;
  const double vtxDecaySigXYZCut_;
  const double cosThetaXYCut_;
  const double cosThetaXYZCut_;
  const double DtkPtCut_;
  const double diTrack2_dca_;
  const double Trk34SigXYCut_;
  const double Ks0_l_xyzSigCut_;
  const double D0_PtCut_;
  const double D0vtxDecaySigXYCut_;
  const double B_PtCut_;
  const double Btrk_dcaSigCut_;
  const int verbose;

};

void BtoD0KstarProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using std::vector;
    
    // -------------------------
    // 1) Get Event Objects
    // -------------------------
    
    edm::Handle<reco::BeamSpot> theBeamSpotHandle;
    iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
    if (!theBeamSpotHandle.isValid()) return;
    const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
    const auto& theMagneticField = iSetup.getData(bFieldToken_);

    const auto& magField = iSetup.getData(bFieldToken_);

    edm::Handle<std::vector<reco::Vertex>> vtxH;
    iEvent.getByToken(token_vertices, vtxH);
    if (vtxH->empty()) return;
    const reco::Vertex& pv = (*vtxH)[0];
    math::XYZPoint pvPos(pv.position());

    edm::Handle<edm::View<pat::PackedCandidate>> tracksH;
    edm::Handle<edm::View<pat::PackedCandidate>> lostTracksH;
    iEvent.getByToken(tracksToken_, tracksH);
    iEvent.getByToken(lostTracksToken_, lostTracksH);
    if (!tracksH.isValid() || !lostTracksH.isValid()) return;

    // -------------------------
    // 2) Preselect good tracks (store index!)
    // -------------------------
    struct GoodTrack {
        const pat::PackedCandidate* pc;
        reco::TransientTrack tt;
        size_t index; // index into tracksH
    };
    vector<GoodTrack> goodTracks;
    goodTracks.reserve(tracksH->size());

    for (size_t i = 0; i < tracksH->size(); ++i) {
        const auto& pc = (*tracksH)[i];
        if (!pc.hasTrackDetails() || !pc.bestTrack()) continue;
        if (pc.charge() == 0) continue;
        if (pc.bestTrack()->hitPattern().numberOfValidHits() < tkNHitsCut_) continue;
        if (pc.pt() < minTrackPt_) continue;
        if (std::abs(pc.eta()) > maxTrackEta_) continue;
        if (!pc.trackHighPurity()) continue;
        if (pc.bestTrack()->normalizedChi2() > tkChi2Cut_) continue;
        const reco::Track& pseu = pc.pseudoTrack();
        double ipsigXY = std::abs(pseu.dxy(pvPos) / pseu.dxyError());
        if (ipsigXY < tkIPSigXYCut_) continue;

        goodTracks.push_back({&pc, reco::TransientTrack(*pc.bestTrack(), &magField), i});
    }

    // -------------------------
    // helpers
    // -------------------------
    auto dca2tracks = [&](const reco::TransientTrack& t1, const reco::TransientTrack& t2, GlobalPoint& cxPt) -> float {
        const auto& s1 = t1.impactPointTSCP();
        const auto& s2 = t2.impactPointTSCP();
        if (!s1.isValid() || !s2.isValid()) return -1;
        ClosestApproachInRPhi cApp;
        cApp.calculate(s1.theState(), s2.theState());
        if (!cApp.status()) return -1;
        cxPt = cApp.crossingPoint();
        return std::abs(cApp.distance());
    };
    auto goodPOCA = [&](const GlobalPoint& p) {
        double r2 = p.x()*p.x() + p.y()*p.y();
        return (r2 < 120.*120. && std::abs(p.z()) < 300.);
    };
    // -------------------------
    // containers for intermediate candidates
    // store indices of original tracks for overlap checks
    // -------------------------
    struct BuiltCand {
        pat::CompositeCandidate cand;
        reco::TransientTrack fittedTT;       // fitted composite transient track (from KinVtxFitter)
        reco::TransientTrack fittedTTMC;       // fitted composite transient track (from KinVtxFitter)
        std::vector<size_t> origTrackIdx;    // indices in goodTracks
    };
    vector<BuiltCand> d0List;    d0List.reserve(32);
    vector<BuiltCand> kstarList; kstarList.reserve(32);

    // --- Output collections ---
    //auto ksCands_out = std::make_unique<pat::CompositeCandidateCollection>();
    auto d0Cands_out = std::make_unique<pat::CompositeCandidateCollection>();
    auto kstarCands_out = std::make_unique<pat::CompositeCandidateCollection>();
    auto bCands_out = std::make_unique<pat::CompositeCandidateCollection>();


    // -------------------------
    // 3) Loop pairs -> K* or D0 candidates
    //    (preselection + unconstrained fit, then mass-constrained fit)
    // -------------------------

    struct D0Hypothesis { double m1; double m2; std::string n1,n2; };
    /*std::vector<D0Hypothesis> D0modes = {
        {K_MASS, PI_MASS, "Kminus","piplus"},
        {K_MASS, K_MASS, "Kminus","Kplus"},
        {PI_MASS, PI_MASS, "piminus","piplus"},
        {PI_MASS, K_MASS, "piminus","Kplus"}
    };*/

    std::vector<D0Hypothesis> D0modes = {
        {K_MASS, PI_MASS, "Kminus","piplus"},
        {PI_MASS, K_MASS, "piminus","Kplus"}
    };
    for (size_t ia = 0; ia < goodTracks.size(); ++ia) {
        for (size_t ib = ia + 1; ib < goodTracks.size(); ++ib) {
            const auto& g1 = goodTracks[ia];
            const auto& g2 = goodTracks[ib];
            if (g1.pc->charge() * g2.pc->charge() != -1) continue;
            // quick DCA and POCA check
            GlobalPoint cxPt;
            float dca = dca2tracks(g1.tt, g2.tt, cxPt);
            if (dca < 0 || !goodPOCA(cxPt)) continue;

            // quick momentum dot product to ensure not back-to-back
            auto ts1 = g1.tt.trajectoryStateClosestToPoint(cxPt);
            auto ts2 = g2.tt.trajectoryStateClosestToPoint(cxPt);
            if (!ts1.isValid() || !ts2.isValid()) continue;
            if (ts1.momentum().dot(ts2.momentum()) < 0) continue;
	
            // basic pi-pi mass (unconstrained) to check for Ks0
            math::PtEtaPhiMLorentzVector p4_1(g1.pc->pt(), g1.pc->eta(), g1.pc->phi(), PI_MASS);
            math::PtEtaPhiMLorentzVector p4_2(g2.pc->pt(), g2.pc->eta(), g2.pc->phi(), PI_MASS);
            auto p4_pi_pi = p4_1 + p4_2;
            bool tryKs = (std::abs(p4_pi_pi.M() - KS0_MASS) < KS0_MASS_window);
            //distance closest approach in x,y wrt beam spot
            std::pair<double, double> DCA_pos_beamspot = computeDCA(g1.tt, *theBeamSpot);
            std::pair<double, double> DCA_neg_beamspot = computeDCA(g2.tt, *theBeamSpot);
            if(DCA_pos_beamspot.second==0 || DCA_neg_beamspot.second==0) continue;
            if (tryKs) {
                if( abs(DCA_pos_beamspot.first/DCA_pos_beamspot.second) < TrkSigXYCut_ && abs(DCA_neg_beamspot.first/DCA_neg_beamspot.second) < TrkSigXYCut_) continue;

                // 1) Unconstrained fit
                KinVtxFitter ksFitUnc({g1.tt, g2.tt}, {PI_MASS, PI_MASS}, {PI_SIGMA, PI_SIGMA});
                if (!ksFitUnc.success()) { continue; }
                else {
                    // Optional: apply pre/post-fit cuts
                    if (ksFitUnc.prob() > 1e-4 && std::abs(ksFitUnc.fitted_p4().mass() - KS0_MASS) < KS0_MASS_window) {
                        // 2) Mass-constrained fit (apply constraint to KS0)
                        KinVtxFitter ksFitMC;
                        try {
                            ksFitMC = KinVtxFitter({g1.tt, g2.tt}, {PI_MASS, PI_MASS}, {PI_SIGMA, PI_SIGMA}, KS0_MASS);
                        } catch (const VertexException &e) {
			    ksFitMC = ksFitUnc;
                        }

			// 
                        // -------------------------
                        // 4) Build K* from Ks0 + a bachelor pion
                        // -------------------------
                        for (size_t k = 0; k < goodTracks.size(); ++k) {
                            // bachelor should not be one of ks daughters
                            if ( ia==k || ib==k) continue;

                            const auto& btrk = goodTracks[k];

                            auto p4_kstar = ksFitUnc.fitted_p4() + btrk.pc->p4();
                            if (std::abs(p4_kstar.M() - KSTAR_MASS) > KSTAR_MASS_window) continue;
			    float ks_mass_unc = sqrt(ksFitUnc.fitted_candidate().kinematicParametersError().matrix()(6,6));
    			    // K* unconstrained fit: (bachelor PI, ks fitted TT)
                            KinVtxFitter kstarUnc({btrk.tt, ksFitUnc.fitted_candidate_ttrk()}, {PI_MASS, KS0_MASS}, {PI_SIGMA, ks_mass_unc});
                            //KinVtxFitter kstarUnc({btrk.tt, ksFitMC.fitted_candidate_ttrk()}, {PI_MASS, KS0_MASS}, {PI_SIGMA, ks_mass_unc});
                            if (!kstarUnc.success()) continue;
                            if (kstarUnc.prob() < 1e-4) continue;
                            // mass-constrained fit for K*
                            KinVtxFitter kstarMC;
                            try {
                                kstarMC = KinVtxFitter({btrk.tt, ksFitMC.fitted_candidate_ttrk()}, {PI_MASS, KS0_MASS}, {PI_SIGMA,ks_mass_unc}, KSTAR_MASS);
                            } catch (const VertexException &e) {
                                kstarMC = kstarUnc;
                            }
                            if (!kstarMC.success()) continue;
                    
	                    auto Ks0_Kstar_lxyz =  l_xyz(ksFitMC, kstarMC.fitted_vtx().x(), kstarMC.fitted_vtx().y(), kstarMC.fitted_vtx().z());
                            double Ks0_Kin_d0_l_xyzSig = Ks0_Kstar_lxyz.error()==0? -99: abs(Ks0_Kstar_lxyz.value()/Ks0_Kstar_lxyz.error());
                            if(Ks0_Kin_d0_l_xyzSig < Ks0_l_xyzSigCut_) continue;

			    pat::CompositeCandidate kstar;
                            kstar.addUserFloat("ks0_pt", ksFitUnc.fitted_p4().pt() );
                            kstar.addUserFloat("ks0_eta", ksFitUnc.fitted_p4().eta() );
                            kstar.addUserFloat("ks0_phi", ksFitUnc.fitted_p4().phi() );
                            kstar.addUserFloat("ks0_mass", ksFitUnc.fitted_candidate().mass());
                            kstar.addUserFloat("ks0_massErr", sqrt(ksFitUnc.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                            kstar.addUserFloat("ks0_chi2", ksFitUnc.chi2());
                            kstar.addUserFloat("ks0_ndof", ksFitUnc.dof());
                            kstar.addUserFloat("ks0_prob", ksFitUnc.prob());
                            kstar.addUserFloat("ks0_vx", ksFitUnc.fitted_vtx().x());
                            kstar.addUserFloat("ks0_vy", ksFitUnc.fitted_vtx().y());
                            kstar.addUserFloat("ks0_vz", ksFitUnc.fitted_vtx().z());
                            kstar.addUserFloat("ks0_MC_pt", ksFitMC.fitted_p4().pt() );
                            kstar.addUserFloat("ks0_MC_eta", ksFitMC.fitted_p4().eta() );
                            kstar.addUserFloat("ks0_MC_phi", ksFitMC.fitted_p4().phi() );
                            kstar.addUserFloat("ks0_MC_mass", ksFitMC.fitted_candidate().mass());
                            kstar.addUserFloat("ks0_MC_massErr", sqrt(ksFitMC.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                            kstar.addUserFloat("ks0_MC_chi2", ksFitMC.chi2());
                            kstar.addUserFloat("ks0_MC_ndof", ksFitMC.dof());
                            kstar.addUserFloat("ks0_MC_prob", ksFitMC.prob());
                            kstar.addUserFloat("ks0_MC_vx", ksFitMC.fitted_vtx().x());
                            kstar.addUserFloat("ks0_MC_vy", ksFitMC.fitted_vtx().y());
                            kstar.addUserFloat("ks0_MC_vz", ksFitMC.fitted_vtx().z());

                            kstar.setP4(kstarUnc.fitted_p4());
                            kstar.addUserFloat("mass", kstarUnc.fitted_candidate().mass());
                            kstar.addUserFloat("massErr", sqrt(kstarUnc.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                            kstar.addUserFloat("chi2", kstarUnc.chi2());
                            kstar.addUserFloat("prob", kstarUnc.prob());
                            kstar.addUserInt("ks_trk1_idx", (int)ia);
                            kstar.addUserInt("ks_trk2_idx", (int)ib);
                            kstar.addUserInt("bachelor_idx", (int)k);
                            kstar.addUserFloat("MC_pt",      ksFitMC.fitted_p4().pt() );
                            kstar.addUserFloat("MC_eta",     ksFitMC.fitted_p4().eta() );
                            kstar.addUserFloat("MC_phi",     ksFitMC.fitted_p4().phi() );
			    kstar.addUserFloat("MC_mass",    ksFitMC.fitted_candidate().mass());
                            kstar.addUserFloat("MC_massErr", sqrt(ksFitMC.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                            kstar.addUserFloat("MC_chi2", ksFitMC.chi2());
                            kstar.addUserFloat("MC_prob", ksFitMC.prob());
			    kstarCands_out->push_back(kstar);
                            
			    BuiltCand bkstar;
                            bkstar.cand = kstar;
                            bkstar.fittedTT = ksFitUnc.fitted_candidate_ttrk();
                            bkstar.fittedTTMC = kstarMC.fitted_candidate_ttrk();
                            bkstar.origTrackIdx = {ia, ib, k};
                            kstarList.push_back(std::move(bkstar));
                        }

                    }
                    if(verbose>=2) std::cout<< "Pass K*" << std::endl;
                }
            } // end tryKs
	else if( g1.pc->pt() > DtkPtCut_ && g2.pc->pt() > DtkPtCut_ && abs(DCA_pos_beamspot.first/DCA_pos_beamspot.second) > Trk34SigXYCut_ && abs(DCA_neg_beamspot.first/DCA_neg_beamspot.second) > Trk34SigXYCut_) {
	        // -----------------
                // D0 hypotheses (K pi and pi K) - keep minimal modes for speed
                // -----------------
                for (const auto &hyp : D0modes) {
    		    math::PtEtaPhiMLorentzVector p4_d1(g1.pc->pt(), g1.pc->eta(), g1.pc->phi(), hyp.m1);
                    math::PtEtaPhiMLorentzVector p4_d2(g2.pc->pt(), g2.pc->eta(), g2.pc->phi(), hyp.m2);
                    auto p4_d = p4_d1 + p4_d2;
                    if (std::abs(p4_d.M() - D0_MASS) > D0_MASS_window) continue;

                    // apply simple kinematic pre-cuts (D0 daughters pT etc)
                    if (g1.pc->pt() < DtkPtCut_ || g2.pc->pt() < DtkPtCut_) continue;
                    // optionally require dca small enough
                    if (dca > diTrack2_dca_) continue;
                    // unconstrained fit
                    KinVtxFitter d0Unc({g1.tt, g2.tt}, {hyp.m1, hyp.m2}, {PI_SIGMA, PI_SIGMA});
                    if (!d0Unc.success()) continue;
                    if (d0Unc.prob() < 1e-4) continue;
                    if (std::abs(d0Unc.fitted_p4().mass() - D0_MASS) > D0_MASS_window) continue;
                    // mass-constrained fit to D0 mass
                    KinVtxFitter d0MC;
                    try {
                        d0MC = KinVtxFitter({g1.tt, g2.tt}, {hyp.m1, hyp.m2}, {PI_SIGMA, PI_SIGMA}, D0_MASS);
                    } catch (const VertexException &e) {
                        d0MC = d0Unc;
                    }
                    //if (!d0MC.success()) continue;
                    pat::CompositeCandidate d0;
                    d0.setP4(d0Unc.fitted_p4());
                    d0.addUserFloat("mass", d0Unc.fitted_candidate().mass());
                    d0.addUserFloat("massErr", sqrt(d0Unc.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                    d0.addUserFloat("chi2", d0Unc.chi2());
                    d0.addUserFloat("ndof", d0Unc.dof());
                    d0.addUserFloat("prob", d0Unc.prob());
                    d0.addUserFloat("MC_pt", d0MC.fitted_p4().pt() );
                    d0.addUserFloat("MC_eta", d0MC.fitted_p4().eta() );
                    d0.addUserFloat("MC_phi", d0MC.fitted_p4().phi() );
                    d0.addUserFloat("MC_mass", d0MC.fitted_candidate().mass());
                    d0.addUserFloat("MC_massErr", sqrt(d0MC.fitted_candidate().kinematicParametersError().matrix()(6,6)));
                    d0.addUserFloat("MC_chi2", d0MC.chi2());
                    d0.addUserFloat("MC_ndof", d0MC.dof());
                    d0.addUserFloat("MC_prob", d0MC.prob());
                    d0.addUserInt("hyp_m1_idx", (int)ia);
                    d0.addUserInt("hyp_m2_idx", (int)ib);
                    d0Cands_out->push_back(d0);

                    BuiltCand bd;
                    bd.cand = d0;
                    bd.fittedTT = d0Unc.fitted_candidate_ttrk();
                    bd.fittedTTMC = d0MC.fitted_candidate_ttrk();
                    bd.origTrackIdx = {ia, ib};
                    d0List.push_back(std::move(bd));
                } // end D0 modes
	    }
        } // ib
    } // ia

    

    // -------------------------
    // 5) Build B from D0 + K* (check overlap)
    // -------------------------
    // Combine every d0 with every kstar, skipping overlap
    for (size_t id0 = 0; id0 < d0List.size(); ++id0) {
        const auto& bd = d0List[id0];
        for (size_t ik = 0; ik < kstarList.size(); ++ik) {
            const auto& bk = kstarList[ik];

            // overlap check
            bool overlap = false;
            for (size_t idx_d : bd.origTrackIdx)
                if (std::find(bk.origTrackIdx.begin(), bk.origTrackIdx.end(), idx_d) != bk.origTrackIdx.end()) {
                    overlap = true;
                    break;
                }
            if (overlap) continue;

            // quick mass preselection
            math::XYZTLorentzVector p4B = bd.cand.p4() + bk.cand.p4();
            if (std::abs(p4B.M() - BPLUS_MASS) > BPLUS_MASS_window) continue;

            // B fit: combine fitted transient tracks of K* and D0
            KinVtxFitter bUnc({bk.fittedTTMC, bd.fittedTTMC}, {KSTAR_MASS, D0_MASS}, {bk.cand.userFloat("massErr"), bd.cand.userFloat("massErr")});
            if (!bUnc.success()) continue;
            if (bUnc.prob() < 1e-4) continue;

            // mass-constrained B fit (optional; here we do it to improve mass resolution)
            KinVtxFitter bMC;
            try {
                bMC = KinVtxFitter({bk.fittedTTMC, bd.fittedTTMC}, {KSTAR_MASS, D0_MASS}, {bk.cand.userFloat("massErr"), bd.cand.userFloat("massErr")}, BPLUS_MASS);
            } catch (const VertexException &e) {
                bMC = bUnc;
            }
            //if (!bMC.success()) continue;

            pat::CompositeCandidate B;
            B.setP4(bUnc.fitted_p4());
            B.addUserFloat("rawmass", p4B.M());
            B.addUserFloat("mass", bUnc.fitted_candidate().mass());
            B.addUserFloat("massErr", sqrt(bUnc.fitted_candidate().kinematicParametersError().matrix()(6,6)));
            B.addUserFloat("chi2", bUnc.chi2());
            B.addUserFloat("prob", bUnc.prob());
            B.addUserFloat("lxyz", l_xyz(bUnc, pvPos.x(), pvPos.y(), pvPos.z()).value());
            B.addUserFloat("lxyzSig", l_xyz(bUnc, pvPos.x(), pvPos.y(), pvPos.z()).error()==0 ? -99 :
                           std::abs(l_xyz(bUnc, pvPos.x(), pvPos.y(), pvPos.z()).value()/l_xyz(bUnc, pvPos.x(), pvPos.y(), pvPos.z()).error()));

	    B.addUserFloat("MC_pt", bMC.fitted_p4().pt() );
  	    B.addUserFloat("MC_eta", bMC.fitted_p4().eta() );  
	    B.addUserFloat("MC_phi", bMC.fitted_p4().phi() );
	    B.addUserFloat("MC_mass", bMC.fitted_candidate().mass());
            B.addUserFloat("MC_massErr", sqrt(bMC.fitted_candidate().kinematicParametersError().matrix()(6,6)));
            B.addUserFloat("MC_chi2", bMC.chi2());
            B.addUserFloat("MC_prob", bMC.prob());
            B.addUserFloat("MC_lxyz", l_xyz(bMC, pvPos.x(), pvPos.y(), pvPos.z()).value());
            B.addUserFloat("MC_lxyzSig", l_xyz(bMC, pvPos.x(), pvPos.y(), pvPos.z()).error()==0 ? -99 :
                           std::abs(l_xyz(bMC, pvPos.x(), pvPos.y(), pvPos.z()).value()/l_xyz(bMC, pvPos.x(), pvPos.y(), pvPos.z()).error()));

            B.addUserInt("d0_idx", (int)id0);
            B.addUserInt("kstar_idx", (int)ik);
            // store also some intermediate masses
            B.addUserFloat("D0_mass", bd.cand.userFloat("mass"));
            B.addUserFloat("Kstar_mass", bk.cand.userFloat("mass"));
            bCands_out->push_back(B);
        }
    }

    iEvent.put(std::move(d0Cands_out), "D0");
    iEvent.put(std::move(kstarCands_out), "Kstar");
    iEvent.put(std::move(bCands_out), "B");

}

DEFINE_FWK_MODULE(BtoD0KstarProducer);

