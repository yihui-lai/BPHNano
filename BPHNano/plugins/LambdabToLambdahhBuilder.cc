////////////////////////////// LambdabToLambdahhBuilder //////////////////////////////
/// authors: Y Lai (Princeton)
// takes the ditrack collection and a dilepton collection and produces Lambda_b moth
// - ers using a four-track vertex

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"

class LambdabToLambdahhBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit LambdabToLambdahhBuilder(const edm::ParameterSet &cfg):
    bFieldToken_{esConsumes<MagneticField, IdealMagneticFieldRecord>()},
    // selections
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    //inputs
    lambda0_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("lambda0") )},
    v0_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("v0TransientTracks") )},
    tracks_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    transientTracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )}
  {
    //output
    produces<pat::CompositeCandidateCollection>();
  }

  ~LambdabToLambdahhBuilder() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fihhDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;

  // inputs
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> lambda0_;
  const edm::EDGetTokenT<TransientTrackCollection> v0_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracks_;
  const edm::EDGetTokenT<TransientTrackCollection> transientTracks_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;

};

void LambdabToLambdahhBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> lambda;
  evt.getByToken(lambda0_, lambda);
  edm::Handle<TransientTrackCollection> v0_ttracks;
  evt.getByToken(v0_ttracks_, v0_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> tracks;
  evt.getByToken(tracks_, tracks);
  edm::Handle<TransientTrackCollection> transientTracks;
  evt.getByToken(transientTracks_, transientTracks);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);

  edm::ESHandle<MagneticField> fieldHandle;
  const auto& bField = iSetup.getData(bFieldToken_);
  AnalyticalImpactPointExtrapolator extrapolator(&bField);
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());

  for (size_t lambda_idx = 0; lambda_idx < lambda->size(); ++lambda_idx) {
    edm::Ptr<pat::CompositeCandidate> lambda_ptr(lambda, lambda_idx);
    edm::Ptr<reco::Candidate> trk1_ptr = lambda_ptr->userCand("trk1");
    edm::Ptr<reco::Candidate> trk2_ptr = lambda_ptr->userCand("trk2");
    int trk1_idx = lambda_ptr->userInt("trk1_idx");
    int trk2_idx = lambda_ptr->userInt("trk2_idx");


    for (size_t trk3_idx = 0; trk3_idx < tracks->size(); ++trk3_idx ) {
      edm::Ptr<pat::CompositeCandidate> trk3_ptr( tracks, trk3_idx );

      for (size_t trk4_idx = trk3_idx + 1; trk4_idx < tracks->size(); ++trk4_idx) {
        edm::Ptr<pat::CompositeCandidate> trk4_ptr( tracks, trk4_idx );

        if(trk1_idx == int(trk3_idx) || trk2_idx == int(trk3_idx) || trk2_idx == int(trk4_idx) || trk2_idx == int(trk4_idx) ) continue;

        // more cut on trk 3,4
        if (trk3_ptr->charge() + trk4_ptr->charge() != 0) continue;
	const auto& tPosTrk = (trk3_ptr->charge() > 0) ? transientTracks->at(trk3_idx) : transientTracks->at(trk4_idx);
        const auto& tNegTrk = (trk3_ptr->charge() < 0) ? transientTracks->at(trk3_idx) : transientTracks->at(trk4_idx);
        const auto& posImpact = tPosTrk.impactPointTSCP();
        const auto& negImpact = tNegTrk.impactPointTSCP();
        if (!posImpact.isValid() || !negImpact.isValid()) continue;
        FreeTrajectoryState const& posState = posImpact.theState();
        FreeTrajectoryState const& negState = negImpact.theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(posState, negState);
        if (!cApp.status()) continue;
        float dca = std::abs(cApp.distance());
        GlobalPoint cxPt = cApp.crossingPoint();
        double cxPtR = sqrt(cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y());
        if (cxPtR > 120.  || std::abs(cxPt.z()) > 300.) continue;

        // lambda_b candidate
        pat::CompositeCandidate cand;
        math::PtEtaPhiMLorentzVector trk3_p4(
            trk3_ptr->pt(),
            trk3_ptr->eta(),
            trk3_ptr->phi(),
            K_MASS
        );
        math::PtEtaPhiMLorentzVector trk4_p4(
            trk4_ptr->pt(),
            trk4_ptr->eta(),
            trk4_ptr->phi(),
            K_MASS
        );
        math::PtEtaPhiMLorentzVector trk3_pion_p4(
            trk3_ptr->pt(),
            trk3_ptr->eta(),
            trk3_ptr->phi(),
            PI_MASS
        );
        math::PtEtaPhiMLorentzVector trk4_pion_p4(
            trk4_ptr->pt(),
            trk4_ptr->eta(),
            trk4_ptr->phi(),
            PI_MASS
        );
        math::PtEtaPhiMLorentzVector lambda_p4(
            lambda_ptr->pt(),
            lambda_ptr->eta(),
            lambda_ptr->phi(),
            LAMBDA_MASS
        );
	auto lambdab_kk = trk3_p4 + trk4_p4 + lambda_p4;
	auto lambdab_kp = trk3_pion_p4 + trk4_p4 + lambda_p4;
	auto lambdab_pk = trk3_p4 + trk4_pion_p4 + lambda_p4;
	auto lambdab_pp = trk3_pion_p4 + trk4_pion_p4 + lambda_p4;
        cand.addUserFloat("alt_pk_pt"  , lambdab_pk.pt());
        cand.addUserFloat("alt_pk_eta" , lambdab_pk.eta());
        cand.addUserFloat("alt_pk_phi" , lambdab_pk.phi());
        cand.addUserFloat("alt_pk_mass", lambdab_pk.mass());
        cand.addUserFloat("alt_kp_pt"  , lambdab_kp.pt());
        cand.addUserFloat("alt_kp_eta" , lambdab_kp.eta());
        cand.addUserFloat("alt_kp_phi" , lambdab_kp.phi());
        cand.addUserFloat("alt_kp_mass", lambdab_kp.mass());
        cand.addUserFloat("alt_pp_pt"  , lambdab_pp.pt());
        cand.addUserFloat("alt_pp_eta" , lambdab_pp.eta());
        cand.addUserFloat("alt_pp_phi" , lambdab_pp.phi());
        cand.addUserFloat("alt_pp_mass", lambdab_pp.mass());

	cand.setP4(lambdab_kk);
	cand.setCharge(trk3_ptr->charge() + trk4_ptr->charge());

        // save daughters - unfitted
        cand.addUserCand("trk1", trk1_ptr);
        cand.addUserCand("trk2", trk2_ptr);
        cand.addUserCand("lambda", lambda_ptr);
        cand.addUserCand("trk3", trk3_ptr);
        cand.addUserCand("trk4", trk3_ptr);
  
        // save indices
        cand.addUserInt("trk1_idx", trk1_idx);
        cand.addUserInt("trk2_idx", trk2_idx);
        cand.addUserInt("lambda_idx" , lambda_idx);
        cand.addUserInt("trk3_idx", trk3_idx);
        cand.addUserInt("trk4_idx", trk4_idx);
        cand.addUserFloat("trk1_mass", lambda_ptr->userFloat("trk1_mass"));
        cand.addUserFloat("trk2_mass", lambda_ptr->userFloat("trk2_mass"));
        cand.addUserFloat("trk12_dr", reco::deltaR(*trk1_ptr, *trk2_ptr));
        cand.addUserFloat("trk34_dr", reco::deltaR(*trk3_ptr, *trk4_ptr));
        cand.addUserFloat("trk34_dca", dca);
        cand.addUserFloat("trk34_cxPtR", cxPtR);
        auto dr_info = min_max_dr({trk3_ptr, trk4_ptr, trk1_ptr, trk2_ptr});
        cand.addUserFloat("trk1234_min_dr", dr_info.first);
        cand.addUserFloat("trk1234_max_dr", dr_info.second);
  
  
        // check if pass pre vertex cut
        if ( !pre_vtx_selection_(cand) ) continue;
 
        KinVtxFitter fitter_kk(
          { transientTracks->at(trk3_idx), transientTracks->at(trk4_idx), v0_ttracks->at(lambda_idx) },
          { K_MASS, K_MASS, LAMBDA_MASS },
          { K_SIGMA, K_SIGMA, lambda_ptr->userFloat("fitted_massErr") }
          );
	if (!fitter_kk.success()) continue;

        // vertex position
        cand.setVertex(
          reco::Candidate::Point(
            fitter_kk.fitted_vtx().x(),
            fitter_kk.fitted_vtx().y(),
            fitter_kk.fitted_vtx().z()
          )
        );
  
        // vertex vars
        cand.addUserFloat("sv_chi2", fitter_kk.chi2());
        cand.addUserFloat("sv_ndof", fitter_kk.dof());
        cand.addUserFloat("sv_prob", fitter_kk.prob());
  
        // refitted kinematic vars
        cand.addUserFloat("fitted_mhh", (fitter_kk.daughter_p4(0) + fitter_kk.daughter_p4(1)).mass());
        cand.addUserFloat("fitted_mlambda_h1", (fitter_kk.daughter_p4(2) + fitter_kk.daughter_p4(0)).mass());
        cand.addUserFloat("fitted_mlambda_h2", (fitter_kk.daughter_p4(2) + fitter_kk.daughter_p4(1)).mass());

  
        auto fit_p4 = fitter_kk.fitted_p4();
        cand.addUserFloat("fitted_pt"  , fit_p4.pt());
        cand.addUserFloat("fitted_eta" , fit_p4.eta());
        cand.addUserFloat("fitted_phi" , fit_p4.phi());
        cand.addUserFloat("fitted_mass", fit_p4.mass());
        cand.addUserFloat("fitted_massErr",
                          sqrt(fitter_kk.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
  
        // other vars
        cand.addUserFloat("cos_theta_2D",
                          cos_theta_2D(fitter_kk, *beamspot, cand.p4()));
  
        cand.addUserFloat("fitted_cos_theta_2D",
                          cos_theta_2D(fitter_kk, *beamspot, fit_p4));
        //std::cout<< cos_theta_2D(fitter_kk, *beamspot, fit_p4) <<std::endl;
  
        auto lxy = l_xy(fitter_kk, *beamspot);
        cand.addUserFloat("l_xy", lxy.value());
        cand.addUserFloat("l_xy_unc", lxy.error());


	// lambda0 flys
	auto ld0_ldb_lxy = l_xy(fitter_kk, lambda_ptr->userFloat("vtx_x"), lambda_ptr->userFloat("vtx_y"), lambda_ptr->userFloat("vtx_z"));
        auto ld0_ldb_fitted_cos_theta_2D = cos_theta_2D(fitter_kk, lambda_ptr->userFloat("vtx_x"), lambda_ptr->userFloat("vtx_y"), lambda_ptr->userFloat("vtx_z"), fit_p4 - lambda_p4);
        cand.addUserFloat("ld0_ldb_fitted_cos_theta_2D", -ld0_ldb_fitted_cos_theta_2D);
        cand.addUserFloat("ld0_ldb_l_xy", ld0_ldb_lxy.value());
        cand.addUserFloat("ld0_ldb_l_xy_unc", ld0_ldb_lxy.error());

        // post fit selection
        if ( !post_vtx_selection_(cand) ) continue;
  
        cand.addUserFloat("vtx_x", cand.vx());
        cand.addUserFloat("vtx_y", cand.vy());
        cand.addUserFloat("vtx_z", cand.vz());
  
        const auto& covMatrix = fitter_kk.fitted_vtx_uncertainty();
        cand.addUserFloat("vtx_cxx", covMatrix.cxx());
        cand.addUserFloat("vtx_cyy", covMatrix.cyy());
        cand.addUserFloat("vtx_czz", covMatrix.czz());
        cand.addUserFloat("vtx_cyx", covMatrix.cyx());
        cand.addUserFloat("vtx_czx", covMatrix.czx());
        cand.addUserFloat("vtx_czy", covMatrix.czy());
  
        // refitted daughters (leptons/tracks)
        std::vector<std::string> dnames{ "trk3", "trk4", "lambda" };
  
        for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" , fitter_kk.daughter_p4(idaughter).pt() );
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta", fitter_kk.daughter_p4(idaughter).eta() );
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi", fitter_kk.daughter_p4(idaughter).phi() );
        }
  
        // track impact parameter from ditrack SV
        TrajectoryStateOnSurface tsos1 = extrapolator.extrapolate(transientTracks->at(trk1_idx).impactPointState(), fitter_kk.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP1 = signedTransverseImpactParameter(tsos1, fitter_kk.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk1_svip2d" , cur2DIP1.second.value());
        cand.addUserFloat("trk1_svip2d_err" , cur2DIP1.second.error());
  
        TrajectoryStateOnSurface tsos2 = extrapolator.extrapolate(transientTracks->at(trk2_idx).impactPointState(), fitter_kk.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP2 = signedTransverseImpactParameter(tsos2, fitter_kk.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk2_svip2d" , cur2DIP2.second.value());
        cand.addUserFloat("trk2_svip2d_err" , cur2DIP2.second.error());
  
        TrajectoryStateOnSurface tsos3 = extrapolator.extrapolate(transientTracks->at(trk3_idx).impactPointState(), fitter_kk.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP3 = signedTransverseImpactParameter(tsos3, fitter_kk.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk3_svip2d" , cur2DIP3.second.value());
        cand.addUserFloat("trk3_svip2d_err" , cur2DIP3.second.error());
  
        TrajectoryStateOnSurface tsos4 = extrapolator.extrapolate(transientTracks->at(trk4_idx).impactPointState(), fitter_kk.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP4 = signedTransverseImpactParameter(tsos4, fitter_kk.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk4_svip2d" , cur2DIP4.second.value());
        cand.addUserFloat("trk4_svip2d_err" , cur2DIP4.second.error());
  
  
        //compute isolation
        std::vector<float> isos = TrackerIsolation(tracks, cand, dnames );
        for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
          cand.addUserFloat(dnames[idaughter] + "_iso04", isos[idaughter]);
        }

	std::pair<float, int> isos_self = TrackerIsolation_self(tracks, cand);
        cand.addUserFloat("iso04", isos_self.first);
        cand.addUserFloat("iso04_ntrk", isos_self.second);

        ret_val->push_back(cand);

      } // for (size_t trk3_idx = 0; trk3_idx < tracks->size(); ++trk3_idx ) {
    } //for (size_t trk4_idx = trk3_idx + 1; trk4_idx < tracks->size(); ++trk4_idx) {
  } // for(size_t lambda_idx = 0; lambda_idx < lambda->size(); ++lambda_idx)

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LambdabToLambdahhBuilder);
