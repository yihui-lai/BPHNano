////////////////////////////// EtaTo2L2PiBuilder //////////////////////////////
/// authors: Y Lai (Princeton)
// takes the ditrack collection and a dilepton collection and produce eta
// - ers using a four-track vertex

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

class EtaTo2L2PiBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit EtaTo2L2PiBuilder(const edm::ParameterSet &cfg):
    bFieldToken_{esConsumes<MagneticField, IdealMagneticFieldRecord>()},
    // selections
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    //inputs
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    tracks_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )}
  {
    //output
    produces<pat::CompositeCandidateCollection>();
  }

  ~EtaTo2L2PiBuilder() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_;
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_;

  // inputs
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracks_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;

};

void EtaTo2L2PiBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);


  edm::Handle<pat::CompositeCandidateCollection> tracks;
  evt.getByToken(tracks_, tracks);
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);

  edm::ESHandle<MagneticField> fieldHandle;
  const auto& bField = iSetup.getData(bFieldToken_);
  AnalyticalImpactPointExtrapolator extrapolator(&bField);
  
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());


  // all muon pairs
//  std::cout << "Number of diMuons: "<< dileptons->size()<<std::endl;

  for (size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
    edm::Ptr<pat::CompositeCandidate> ll_ptr(dileptons, ll_idx);
    edm::Ptr<reco::Candidate> l1_ptr = ll_ptr->userCand("l1");
    edm::Ptr<reco::Candidate> l2_ptr = ll_ptr->userCand("l2");
    int l1_idx = ll_ptr->userInt("l1_idx");
    int l2_idx = ll_ptr->userInt("l2_idx");

    // two pion tracks
    for (size_t trk1_idx = 0; trk1_idx < tracks->size(); ++trk1_idx ) {
      edm::Ptr<pat::CompositeCandidate> trk1_ptr( tracks, trk1_idx );
      if (!trk1_selection_(*trk1_ptr)) continue;

      for (size_t trk2_idx = trk1_idx + 1; trk2_idx < tracks->size(); ++trk2_idx) {
        edm::Ptr<pat::CompositeCandidate> trk2_ptr( tracks, trk2_idx );
        if (!trk2_selection_(*trk2_ptr)) continue;

        if ((trk1_ptr->charge() + trk2_ptr->charge())!=0) continue;

//	std::cout << "ditrack with opposite charge: "<<std::endl;

        // lambda_b candidate
        pat::CompositeCandidate cand;
        math::PtEtaPhiMLorentzVector l1_p4(
            l1_ptr->pt(),
            l1_ptr->eta(),
            l1_ptr->phi(),
            MUON_MASS
        );
        math::PtEtaPhiMLorentzVector l2_p4(
            l2_ptr->pt(),
            l2_ptr->eta(),
            l2_ptr->phi(),
            MUON_MASS
        );
        math::PtEtaPhiMLorentzVector trk1_p4(
            trk1_ptr->pt(),
            trk1_ptr->eta(),
            trk1_ptr->phi(),
            PI_MASS
        );
        math::PtEtaPhiMLorentzVector trk2_p4(
            trk2_ptr->pt(),
            trk2_ptr->eta(),
            trk2_ptr->phi(),
            PI_MASS
        );
	auto etap4 = l1_p4 + l2_p4 + trk1_p4 + trk2_p4;
	cand.setP4(etap4);
        cand.setCharge( l1_ptr->charge() + l2_ptr->charge() + trk1_ptr->charge() + trk2_ptr->charge() );
        // save daughters - unfitted
        cand.addUserCand("l1", l1_ptr);
        cand.addUserCand("l2", l2_ptr);
        cand.addUserCand("trk1", trk1_ptr);
        cand.addUserCand("trk2", trk2_ptr);
        cand.addUserCand("dilepton", ll_ptr);

        // save indices
        cand.addUserInt("l1_idx", l1_idx);
        cand.addUserInt("l2_idx", l2_idx);
        cand.addUserInt("ll_idx" , ll_idx);
        cand.addUserInt("trk1_idx", trk1_idx);
        cand.addUserInt("trk2_idx", trk2_idx);
        cand.addUserFloat("trk1_mass", trk1_ptr->mass());
        cand.addUserFloat("trk2_mass", trk2_ptr->mass());

        auto dr_info = min_max_dr({l1_ptr, l2_ptr, trk1_ptr, trk2_ptr});
        cand.addUserFloat("min_dr", dr_info.first);
        cand.addUserFloat("max_dr", dr_info.second);


        // check if pass pre vertex cut
        if ( !pre_vtx_selection_(cand) ) continue;
//	std::cout << "Pass presel: "<<std::endl;

        KinVtxFitter fitter(
          { leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), ttracks->at(trk1_idx), ttracks->at(trk2_idx) },
          { MUON_MASS, MUON_MASS, PI_MASS, PI_MASS},
          { LEP_SIGMA, LEP_SIGMA, PI_SIGMA, PI_SIGMA }
          );

        if (!fitter.success()) continue;
//	std::cout << "fitter.success "<<std::endl;

        // B0 position
        cand.setVertex(
          reco::Candidate::Point(
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )
        );

        // vertex vars
        cand.addUserFloat("sv_chi2", fitter.chi2());
        cand.addUserFloat("sv_ndof", fitter.dof());
        cand.addUserFloat("sv_prob", fitter.prob());

        // refitted kinematic vars
        cand.addUserFloat("fitted_ditrack_mass", 
                          (fitter.daughter_p4(2) + fitter.daughter_p4(3)).mass());
        cand.addUserFloat("fitted_mll",
                          (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());

        auto fit_p4 = fitter.fitted_p4();
        cand.addUserFloat("fitted_pt"  , fit_p4.pt());
        cand.addUserFloat("fitted_eta" , fit_p4.eta());
        cand.addUserFloat("fitted_phi" , fit_p4.phi());
        cand.addUserFloat("fitted_mass", fit_p4.mass());
        cand.addUserFloat("fitted_massErr",
                          sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));

        // other vars
        cand.addUserFloat("cos_theta_2D",
                          cos_theta_2D(fitter, *beamspot, cand.p4()));

        cand.addUserFloat("fitted_cos_theta_2D",
                          cos_theta_2D(fitter, *beamspot, fit_p4));
        //std::cout<< cos_theta_2D(fitter, *beamspot, fit_p4) <<std::endl;

        auto lxy = l_xy(fitter, *beamspot);
        cand.addUserFloat("l_xy", lxy.value());
        cand.addUserFloat("l_xy_unc", lxy.error());
                
        // post fit selection
        if ( !post_vtx_selection_(cand) ) continue;

        cand.addUserFloat("vtx_x", cand.vx());
        cand.addUserFloat("vtx_y", cand.vy());
        cand.addUserFloat("vtx_z", cand.vz());

        const auto& covMatrix = fitter.fitted_vtx_uncertainty();
        cand.addUserFloat("vtx_cxx", covMatrix.cxx());
        cand.addUserFloat("vtx_cyy", covMatrix.cyy());
        cand.addUserFloat("vtx_czz", covMatrix.czz());
        cand.addUserFloat("vtx_cyx", covMatrix.cyx());
        cand.addUserFloat("vtx_czx", covMatrix.czx());
        cand.addUserFloat("vtx_czy", covMatrix.czy());

        // refitted daughters (leptons/tracks)
        std::vector<std::string> dnames{ "l1", "l2", "trk1", "trk2" };

        for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" , fitter.daughter_p4(idaughter).pt() );
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta", fitter.daughter_p4(idaughter).eta() );
          cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi", fitter.daughter_p4(idaughter).phi() );
        }

        // track impact parameter from dilepton SV
        TrajectoryStateOnSurface tsos1 = extrapolator.extrapolate(ttracks->at(trk1_idx).impactPointState(), fitter.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP1 = signedTransverseImpactParameter(tsos1, fitter.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk1_svip2d" , cur2DIP1.second.value());
        cand.addUserFloat("trk1_svip2d_err" , cur2DIP1.second.error());

        TrajectoryStateOnSurface tsos2 = extrapolator.extrapolate(ttracks->at(trk2_idx).impactPointState(), fitter.fitted_vtx());
        std::pair<bool, Measurement1D> cur2DIP2 = signedTransverseImpactParameter(tsos2, fitter.fitted_refvtx(), *beamspot);
        cand.addUserFloat("trk2_svip2d" , cur2DIP2.second.value());
        cand.addUserFloat("trk2_svip2d_err" , cur2DIP2.second.error());

        //compute isolation
        std::vector<float> isos = TrackerIsolation(tracks, cand, dnames );
        for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
          cand.addUserFloat(dnames[idaughter] + "_iso04", isos[idaughter]);
        }

        ret_val->push_back(cand);

      } // for (size_t trk2_idx = trk2_idx + 1; trk2_idx < tracks->size(); ++trk2_idx) {
    } // for (size_t trk1_idx = 0; trk1_idx < tracks->size(); ++trk1_idx ) {
  } // for (size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {

  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EtaTo2L2PiBuilder);
