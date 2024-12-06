/////////////////////////////// DiTrackBuilder ///////////////////////////////
/// original author: G Karathanasis, CERN
// takes selected track collection and a mass hypothesis and produces ditrack ca
// -ndidates



#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"




class DiTrackBuilder : public edm::global::EDProducer<> {


public:

  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DiTrackBuilder(const edm::ParameterSet &cfg):
    trk1_selection_{cfg.getParameter<std::string>("trk1Selection")},
    trk2_selection_{cfg.getParameter<std::string>("trk2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    pfcands_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("tracks") )},
    ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracks") )},
    trk1_mass_{ cfg.getParameter<double>("trk1Mass")},
    trk2_mass_{ cfg.getParameter<double>("trk2Mass")}
  {

    //output
    produces<pat::CompositeCandidateCollection>();

  }

  ~DiTrackBuilder() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  const StringCutObjectSelector<pat::CompositeCandidate> trk1_selection_; // cuts on leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> trk2_selection_; // sub-leading cand
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> pfcands_; //input PF cands this is sorted in pT in previous step
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_; //input TTracks of PF cands
  double trk1_mass_;
  double trk2_mass_;
};


void DiTrackBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //inputs
  edm::Handle<pat::CompositeCandidateCollection> pfcands;
  evt.getByToken(pfcands_, pfcands);
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_, ttracks);


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> kstar_out(new pat::CompositeCandidateCollection());



  // main loop
  for (size_t trk1_idx = 0; trk1_idx < pfcands->size(); ++trk1_idx ) {

    edm::Ptr<pat::CompositeCandidate> trk1_ptr( pfcands, trk1_idx );
    if (!trk1_selection_(*trk1_ptr)) continue;

    for (size_t trk2_idx = trk1_idx + 1; trk2_idx < pfcands->size(); ++trk2_idx) {

      edm::Ptr<pat::CompositeCandidate> trk2_ptr( pfcands, trk2_idx );
      //if (trk1_ptr->charge() == trk2_ptr->charge()) continue;
      if (!trk2_selection_(*trk2_ptr)) continue;

      bool UsedAgain = false;
      // Loop in all possible hypothesis
      for ( std::pair<double, double> masses : { std::pair<double, double>(trk1_mass_, trk2_mass_), std::pair<double, double>(trk2_mass_, trk1_mass_) } ) {
        // create a K* candidate; add first quantities that can be used for pre fit selection
        pat::CompositeCandidate kstar_cand;
        auto trk1_p4 = trk1_ptr->polarP4();
        auto trk2_p4 = trk2_ptr->polarP4();
        trk1_p4.SetM(masses.first);
        trk2_p4.SetM(masses.second);
        //adding stuff for pre fit selection
        kstar_cand.setP4(trk1_p4 + trk2_p4);
        kstar_cand.setCharge(trk1_ptr->charge() + trk2_ptr->charge());
        kstar_cand.addUserFloat("trk_deltaR", reco::deltaR(*trk1_ptr, *trk2_ptr));
        // save indices
        kstar_cand.addUserInt("trk1_idx", trk1_idx );
        kstar_cand.addUserInt("trk2_idx", trk2_idx );
        kstar_cand.addUserFloat("trk1_mass", masses.first);
        kstar_cand.addUserFloat("trk2_mass", masses.second);
        // save cands
        kstar_cand.addUserCand("trk1", trk1_ptr );
        kstar_cand.addUserCand("trk2", trk2_ptr );
        // selection before fit
        if ( !pre_vtx_selection_(kstar_cand) ) continue;

        KinVtxFitter fitter(
              {ttracks->at(trk1_idx), ttracks->at(trk2_idx)},
            { masses.first, masses.second },
            {K_SIGMA, K_SIGMA} //K and PI sigma equal...
                           );

        if ( !fitter.success() ) continue;
        kstar_cand.setVertex(
          reco::Candidate::Point(
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          )
        );
        // save quantities after fit
        kstar_cand.addUserInt("sv_ok", fitter.success() ? 1 : 0);
        kstar_cand.addUserFloat("sv_chi2", fitter.chi2());
        kstar_cand.addUserFloat("sv_ndof", fitter.dof());
        kstar_cand.addUserFloat("sv_prob", fitter.prob());
        kstar_cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass() );
        kstar_cand.addUserFloat("fitted_pt",
                                fitter.fitted_candidate().globalMomentum().perp() );

        kstar_cand.addUserFloat("fitted_eta",
                                fitter.fitted_candidate().globalMomentum().eta() );

        kstar_cand.addUserFloat("fitted_phi",
                                fitter.fitted_candidate().globalMomentum().phi() );

        kstar_cand.addUserInt("second_mass_hypothesis", UsedAgain );
        kstar_cand.addUserFloat("vtx_x", kstar_cand.vx());
        kstar_cand.addUserFloat("vtx_y", kstar_cand.vy());
        kstar_cand.addUserFloat("vtx_z", kstar_cand.vz());

        // after fit selection
        if ( !post_vtx_selection_(kstar_cand) ) continue;
        kstar_out->emplace_back(kstar_cand);
        UsedAgain = true;
        if (masses.first == masses.second) break;

      } // end for ( auto & masses:
    } // end for(size_t trk2_idx = trk1_idx + 1
  } //for(size_t trk1_idx = 0

  evt.put(std::move(kstar_out));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DiTrackBuilder);