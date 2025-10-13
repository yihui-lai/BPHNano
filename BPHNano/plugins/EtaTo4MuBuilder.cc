/// original authors: RK18 team
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

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
#include <numeric>
#include "KinVtxFitter.h"

template<typename Lepton>
class EtaTo4LepBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<Lepton> LeptonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit EtaTo4LepBuilder(const edm::ParameterSet &cfg):
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    l3_selection_{cfg.getParameter<std::string>("lep3Selection")},
    l4_selection_{cfg.getParameter<std::string>("lep4Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<LeptonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )} {
    produces<pat::CompositeCandidateCollection>("Selected4Leptons");
  }

  ~EtaTo4LepBuilder() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}

private:
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<Lepton> l3_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l4_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
  const edm::EDGetTokenT<LeptonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
};

template<typename Lepton>
void EtaTo4LepBuilder<Lepton>::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<LeptonCollection> leptons;
  evt.getByToken(src_, leptons);

  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  std::unique_ptr<std::vector<KinVtxFitter> > kinVtx_out( new std::vector<KinVtxFitter> );

  // std::cout<<"leptons->size: "<< leptons->size()<<std::endl;
  const size_t nLep = leptons->size();
  if (nLep < 4) {
    evt.put(std::move(ret_value), "Selected4Leptons");
    return;
  }

  for (size_t l1_idx = 0; l1_idx < nLep; ++l1_idx) {
    edm::Ptr<Lepton> l1_ptr(leptons, l1_idx);
    if (!l1_selection_(*l1_ptr)) continue;

    for (size_t l2_idx = l1_idx + 1; l2_idx < nLep; ++l2_idx) {
      edm::Ptr<Lepton> l2_ptr(leptons, l2_idx);
      if (!l2_selection_(*l2_ptr)) continue;
      
      for (size_t l3_idx = l2_idx + 1; l3_idx < nLep; ++l3_idx) {
        edm::Ptr<Lepton> l3_ptr(leptons, l3_idx);
        if (!l3_selection_(*l3_ptr)) continue;

        for (size_t l4_idx = l3_idx + 1; l4_idx < nLep; ++l4_idx) {
          edm::Ptr<Lepton> l4_ptr(leptons, l4_idx);
          if (!l4_selection_(*l4_ptr)) continue;

          pat::CompositeCandidate four_lepton;
          four_lepton.setP4(l1_ptr->p4() + l2_ptr->p4() + l3_ptr->p4() + l4_ptr->p4());
          four_lepton.setCharge(l1_ptr->charge() + l2_ptr->charge() + l3_ptr->charge() + l4_ptr->charge());

          four_lepton.addUserInt("l1_idx", l1_idx);
          four_lepton.addUserInt("l2_idx", l2_idx);
          four_lepton.addUserInt("l3_idx", l3_idx);
          four_lepton.addUserInt("l4_idx", l4_idx);

          four_lepton.addUserCand("l1", l1_ptr);
          four_lepton.addUserCand("l2", l2_ptr);
          four_lepton.addUserCand("l3", l3_ptr);
          four_lepton.addUserCand("l4", l4_ptr);

          if (!pre_vtx_selection_(four_lepton)) continue;

	  // Vertex fit
	  KinVtxFitter fitter(
          {ttracks->at(l1_idx), ttracks->at(l2_idx), ttracks->at(l3_idx), ttracks->at(l4_idx)},
          {l1_ptr->mass(), l2_ptr->mass(), l3_ptr->mass(), l4_ptr->mass()},
          {LEP_SIGMA, LEP_SIGMA, LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
          );
          if (!fitter.success()) continue;
          
	  four_lepton.setVertex(reco::Candidate::Point(
            fitter.fitted_vtx().x(),
            fitter.fitted_vtx().y(),
            fitter.fitted_vtx().z()
          ));

          std::vector<std::string> dnames{ "l1", "l2", "l3", "l4" };
          for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
            four_lepton.addUserFloat("fitted_" + dnames[idaughter] + "_pt" , fitter.daughter_p4(idaughter).pt() );
            four_lepton.addUserFloat("fitted_" + dnames[idaughter] + "_eta", fitter.daughter_p4(idaughter).eta() );
            four_lepton.addUserFloat("fitted_" + dnames[idaughter] + "_phi", fitter.daughter_p4(idaughter).phi() );
          }

          //four_lepton.addUserFloat("lep_deltaR", reco::deltaR(*l1_ptr, *l2_ptr));
          four_lepton.addUserInt("sv_ok", fitter.success() ? 1 : 0);
          four_lepton.addUserFloat("sv_chi2", fitter.chi2());
          four_lepton.addUserFloat("sv_ndof", fitter.dof());
          four_lepton.addUserFloat("sv_prob", fitter.prob());
          four_lepton.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
          four_lepton.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)) : -1);
          four_lepton.addUserFloat("vtx_x", four_lepton.vx());
          four_lepton.addUserFloat("vtx_y", four_lepton.vy());
          four_lepton.addUserFloat("vtx_z", four_lepton.vz());

          if (!post_vtx_selection_(four_lepton)) continue;

	  // You can compute pairwise dca if needed, e.g. average over 6 pairs:
          float dca_sum = 0; int dca_n = 0;
          std::array<size_t, 4> idxs = {{l1_idx, l2_idx, l3_idx, l4_idx}};
          for (size_t i = 0; i < 4; ++i) {
            for (size_t j = i+1; j < 4; ++j) {
              const auto& imp1 = ttracks->at(idxs[i]).impactPointTSCP();
              const auto& imp2 = ttracks->at(idxs[j]).impactPointTSCP();
              if (!imp1.isValid() || !imp2.isValid()) continue;
              ClosestApproachInRPhi cApp;
              cApp.calculate(imp1.theState(), imp2.theState());
              if (!cApp.status()) continue;
              dca_sum += std::abs(cApp.distance());
              dca_n++;
            }
          }
          four_lepton.addUserFloat("dca_avg", dca_n > 0 ? dca_sum/dca_n : -1);

          std::vector<float> dRs;
          dRs.reserve(6);
          dRs.push_back(reco::deltaR(*l1_ptr, *l2_ptr));
          dRs.push_back(reco::deltaR(*l1_ptr, *l3_ptr));
          dRs.push_back(reco::deltaR(*l1_ptr, *l4_ptr));
          dRs.push_back(reco::deltaR(*l2_ptr, *l3_ptr));
          dRs.push_back(reco::deltaR(*l2_ptr, *l4_ptr));
          dRs.push_back(reco::deltaR(*l3_ptr, *l4_ptr));
          
          float minDR = *std::min_element(dRs.begin(), dRs.end());
          float maxDR = *std::max_element(dRs.begin(), dRs.end());
          float avgDR = std::accumulate(dRs.begin(), dRs.end(), 0.0f) / dRs.size();
          
          four_lepton.addUserFloat("lep_min_deltaR", minDR);
          four_lepton.addUserFloat("lep_max_deltaR", maxDR);
          four_lepton.addUserFloat("lep_avg_deltaR", avgDR);

          ret_value->push_back(four_lepton);
	}
      }
    }
  }

  evt.put(std::move(ret_value), "Selected4Leptons");
}

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
typedef EtaTo4LepBuilder<pat::Muon> EtaTo4MuBuilder;
typedef EtaTo4LepBuilder<pat::Electron> EtaTo4ElBuilder;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EtaTo4MuBuilder);
DEFINE_FWK_MODULE(EtaTo4ElBuilder);
