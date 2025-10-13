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

#include "helper.h"
#include <memory>
#include <tuple>  // For std::tuple

class BDKstarGen : public edm::global::EDProducer<> {

public:
  explicit BDKstarGen(const edm::ParameterSet &theParameters):
      genParticleToken_(consumes<reco::GenParticleCollection>(theParameters.getParameter<edm::InputTag>("genParticle")))
	{
          produces<pat::CompositeCandidateCollection>("BDKstarGenmatch");
  }

  ~BDKstarGen() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
};


int GetBDKstarGlobalIndex(const reco::Candidate *cand, const edm::Handle<reco::GenParticleCollection> &handle) {
    for (size_t j = 0; j < handle->size(); j++) {
        if (&(*handle)[j] == cand) return j;
    }
    return -1;
}

void BDKstarGen::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {


  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);
  if (!genParticles.isValid()) {
    edm::LogError("BplusToD0KstarGen") << "GenParticles not found!";
    return;
  }

  std::unique_ptr<pat::CompositeCandidateCollection> out(new pat::CompositeCandidateCollection());


  for (size_t i = 0; i < genParticles->size(); i++) {

    const reco::Candidate *bCand = &(*genParticles)[i];
    if (abs(bCand->pdgId()) != 521) continue; // B+ or B-

    pat::CompositeCandidate Bmatch;
    Bmatch.setP4(bCand->p4());
    Bmatch.addUserInt("idx_B", GetBDKstarGlobalIndex(bCand, genParticles));
    Bmatch.addUserInt("B_charge", bCand->charge());

    int idx_D0 = -1;
    int idx_Kstar = -1;
    int idx_Ks = -1;
    int idx_D0_dau1 = -1;
    int idx_D0_dau2 = -1;
    int idx_Ks_pi1 = -1;
    int idx_Ks_pi2 = -1;
    int idx_Kstar_pi = -1;

    const reco::Candidate *D0Cand = nullptr;
    const reco::Candidate *KstarCand = nullptr;
    const reco::Candidate *KsCand = nullptr;

    // Loop over B+ daughters
    for (size_t j=0; j < bCand->numberOfDaughters(); j++) {
      const reco::Candidate *dau = bCand->daughter(j);
      int dauId = abs(dau->pdgId());

      if (dauId == 421) { // D0
        idx_D0 = GetBDKstarGlobalIndex(dau, genParticles);
        D0Cand = dau;
      } else if (dauId == 323) { // K*+
        idx_Kstar = GetBDKstarGlobalIndex(dau, genParticles);
        KstarCand = dau;
      }
    }

    if (!D0Cand || !KstarCand) continue; // skip if missing

    // D0 daughters
    if (D0Cand->numberOfDaughters() != 2) continue;
    int pdg0 = abs(D0Cand->daughter(0)->pdgId());
    int pdg1 = abs(D0Cand->daughter(1)->pdgId());
    // Allowed channels: K pi, K K, pi pi, pi K
    if (!((pdg0==321 && pdg1==211) || (pdg0==321 && pdg1==321) || (pdg0==211 && pdg1==211) || (pdg0==211 && pdg1==321))) continue;

    idx_D0_dau1 = GetBDKstarGlobalIndex(D0Cand->daughter(0), genParticles);
    idx_D0_dau2 = GetBDKstarGlobalIndex(D0Cand->daughter(1), genParticles);

    // K*+ daughters: Ks0 + pi+
    for (size_t j=0; j < KstarCand->numberOfDaughters(); j++) {
      const reco::Candidate *kstarDau = KstarCand->daughter(j);
      int dauId = abs(kstarDau->pdgId());
      if (dauId == 310) { // Ks0
        idx_Ks = GetBDKstarGlobalIndex(kstarDau, genParticles);
        KsCand = kstarDau;
      } else if (dauId == 211) { // pi+
        idx_Kstar_pi = GetBDKstarGlobalIndex(kstarDau, genParticles);
      }
    }
    if (!KsCand || idx_Kstar_pi==-1) continue;

    // Ks0 daughters: pi+ pi-
    if (KsCand->numberOfDaughters()!=2) continue;
    idx_Ks_pi1 = GetBDKstarGlobalIndex(KsCand->daughter(0), genParticles);
    idx_Ks_pi2 = GetBDKstarGlobalIndex(KsCand->daughter(1), genParticles);

    // Fill candidate
    Bmatch.addUserInt("idx_D0", idx_D0);
    Bmatch.addUserInt("idx_D0_dau1", idx_D0_dau1);
    Bmatch.addUserInt("idx_D0_dau2", idx_D0_dau2);
    Bmatch.addUserInt("idx_Kstar", idx_Kstar);
    Bmatch.addUserInt("idx_Kstar_pi", idx_Kstar_pi);
    Bmatch.addUserInt("idx_Ks", idx_Ks);
    Bmatch.addUserInt("idx_Ks_pi1", idx_Ks_pi1);
    Bmatch.addUserInt("idx_Ks_pi2", idx_Ks_pi2);

    out->push_back(Bmatch);
  }

  iEvent.put(std::move(out), "BDKstarGenmatch"); // Changed output label

}

DEFINE_FWK_MODULE(BDKstarGen);

