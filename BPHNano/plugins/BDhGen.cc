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

class BDhGen : public edm::global::EDProducer<> {

public:
  explicit BDhGen(const edm::ParameterSet &theParameters):
      genParticleToken_(consumes<reco::GenParticleCollection>(theParameters.getParameter<edm::InputTag>("genParticle")))
	{
          produces<pat::CompositeCandidateCollection>("Genmatch");
  }

  ~BDhGen() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
};

int GetGlobalIndex(const reco::Candidate *cand, const edm::Handle<reco::GenParticleCollection> &handle) {
    for (size_t j = 0; j < handle->size(); j++) {
        if (&(*handle)[j] == cand) return j;
    }
    return -1;
}


void BDhGen::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {

  using std::vector;
  edm::Handle<reco::GenParticleCollection> genParticleHandle_;
  iEvent.getByToken(genParticleToken_, genParticleHandle_);
  if (!genParticleHandle_.isValid()) {
    edm::LogError("BDhGen") << "GenParticles not found!";
    return;
  }

  std::unique_ptr<pat::CompositeCandidateCollection> Genmatch_out(new pat::CompositeCandidateCollection());

  for (size_t i = 0; i < genParticleHandle_->size(); i++) {

      const reco::Candidate *bCand = &(*genParticleHandle_)[i];
      int b_pdgId = bCand->pdgId();

      if (abs(b_pdgId) != 521) continue;  // Only B+/B-

      // Init candidate
      pat::CompositeCandidate Genmatch;
      Genmatch.setP4(bCand->p4());
      Genmatch.addUserInt("idx_b", GetGlobalIndex(bCand, genParticleHandle_));
      Genmatch.addUserInt("b_charge", bCand->charge());
      int channelFlag = -1;
      int idx_d0 = -1;          // D0 from B or D*0
      int idx_dstar0 = -1;      // D*0 from B
      int idx_dstar_decay = -1; // D*0 decay product, gamma or pi0
      int idx_b_pion = -1;         // pion from B
      int idx_b_kaon = -1;         // kaon from B
      int idx_ks = -1;
      int idx_pip = -1;
      int idx_pim = -1;
      int idx_pi0_1 = -1;
      int idx_pi0_2 = -1;
      int idx_kspip = -1;
      int idx_kspim = -1;
      int ndt_b = -1;
      int ndt_d0 = -1;
      int ndt_ks0 = -1;
      double ks0_flight_distance = -1;
      double ks0_flight_distance_2D = -1;

      const reco::Candidate *d0Cand = nullptr;

      // Loop over B daughters
      ndt_b = bCand->numberOfDaughters();
      for (size_t j = 0; j < bCand->numberOfDaughters(); j++) {
          const reco::Candidate *bDau = bCand->daughter(j);
          int pdgId_bDau = bDau->pdgId();

          // D0
          if (abs(pdgId_bDau) == 421) {
              idx_d0 = GetGlobalIndex(bDau, genParticleHandle_);
              d0Cand = bDau;
          }

          // D*0
          if (abs(pdgId_bDau) == 423) {
              idx_dstar0 = GetGlobalIndex(bDau, genParticleHandle_);

              // Loop over D*0 daughters
              for (size_t k = 0; k < bDau->numberOfDaughters(); k++) {
                  const reco::Candidate *dstarDau = bDau->daughter(k);
                  int pdgId_dstarDau = dstarDau->pdgId();

                  if (abs(pdgId_dstarDau) == 421) {
                      d0Cand = dstarDau;
                      idx_d0 = GetGlobalIndex(d0Cand, genParticleHandle_);
                  }

                  if (pdgId_dstarDau == 111 || pdgId_dstarDau == 22) {
                      idx_dstar_decay = GetGlobalIndex(dstarDau, genParticleHandle_);
                  }
              }
          }

          if (pdgId_bDau == 211 || pdgId_bDau == -211) {
              idx_b_pion = GetGlobalIndex(bDau, genParticleHandle_);
          }
          if (pdgId_bDau == 321 || pdgId_bDau == -321) {
              idx_b_kaon = GetGlobalIndex(bDau, genParticleHandle_);
          }
      }

      if (!d0Cand) continue;  // No D0 found, skip this B

      // Loop over D0 daughters
      // D0 could decay to Ks pi- pi+ pi0, also try to look for more pi0
      int found_ks = 0, found_pip = 0, found_pim = 0, found_pi0 = 0;
      int found_kspip = 0, found_kspim = 0;

      ndt_d0 = d0Cand->numberOfDaughters();
      for (size_t j = 0; j < d0Cand->numberOfDaughters(); j++) {
          const reco::Candidate *d0Dau = d0Cand->daughter(j);
          int pdgId_d0Dau = d0Dau->pdgId();

          if (pdgId_d0Dau == 310) {
      	reco::Candidate::Point ks0_production_vtx(d0Dau->vx(), d0Dau->vy(), d0Dau->vz());
      	if (d0Dau->numberOfDaughters() > 0) {
                  const reco::Candidate* dau0 = d0Dau->daughter(0);
                  reco::Candidate::Point ks0_decay_vtx(dau0->vx(), dau0->vy(), dau0->vz());
                  math::XYZVector ks0_flight = ks0_decay_vtx - ks0_production_vtx;
                  ks0_flight_distance = ks0_flight.R();
                  ks0_flight_distance_2D = sqrt(ks0_flight.x() * ks0_flight.x() + ks0_flight.y() * ks0_flight.y());
              }
              idx_ks = GetGlobalIndex(d0Dau, genParticleHandle_);
              found_ks++;
      	ndt_ks0 = d0Dau->numberOfDaughters();
      	for (size_t k = 0; k < d0Dau->numberOfDaughters(); k++) {
      	    const reco::Candidate *ks0Dau = d0Dau->daughter(k);
      	    int pdgId_ks0Dau = ks0Dau->pdgId();
      	    if (pdgId_ks0Dau == 211) {
                      idx_kspip = GetGlobalIndex(ks0Dau, genParticleHandle_);
                      found_kspip++;
                  }
                  if (pdgId_ks0Dau == -211) {
                      idx_kspim = GetGlobalIndex(ks0Dau, genParticleHandle_);
                      found_kspim++;
                  }
      	}
          }
          if (pdgId_d0Dau == 211) {
              idx_pip = GetGlobalIndex(d0Dau, genParticleHandle_);
              found_pip++;
          }
          if (pdgId_d0Dau == -211) {
              idx_pim = GetGlobalIndex(d0Dau, genParticleHandle_);
              found_pim++;
          }
          if (pdgId_d0Dau == 111) {
              if (found_pi0 == 0) idx_pi0_1 = GetGlobalIndex(d0Dau, genParticleHandle_);
              else if (found_pi0 == 1) idx_pi0_2 = GetGlobalIndex(d0Dau, genParticleHandle_);
              found_pi0++;
          }
      }

      // Define decay mode
      // 0, 1
      if(idx_b_pion!=-1){
          if (found_ks && found_pip && found_pim && found_pi0 == 0) {
              if (idx_dstar0 != -1) channelFlag = 1;  // B -> D*0 pi+, D*0 -> D0 pi0/γ, D0->  Ks pi+ pi-
              else channelFlag = 0;                   // B -> D0 pi+, D0 -> Ks pi+ pi-
          } else if (found_ks && found_pip && found_pim && found_pi0 == 1) {
              if (idx_dstar0 != -1) channelFlag = 3;  // B -> D*0 pi+, D*0 -> D0 pi0/γ, D0->  Ks pi+ pi- pi0
              else channelFlag = 2;                   // B -> D0 pi+, D0 -> Ks pi+ pi- pi0
          } else if (found_ks && found_pip && found_pim && found_pi0 == 2) {
              if (idx_dstar0 != -1) channelFlag = 5;  // B -> D*0 pi+, D*0 -> D0 pi0/γ, D0 ->  Ks pi+ pi- 2pi0
              else channelFlag = 4;                   // B -> D0 pi+, D0 -> Ks pi+ pi- 2pi0
          } else {
              continue;  // Skip if not matching any decay chain
          }
      }else if(idx_b_kaon!=-1){
          if (found_ks && found_pip && found_pim && found_pi0 == 0) {
              if (idx_dstar0 != -1) channelFlag = 7;  // B -> D*0 K+, D*0 -> D0 pi0/γ, D0->  Ks pi+ pi-
              else channelFlag = 6;                   // B -> D0 K+, D0 -> Ks pi+ pi-
          } else if (found_ks && found_pip && found_pim && found_pi0 == 1) {
              if (idx_dstar0 != -1) channelFlag = 9;  // B -> D*0 K+, D*0 -> D0 pi0/γ, D0->  Ks pi+ pi- pi0
              else channelFlag = 8;                   // B -> D0 K+, D0 -> Ks pi+ pi- pi0
          } else if (found_ks && found_pip && found_pim && found_pi0 == 2) {
              if (idx_dstar0 != -1) channelFlag = 11;  // B -> D*0 K+, D*0 -> D0 pi0/γ, D0 ->  Ks pi+ pi- 2pi0
              else channelFlag = 10;                   // B -> D0 K+, D0 -> Ks pi+ pi- 2pi0
          } else {
              continue;  // Skip if not matching any decay chain
          }
      }

      // Fill the CompositeCandidate
      Genmatch.addUserInt("channelFlag", channelFlag);
      Genmatch.addUserInt("b_numberOfDaughters", ndt_b);
      Genmatch.addUserInt("d0_numberOfDaughters", ndt_d0);
      Genmatch.addUserInt("ks0_numberOfDaughters", ndt_ks0);
      Genmatch.addUserFloat("ks0_flight_distance", ks0_flight_distance);
      Genmatch.addUserFloat("ks0_flight_distance_2D", ks0_flight_distance_2D);
      Genmatch.addUserInt("idx_b_pion", idx_b_pion);
      Genmatch.addUserInt("idx_b_kaon", idx_b_kaon);
      Genmatch.addUserInt("idx_dstar0", idx_dstar0);
      Genmatch.addUserInt("idx_dstar_decay", idx_dstar_decay);
      Genmatch.addUserInt("idx_d0", idx_d0);
      Genmatch.addUserInt("idx_ks", idx_ks);
      Genmatch.addUserInt("idx_kspip", idx_kspip);
      Genmatch.addUserInt("idx_kspim", idx_kspim);
      Genmatch.addUserInt("idx_pip", idx_pip);
      Genmatch.addUserInt("idx_pim", idx_pim);
      Genmatch.addUserInt("idx_pi0_1", idx_pi0_1);
      Genmatch.addUserInt("idx_pi0_2", idx_pi0_2);
      Genmatch_out->push_back(Genmatch);
  }

  iEvent.put(std::move(Genmatch_out), "Genmatch");

}

DEFINE_FWK_MODULE(BDhGen);

