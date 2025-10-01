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

class Lambdab0Gen : public edm::global::EDProducer<> {

public:
  explicit Lambdab0Gen(const edm::ParameterSet &theParameters):
      genParticleToken_(consumes<reco::GenParticleCollection>(theParameters.getParameter<edm::InputTag>("genParticle")))
	{
          produces<pat::CompositeCandidateCollection>("LambdaBGenmatch");
  }

  ~Lambdab0Gen() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
};


int GetLambdab0GlobalIndex(const reco::Candidate *cand, const edm::Handle<reco::GenParticleCollection> &handle) {
    for (size_t j = 0; j < handle->size(); j++) {
        if (&(*handle)[j] == cand) return j;
    }
    return -1;
}

void Lambdab0Gen::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {

  using std::vector;
  edm::Handle<reco::GenParticleCollection> genParticleHandle_;
  iEvent.getByToken(genParticleToken_, genParticleHandle_);
  if (!genParticleHandle_.isValid()) {
    edm::LogError("Lambdab0Gen") << "GenParticles not found!";
    return;
  }

  std::unique_ptr<pat::CompositeCandidateCollection> Genmatch_out(new pat::CompositeCandidateCollection());

  for (size_t i = 0; i < genParticleHandle_->size(); i++) {

    const reco::Candidate *lambdaBCand = &(*genParticleHandle_)[i];
    int lambdaB_pdgId = lambdaBCand->pdgId();

    // Only Lambda_b0 (PDG ID 5122 or -5122)
    if (abs(lambdaB_pdgId) != 5122) continue;

    // Init candidate
    pat::CompositeCandidate Genmatch;
    Genmatch.setP4(lambdaBCand->p4());
    Genmatch.addUserInt("idx_lambdaB", GetLambdab0GlobalIndex(lambdaBCand, genParticleHandle_));
    Genmatch.addUserInt("lambdaB_charge", lambdaBCand->charge()); // Lambda_b0 is neutral, so charge will be 0

    int channelFlag = -1;
    int idx_lambda0 = -1;
    int idx_b_kaon1 = -1;
    int idx_b_kaon2 = -1;
    int idx_b_pion1 = -1;
    int idx_b_pion2 = -1;
    int idx_proton = -1;
    int idx_pion_from_lambda0 = -1;
    int ndt_lambdaB = -1;
    int ndt_lambda0 = -1;

    const reco::Candidate *lambda0Cand = nullptr;

    // Loop over Lambda_b0 daughters
    ndt_lambdaB = lambdaBCand->numberOfDaughters();
    for (size_t j = 0; j < lambdaBCand->numberOfDaughters(); j++) {
      const reco::Candidate *lambdaBDau = lambdaBCand->daughter(j);
      int pdgId_lambdaBDau = lambdaBDau->pdgId();

      // Lambda_0 (PDG ID 3122 or -3122)
      if (abs(pdgId_lambdaBDau) == 3122) {
        idx_lambda0 = GetLambdab0GlobalIndex(lambdaBDau, genParticleHandle_);
        lambda0Cand = lambdaBDau;
      }
      // Kaon (K+ or K-)
      if (abs(pdgId_lambdaBDau) == 321) {
        if (idx_b_kaon1 == -1) {
          idx_b_kaon1 = GetLambdab0GlobalIndex(lambdaBDau, genParticleHandle_);
        } else {
          idx_b_kaon2 = GetLambdab0GlobalIndex(lambdaBDau, genParticleHandle_);
        }
      }
      // Pion (pi+ or pi-)
      if (abs(pdgId_lambdaBDau) == 211) {
        if (idx_b_pion1 == -1) {
          idx_b_pion1 = GetLambdab0GlobalIndex(lambdaBDau, genParticleHandle_);
        } else {
          idx_b_pion2 = GetLambdab0GlobalIndex(lambdaBDau, genParticleHandle_);
        }
      }
    }

    if (!lambda0Cand) continue; // No Lambda_0 found, skip this Lambda_b0

    // Loop over Lambda_0 daughters
    int found_proton = 0, found_pion_from_lambda0 = 0;
    ndt_lambda0 = lambda0Cand->numberOfDaughters();
    for (size_t j = 0; j < lambda0Cand->numberOfDaughters(); j++) {
      const reco::Candidate *lambda0Dau = lambda0Cand->daughter(j);
      int pdgId_lambda0Dau = lambda0Dau->pdgId();

      // Proton (PDG ID 2212 or -2212)
      if (abs(pdgId_lambda0Dau) == 2212) {
        idx_proton = GetLambdab0GlobalIndex(lambda0Dau, genParticleHandle_);
        found_proton++;
      }
      // Pion (pi+ or pi-)
      if (abs(pdgId_lambda0Dau) == 211) {
        idx_pion_from_lambda0 = GetLambdab0GlobalIndex(lambda0Dau, genParticleHandle_);
        found_pion_from_lambda0++;
      }
    }

    // Check if Lambda_0 decayed to proton and pion
    if (!(found_proton == 1 && found_pion_from_lambda0 == 1)) {
        continue; // Lambda_0 did not decay as expected
    }

    // Define decay mode based on the remaining Lambda_b0 daughters
    int num_kaons = 0;
    if (idx_b_kaon1 != -1) num_kaons++;
    if (idx_b_kaon2 != -1) num_kaons++;

    int num_pions = 0;
    if (idx_b_pion1 != -1) num_pions++;
    if (idx_b_pion2 != -1) num_pions++;

    // Lambda_b0 -> Lambda_0 + Kaon + Kaon
    if (num_kaons == 2 && num_pions == 0) {
      channelFlag = 0;
    }
    // Lambda_b0 -> Lambda_0 + Kaon + Pion
    else if (num_kaons == 1 && num_pions == 1) {
      channelFlag = 1;
    }
    // Lambda_b0 -> Lambda_0 + Pion + Pion
    else if (num_kaons == 0 && num_pions == 2) {
      channelFlag = 2;
    }
    else {
      continue; // Skip if not matching any of the desired decay chains
    }


    // Fill the CompositeCandidate
    Genmatch.addUserInt("channelFlag", channelFlag);
    Genmatch.addUserInt("lambdaB_numberOfDaughters", ndt_lambdaB);
    Genmatch.addUserInt("lambda0_numberOfDaughters", ndt_lambda0);
    Genmatch.addUserInt("idx_lambda0", idx_lambda0);
    Genmatch.addUserInt("idx_b_kaon1", idx_b_kaon1);
    Genmatch.addUserInt("idx_b_kaon2", idx_b_kaon2);
    Genmatch.addUserInt("idx_b_pion1", idx_b_pion1);
    Genmatch.addUserInt("idx_b_pion2", idx_b_pion2);
    Genmatch.addUserInt("idx_proton", idx_proton);
    Genmatch.addUserInt("idx_pion_from_lambda0", idx_pion_from_lambda0);

    Genmatch_out->push_back(Genmatch);
  }

  iEvent.put(std::move(Genmatch_out), "LambdaBGenmatch"); // Changed output label

}

DEFINE_FWK_MODULE(Lambdab0Gen);

