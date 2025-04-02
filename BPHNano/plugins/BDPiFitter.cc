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
#include <tuple>  // For std::tuple

// pdg mass constants
namespace {
  const double piMass = 0.13957018;
  const double piMassSquared = piMass * piMass;
  const double kShortMass = 0.497614;
  const double kShortMassSquared = kShortMass * kShortMass;
  const double D0Mass = 1.86484;
  const float pion_sigma = piMass*1.e-6;

}  // namespace

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

int FindGlobalIndex(const reco::Candidate *cand, const edm::Handle<reco::GenParticleCollection> &handle) {
    for (size_t j = 0; j < handle->size(); j++) {
        if (&(*handle)[j] == cand) return j;
    }
    return -1;
}

class BDPiFitter : public edm::global::EDProducer<> {

public:
  explicit BDPiFitter(const edm::ParameterSet &theParameters):
      bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      token_beamSpot(consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"))),
      token_vertices(consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("tracks"))),
      lostTracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("lostTracks"))),
      genParticleToken_(consumes<reco::GenParticleCollection>(theParameters.getParameter<edm::InputTag>("genParticle"))),
      tkChi2Cut_( theParameters.getParameter<double>("tkChi2Cut")),
      tkNHitsCut_( theParameters.getParameter<int>("tkNHitsCut")),
      tkPtCut_( theParameters.getParameter<double>("tkPtCut")),
      tkIPSigXYCut_( theParameters.getParameter<double>("tkIPSigXYCut")),
      tkIPSigZCut_( theParameters.getParameter<double>("tkIPSigZCut")),
      vtxChi2Cut_( theParameters.getParameter<double>("vtxChi2Cut")),
      vtxDecaySigXYZCut_( theParameters.getParameter<double>("vtxDecaySigXYZCut")),
      vtxDecaySigXYCut_( theParameters.getParameter<double>("vtxDecaySigXYCut")),
      vtxDecayXYCut_( theParameters.getParameter<double>("vtxDecayXYCut")),
      ssVtxDecayXYCut_( theParameters.getParameter<double>("ssVtxDecayXYCut")),
      allowSS_( theParameters.getParameter<bool>("allowSS")),
      innerOuterTkDCAThreshold_( theParameters.getParameter<double>("innerOuterTkDCAThreshold")),
      innerTkDCACut_( theParameters.getParameter<double>("innerTkDCACut")),
      outerTkDCACut_( theParameters.getParameter<double>("outerTkDCACut")),
      allowWideAngleVtx_( theParameters.getParameter<bool>("allowWideAngleVtx")),
      mPiPiCut_( theParameters.getParameter<double>("mPiPiCut")),
      innerHitPosCut_( theParameters.getParameter<double>("innerHitPosCut")),
      cosThetaXYCut_( theParameters.getParameter<double>("cosThetaXYCut")),
      cosThetaXYZCut_( theParameters.getParameter<double>("cosThetaXYZCut")),
      kShortMassCut_( theParameters.getParameter<double>("kShortMassCut"))
  {
          produces<pat::CompositeCandidateCollection>("Genmatch");
      	  produces<pat::CompositeCandidateCollection>("SelectedTracks");
	  produces<reco::VertexCompositePtrCandidateCollection>("LooseKshort");
          produces<pat::CompositeCandidateCollection>("SelectedV0Collection");
          produces<TransientTrackCollection>("SelectedV0TransientCollection");
	  produces<pat::CompositeCandidateCollection>("D0");
  }

  ~BDPiFitter() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTracksToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;

  const double tkChi2Cut_;
  const int tkNHitsCut_;
  const double tkPtCut_;
  const double tkIPSigXYCut_;
  const double tkIPSigZCut_;
  const double vtxChi2Cut_;
  const double vtxDecaySigXYZCut_;
  const double vtxDecaySigXYCut_;
  const double vtxDecayXYCut_;
  const double ssVtxDecayXYCut_;
  const bool allowSS_;
  const double innerOuterTkDCAThreshold_;
  const double innerTkDCACut_;
  const double outerTkDCACut_;
  const bool allowWideAngleVtx_;
  const double mPiPiCut_;
  const double innerHitPosCut_;
  const double cosThetaXYCut_;
  const double cosThetaXYZCut_;
  const double kShortMassCut_;
  const bool isSignal=true;
  const bool onlymatched=true;
  const bool onlygen=false;
};

std::tuple<int, float, float, float> computeIsoAndMatch(
    const reco::Candidate* genCand,
    std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vectrk_ttrk)
{
    int matchedIndex = -1;
    float thematch_pt = 0;
    float thematch_dr = 99;
    float thematch_iso = 0;
    if (genCand->pt()==0) return std::make_tuple(matchedIndex, thematch_pt, thematch_iso, thematch_dr);
    for (unsigned int iSort = 0; iSort < vectrk_ttrk.size(); ++iSort) {
        const auto& track = vectrk_ttrk[iSort].first;
        float dr = deltaR(genCand->eta(), genCand->phi(), track.eta(), track.phi());
        if (dr > 0 && dr < 0.3) {
            thematch_iso += track.pt();
            float relPtDiff = std::abs(track.pt() - genCand->pt()) / genCand->pt();
            if (relPtDiff < 0.5 && dr < thematch_dr) {
                thematch_dr = dr;
		thematch_pt = track.pt();
                matchedIndex = int(iSort);
            }
        }
    }

    return std::make_tuple(matchedIndex, thematch_pt, thematch_iso, thematch_dr);
}

std::vector<int> FindBDecays(const edm::Handle<reco::GenParticleCollection> &genParticleHandle_,
                 std::vector<pat::CompositeCandidate> *Genmatch_out,
                 std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vectrk_ttrk
		) {
    
    std::vector<int> genpion_index;
    for (size_t i = 0; i < genParticleHandle_->size(); i++) {

        const reco::Candidate *bCand = &(*genParticleHandle_)[i];
        int b_pdgId = bCand->pdgId();

        if (abs(b_pdgId) != 521) continue;  // Only B+/B-

        // Init candidate
        pat::CompositeCandidate Genmatch;
        Genmatch.setP4(bCand->p4());
        Genmatch.addUserInt("idx_b", FindGlobalIndex(bCand, genParticleHandle_));
        Genmatch.addUserInt("b_charge", bCand->charge());
        int channelFlag = -1;
        int idx_bpion = -1;
        int idx_d0 = -1;
        int idx_dstar0 = -1;
        int idx_dstar_decay = -1;
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
	double Ks_flight_distance = -1;
        double Ks_flight_distance_2D = -1;

        const reco::Candidate *d0Cand = nullptr;

        // Loop over B daughters
	ndt_b = bCand->numberOfDaughters();
        for (size_t j = 0; j < bCand->numberOfDaughters(); j++) {
            const reco::Candidate *bDau = bCand->daughter(j);
            int pdgId_bDau = bDau->pdgId();

            // D0
            if (abs(pdgId_bDau) == 421) {
                idx_d0 = FindGlobalIndex(bDau, genParticleHandle_);
                d0Cand = bDau;
            }

            // D*0
            if (abs(pdgId_bDau) == 423) {
                idx_dstar0 = FindGlobalIndex(bDau, genParticleHandle_);

                // Loop over D*0 daughters
                for (size_t k = 0; k < bDau->numberOfDaughters(); k++) {
                    const reco::Candidate *dstarDau = bDau->daughter(k);
                    int pdgId_dstarDau = dstarDau->pdgId();

                    if (abs(pdgId_dstarDau) == 421) {
                        d0Cand = dstarDau;
                        idx_d0 = FindGlobalIndex(d0Cand, genParticleHandle_);
                    }

                    if (pdgId_dstarDau == 111 || pdgId_dstarDau == 22) {
                        idx_dstar_decay = FindGlobalIndex(dstarDau, genParticleHandle_);
                    }
                }
            }

            if (pdgId_bDau == 211 || pdgId_bDau == -211) {
                idx_bpion = FindGlobalIndex(bDau, genParticleHandle_);
            }
        }

        if (!d0Cand) continue;  // No D0 found, skip this B

        // Loop over D0 daughters
        int found_ks = 0, found_pip = 0, found_pim = 0, found_pi0 = 0;
	int found_kspip = 0, found_kspim = 0;
        ndt_d0 = d0Cand->numberOfDaughters();
        for (size_t j = 0; j < d0Cand->numberOfDaughters(); j++) {
            const reco::Candidate *d0Dau = d0Cand->daughter(j);
            int pdgId_d0Dau = d0Dau->pdgId();

            if (pdgId_d0Dau == 310) {
		reco::Candidate::Point Ks_production_vtx(d0Dau->vx(), d0Dau->vy(), d0Dau->vz());
		if (d0Dau->numberOfDaughters() > 0) {
                    const reco::Candidate* dau0 = d0Dau->daughter(0);
                    reco::Candidate::Point Ks_decay_vtx(dau0->vx(), dau0->vy(), dau0->vz());
                    math::XYZVector Ks_flight = Ks_decay_vtx - Ks_production_vtx;
                    Ks_flight_distance = Ks_flight.R();
                    Ks_flight_distance_2D = sqrt(Ks_flight.x() * Ks_flight.x() + Ks_flight.y() * Ks_flight.y());
                }
                idx_ks = FindGlobalIndex(d0Dau, genParticleHandle_);
                found_ks++;
		ndt_ks0 = d0Dau->numberOfDaughters();
		for (size_t k = 0; k < d0Dau->numberOfDaughters(); k++) {
  		    const reco::Candidate *ks0Dau = d0Dau->daughter(k);
		    int pdgId_ks0Dau = ks0Dau->pdgId();
		    if (pdgId_ks0Dau == 211) {
                        idx_kspip = FindGlobalIndex(ks0Dau, genParticleHandle_);
                        found_kspip++;
                    }
                    if (pdgId_ks0Dau == -211) {
                        idx_kspim = FindGlobalIndex(ks0Dau, genParticleHandle_);
                        found_kspim++;
                    }
		}
            }
            if (pdgId_d0Dau == 211) {
                idx_pip = FindGlobalIndex(d0Dau, genParticleHandle_);
                found_pip++;
            }
            if (pdgId_d0Dau == -211) {
                idx_pim = FindGlobalIndex(d0Dau, genParticleHandle_);
                found_pim++;
            }
            if (pdgId_d0Dau == 111) {
                if (found_pi0 == 0) idx_pi0_1 = FindGlobalIndex(d0Dau, genParticleHandle_);
                else if (found_pi0 == 1) idx_pi0_2 = FindGlobalIndex(d0Dau, genParticleHandle_);
                found_pi0++;
            }
        }

        // Define decay mode
	// 0, 1
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

        // Fill the CompositeCandidate
        Genmatch.addUserInt("channelFlag", channelFlag);
        Genmatch.addUserInt("b_numberOfDaughters", ndt_b);
        Genmatch.addUserInt("d0_numberOfDaughters", ndt_d0);
        Genmatch.addUserInt("ks0_numberOfDaughters", ndt_ks0);
        Genmatch.addUserFloat("Ks_flight_distance", Ks_flight_distance);
        Genmatch.addUserFloat("Ks_flight_distance_2D", Ks_flight_distance_2D);
	Genmatch.addUserInt("idx_bpion", idx_bpion);
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

        int idx_kspip_thematch = -1, idx_kspim_thematch = -1, idx_pip_thematch = -1, idx_pim_thematch = -1;
        float idx_kspip_pt = -1, idx_kspim_pt = -1, idx_pip_pt = -1, idx_pim_pt = -1;
        float idx_kspip_iso03 = -1, idx_kspim_iso03 = -1, idx_pip_iso03 = -1, idx_pim_iso03 = -1;
        float idx_kspip_dr = -1, idx_kspim_dr = -1, idx_pip_dr = -1, idx_pim_dr = -1;

        if (idx_kspip >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_kspip];
            std::tie(idx_kspip_thematch, idx_kspip_pt, idx_kspip_iso03, idx_kspip_dr) = computeIsoAndMatch(genCand, vectrk_ttrk);
        }
        Genmatch.addUserFloat("idx_kspip_pt", idx_kspip_pt);
        Genmatch.addUserFloat("idx_kspip_iso03", idx_kspip_iso03);
        Genmatch.addUserFloat("idx_kspip_dr", idx_kspip_dr);
        if (idx_kspim >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_kspim];
            std::tie(idx_kspim_thematch, idx_kspim_pt, idx_kspim_iso03, idx_kspim_dr) = computeIsoAndMatch(genCand, vectrk_ttrk);
        }
        Genmatch.addUserFloat("idx_kspim_pt", idx_kspim_pt);
        Genmatch.addUserFloat("idx_kspim_iso03", idx_kspim_iso03);
        Genmatch.addUserFloat("idx_kspim_dr", idx_kspim_dr);
	if (idx_pip >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_pip];
            std::tie(idx_pip_thematch, idx_pip_pt, idx_pip_iso03, idx_pip_dr) = computeIsoAndMatch(genCand, vectrk_ttrk);
        }
        Genmatch.addUserFloat("idx_pip_pt", idx_pip_pt);
        Genmatch.addUserFloat("idx_pip_iso03", idx_pip_iso03);
        Genmatch.addUserFloat("idx_pip_dr", idx_pip_dr);
        if (idx_pim >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_pim];
            std::tie(idx_pim_thematch, idx_pim_pt, idx_pim_iso03, idx_pim_dr) = computeIsoAndMatch(genCand, vectrk_ttrk);
        }
        Genmatch.addUserFloat("idx_pim_pt", idx_pim_pt);
        Genmatch.addUserFloat("idx_pim_iso03", idx_pim_iso03);
        Genmatch.addUserFloat("idx_pim_dr", idx_pim_dr);

        Genmatch_out->push_back(Genmatch);
        // Push the matched indices into your output
        genpion_index.push_back(idx_kspip_thematch);
        genpion_index.push_back(idx_kspim_thematch);
        genpion_index.push_back(idx_pip_thematch);
        genpion_index.push_back(idx_pim_thematch);
    }
    return genpion_index;
}


void BDPiFitter::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {

  using std::vector;

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
  if (!theBeamSpotHandle.isValid()) {
    edm::LogError("BDPiFitter") << "No BeamSpot found!";
    return;
  }
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();

  math::XYZPoint theBeamSpotPos(theBeamSpot->position());

  const auto& theMagneticField = iSetup.getData(bFieldToken_);
  edm::Handle<std::vector<reco::Vertex>> vertices;
  iEvent.getByToken(token_vertices, vertices);
  reco::Vertex referenceVtx;
  referenceVtx = vertices->at(0);
  math::XYZPoint referencePos = referenceVtx.position();

  edm::Handle<edm::View<pat::PackedCandidate>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<edm::View<pat::PackedCandidate>> lostTracks;
  iEvent.getByToken(lostTracksToken_, lostTracks);
  if (!tracks.isValid() || !lostTracks.isValid()) {
    edm::LogError("BDPiFitter") << "Track collections not found!";
    return;
  }

  edm::Handle<reco::GenParticleCollection> genParticleHandle_;
  iEvent.getByToken(genParticleToken_, genParticleHandle_);
  if (!genParticleHandle_.isValid()) {
    edm::LogError("BDPiFitter") << "GenParticles not found!";
    return;
  }

  // Filter track
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();
  std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vectrk_ttrk;
  std::vector<std::pair<bool, int>> vectrk_idx;
  std::unique_ptr<pat::CompositeCandidateCollection> tracks_earlyout(new pat::CompositeCandidateCollection());
  for (unsigned int iTrk = 0; iTrk < totalTracks; ++iTrk) {
    const pat::PackedCandidate& trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk - nTracks];
    if (!trk.hasTrackDetails() || !trk.bestTrack()) continue;
    if (trk.charge() == 0) continue;
    const reco::Track* bestTrack = trk.bestTrack();
    if (!bestTrack) continue;
    //if (bestTrack->hitPattern().numberOfValidPixelHits() < 1) continue;
    if (bestTrack->hitPattern().numberOfValidHits() < tkNHitsCut_) continue;
    if (bestTrack->pt() < tkPtCut_) continue;
    if (fabs(trk.pdgId()) != 211) continue; // pions only
    //if (iTrk < nTracks && !trk.trackHighPurity()) continue;
    if (bestTrack->normalizedChi2() > tkChi2Cut_) continue;
    const reco::Track& tmpTrack = trk.pseudoTrack();
    double ipXY_bs = std::abs(tmpTrack.dxy(theBeamSpotPos));
    double ipZ_bs = std::abs(tmpTrack.dz(theBeamSpotPos));
    double ipXY_pv  = std::abs(tmpTrack.dxy(referencePos));
    double ipZ_pv  = std::abs(tmpTrack.dz(referencePos));
    double ipsigXY_bs = std::abs(tmpTrack.dxy(theBeamSpotPos) / tmpTrack.dxyError());
    double ipsigZ_bs  = std::abs(tmpTrack.dz(theBeamSpotPos) / tmpTrack.dzError());
    double ipsigXY_pv = std::abs(tmpTrack.dxy(referencePos) / tmpTrack.dxyError());
    double ipsigZ_pv  = std::abs(tmpTrack.dz(referencePos) / tmpTrack.dzError());
    const reco::TransientTrack tmpTransient( (*trk.bestTrack()) , &theMagneticField);  
    vectrk_ttrk.emplace_back(trk, std::move(tmpTransient));
    if (iTrk < nTracks) {
        vectrk_idx.push_back(std::make_pair(true, iTrk));
    } else {
        vectrk_idx.push_back(std::make_pair(false, iTrk - nTracks));
    }
    pat::CompositeCandidate pcand;
    pcand.setP4(trk.p4());
    pcand.setCharge(trk.charge());
    pcand.setVertex(trk.vertex());
    pcand.setPdgId(trk.pdgId());
    pcand.addUserFloat("dxy",  trk.dxy());
    pcand.addUserFloat("dxyS", trk.dxy() / trk.dxyError()); // cut 2 in cmssw
    pcand.addUserFloat("dz",   trk.dz());
    pcand.addUserFloat("dzS",  trk.dz() / trk.dzError()); // not cut
    pcand.addUserFloat("bt_pt",         bestTrack->pt());
    if (iTrk < nTracks && trk.trackHighPurity()){
        pcand.addUserFloat("trackHighPurity", 1);
    }else if (iTrk < nTracks){
        pcand.addUserFloat("trackHighPurity", 0);
    }else{
	pcand.addUserFloat("trackHighPurity", -1);
    }
    pcand.addUserFloat("ipXY_bs",       ipXY_bs);
    pcand.addUserFloat("ipZ_bs",        ipZ_bs);
    pcand.addUserFloat("ipXY_pv",       ipXY_pv);
    pcand.addUserFloat("ipZ_pv",        ipZ_pv);
    pcand.addUserFloat("ipsigXY_bs",       ipsigXY_bs);
    pcand.addUserFloat("ipsigZ_bs",        ipsigZ_bs);
    pcand.addUserFloat("ipsigXY_pv",       ipsigXY_pv);
    pcand.addUserFloat("ipsigZ_pv",        ipsigZ_pv);
    pcand.addUserFloat("ptErr",         bestTrack->ptError());
    pcand.addUserFloat("normChi2",      bestTrack->normalizedChi2());
    pcand.addUserInt("nValidPixelHits", bestTrack->hitPattern().numberOfValidPixelHits());
    pcand.addUserInt("nValidHits",      bestTrack->hitPattern().numberOfValidHits());
    tracks_earlyout -> emplace_back(pcand);
  }

  std::unique_ptr<pat::CompositeCandidateCollection> Genmatch_out(new pat::CompositeCandidateCollection());
  auto LooseKshorts_out = std::make_unique<reco::VertexCompositePtrCandidateCollection>();
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  std::unique_ptr<TransientTrackCollection> trans_out( new TransientTrackCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> D0_out(new pat::CompositeCandidateCollection());

  // Gen Study
  auto genpion_index = FindBDecays(genParticleHandle_, Genmatch_out.get(), vectrk_ttrk);

  // If early stop
  if(onlygen){
      if(isSignal)  iEvent.put(std::move(Genmatch_out), "Genmatch");
      iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");
      iEvent.put(std::move(LooseKshorts_out), "LooseKshort");
      iEvent.put(std::move(ret_val), "SelectedV0Collection");
      iEvent.put(std::move(trans_out), "SelectedV0TransientCollection");
      iEvent.put(std::move(D0_out), "D0");
      return;
  }

  // Get extra 2 pions
  //std::cout<<"main loop to find two pion from D0 "<< tracks->size()<<" "<< pu_tracks->size()<<" "<< ttracks->size() <<std::endl;
  std::vector<std::pair<unsigned int, unsigned int>> vectrk_pions;
  // Loop over track pairs
  for (unsigned int trkidx3 = 0; trkidx3 < vectrk_ttrk.size(); ++trkidx3) {
    const auto& Trk3 = vectrk_ttrk[trkidx3].first;
    const auto& tTrk3 = vectrk_ttrk[trkidx3].second;
    if(onlymatched){
        if(int(trkidx3)!= genpion_index.at(2) && int(trkidx3)!= genpion_index.at(3)){
          //std::cout<<"p3 "<<trkidx3<<" "<<genpion_index.at(2)<<" "<<genpion_index.at(3)<<std::endl;
          continue;
        }
    }
    for (unsigned int trkidx4 = trkidx3 + 1; trkidx4 < vectrk_ttrk.size(); ++trkidx4) {
      const auto& Trk4 = vectrk_ttrk[trkidx4].first;
      const auto& tTrk4 = vectrk_ttrk[trkidx4].second;
      if(onlymatched){
          if(int(trkidx4)!= genpion_index.at(2) && int(trkidx4)!= genpion_index.at(3)){
            //std::cout<<trkidx4<<" "<<genpion_index.at(2)<<" "<<genpion_index.at(3)<<std::endl;
            continue;
          }
      }
      // Opposite charge requirement
      if (Trk3.charge() * Trk4.charge() >= 0) continue;
  
      const auto& tPosTrk = (Trk3.charge() > 0) ? tTrk3 : tTrk4;
      const auto& tNegTrk = (Trk3.charge() < 0) ? tTrk3 : tTrk4;
  
      //double dz_diff = std::abs(tPosTrk.track().dz(*theBeamSpot) - tNegTrk.track().dz(*theBeamSpot));
      //if (dz_diff > dzCut_) continue;
      //double dxy_sign = tPosTrk.track().dxy(*theBeamSpot) * tNegTrk.track().dxy(*theBeamSpot);
      //if (dxy_sign > 0) continue;
  
      // Run the vertex fitter
      KinVtxFitter D0pipi_fitter(
        {tPosTrk, tNegTrk},      // transient tracks
        {piMass, piMass},        // masses
        {K_SIGMA, K_SIGMA}       // sigma (uncertainty)
      );
  
      if (!D0pipi_fitter.success()) continue;
      if (D0pipi_fitter.chi2() < 0 || D0pipi_fitter.dof() < 0) continue;
      if (D0pipi_fitter.prob() < 0.01) continue;
  
      // Store the indices of the tracks for later combinations
      vectrk_pions.push_back(std::make_pair(trkidx3, trkidx4));
    }
  }


  //std::cout<<"vectrk_ttrk.size() "<< vectrk_ttrk.size()<<" vectrk_pions.size() "<< vectrk_pions.size()<<std::endl;
  // Reco V0
  int v0_idx=-1;
  for (unsigned int trkidx1 = 0; trkidx1 < vectrk_ttrk.size(); ++trkidx1) {
    for (unsigned int trkidx2 = trkidx1 + 1; trkidx2 < vectrk_ttrk.size(); ++trkidx2) {
	if(onlymatched){
            if( (int(trkidx1)!= genpion_index.at(0) && int(trkidx1)!= genpion_index.at(1)) || (int(trkidx2)!= genpion_index.at(0) && int(trkidx2)!= genpion_index.at(1))){
              //std::cout<<"V0 "<<trkidx1<<" "<<genpion_index.at(0)<<" "<<genpion_index.at(1)<<std::endl;
              continue;
            }
        }
        const auto Trk1 = (vectrk_ttrk[trkidx1].first.charge() > 0) ? vectrk_ttrk[trkidx1].first : vectrk_ttrk[trkidx2].first; 
        const auto Trk2 = (vectrk_ttrk[trkidx1].first.charge() < 0) ? vectrk_ttrk[trkidx1].first : vectrk_ttrk[trkidx2].first;
        const auto& tTrk1 = (vectrk_ttrk[trkidx1].first.charge() > 0) ? vectrk_ttrk[trkidx1].second : vectrk_ttrk[trkidx2].second;
        const auto& tTrk2 = (vectrk_ttrk[trkidx1].first.charge() < 0) ? vectrk_ttrk[trkidx1].second : vectrk_ttrk[trkidx2].second;
        if (Trk1.charge() * Trk2.charge() >= 0) continue;
        ////std::cout<<"OS tracks:"<<trkidx1<<" "<<trkidx2<<std::endl;

        //// step1: DCA requirement before anything else
        const auto& posImpact = tTrk1.impactPointTSCP();
        const auto& negImpact = tTrk2.impactPointTSCP();
        if (!posImpact.isValid() || !negImpact.isValid()) continue;
        FreeTrajectoryState const& posState = posImpact.theState();
        FreeTrajectoryState const& negState = negImpact.theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(posState, negState);
        if (!cApp.status()) continue;
	//  the distance between the two trajectories at their closest approach in R-phi
        float dca = std::abs(cApp.distance());
        // the POCA should at least be in the sensitive volume
        GlobalPoint cxPt = cApp.crossingPoint();
        const double cxPtR2 = cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y();
        //if (cxPtR2 > 120. * 120. || std::abs(cxPt.z()) > 300.) continue;
        //// allow for different DCA cuts depending on position of POCA
        //if (cxPtR2 < innerOuterTkDCAThreshold_ * innerOuterTkDCAThreshold_) {
        //  if (dca > innerTkDCACut_) continue;
        //} else {
        //  if (dca > outerTkDCACut_) continue;
        //}
        // require the tracks point in the same quadrant
	TrajectoryStateClosestToPoint posTSCP = tTrk1.trajectoryStateClosestToPoint(cxPt);
        TrajectoryStateClosestToPoint negTSCP = tTrk2.trajectoryStateClosestToPoint(cxPt);
        if (!posTSCP.isValid() || !negTSCP.isValid())          continue;
        float trk1_dot_trk2 = posTSCP.momentum().dot(negTSCP.momentum());

	// step2 seperate dca requirement
	//distance closest approach in x,y wrt beam spot
        std::pair<double, double> DCA_trk1_beamspot = computeDCA(tTrk1, *theBeamSpot);
        std::pair<double, double> DCA_trk2_beamspot = computeDCA(tTrk2, *theBeamSpot);
        //if (DCASig >  dcaSig_  && dcaSig_ > 0) continue;

  
        // step3: Pre-fit mPiPi mass cut
        double p1Mag2 = posTSCP.momentum().mag2();
        double p2Mag2 = negTSCP.momentum().mag2();
        double totalE = std::sqrt(p1Mag2 + piMassSquared) + std::sqrt(p2Mag2 + piMassSquared);
        double totalESq = totalE * totalE;
        double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
        double massSquared = totalESq - totalPSq;
        //if (massSquared > mPiPiCut_ * mPiPiCut_)         continue;
        //std::cout << "mPiPi^2 = " << massSquared << " passes pre-fit cut." << std::endl;

	// step4: Fit! 
        // Create dummy covariance (error matrix) for the initial vertex position
        GlobalError dummyError(1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-3);
        std::vector<reco::TransientTrack> transTracks{tTrk1, tTrk2};
        TransientVertex theRecoVertex(cxPt, dummyError, transTracks, 1.0e-3);
        KalmanVertexFitter theKalmanFitter(true);
        theRecoVertex = theKalmanFitter.vertex(transTracks);
        if (!theRecoVertex.isValid()) continue;
        reco::Vertex theVtx_KLM = reco::Vertex(theRecoVertex);
        // if (theVtx_KLM.normalizedChi2() > vtxChi2Cut_ || theVtx_KLM.normalizedChi2() < 0)  continue;
        // Get fitted vertex position
        GlobalPoint vtxPos(theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z());
        //std::cout << "A valid vertex found at ("
        //          << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z() << ")"
        //          << std::endl;
	        
	// step5: check the fit vertex
        // 2D decay significance
        SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx_KLM.covariance();
        SVector3 distVecXY(vtxPos.x() - theBeamSpotPos.x(), vtxPos.y() - theBeamSpotPos.y(), 0.);
        double distMagXY = ROOT::Math::Mag(distVecXY);
        //if (distMagXY < vtxDecayXYCut_) continue;
        double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
        //if (distMagXY < vtxDecaySigXYCut_ * sigmaDistMagXY) continue;
        // 3D decay significance
        SVector3 distVecXYZ(vtxPos.x() - theBeamSpotPos.x(), vtxPos.y() - theBeamSpotPos.y(), vtxPos.z() - theBeamSpotPos.z());
        double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
        double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
        //if (distMagXYZ / sigmaDistMagXYZ < vtxDecaySigXYZCut_) continue;
        //std::cout<<"Good 2D / 3D sig "<<std::endl;

        std::unique_ptr<TrajectoryStateClosestToPoint> trajPlus;
        std::unique_ptr<TrajectoryStateClosestToPoint> trajMins;
        std::vector<reco::TransientTrack> theRefTracks;
        if (theRecoVertex.hasRefittedTracks()) {
          theRefTracks = theRecoVertex.refittedTracks();
        }

        if (theRefTracks.size() > 1) {
          reco::TransientTrack* thePositiveRefTrack = nullptr;
          reco::TransientTrack* theNegativeRefTrack = nullptr;
          for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end();
               ++iTrack) {
            if (iTrack->track().charge() > 0.) {
              thePositiveRefTrack = &*iTrack;
            } else if (iTrack->track().charge() < 0.) {
              theNegativeRefTrack = &*iTrack;
            }
          }
          if (thePositiveRefTrack == nullptr || theNegativeRefTrack == nullptr)
            continue;
          trajPlus =
              std::make_unique<TrajectoryStateClosestToPoint>(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos));
          trajMins =
              std::make_unique<TrajectoryStateClosestToPoint>(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos));
        } else {
          trajPlus =
              std::make_unique<TrajectoryStateClosestToPoint>(tTrk1.trajectoryStateClosestToPoint(vtxPos));
          trajMins =
              std::make_unique<TrajectoryStateClosestToPoint>(tTrk2.trajectoryStateClosestToPoint(vtxPos));
        }

        if (trajPlus.get() == nullptr || trajMins.get() == nullptr || !trajPlus->isValid() || !trajMins->isValid())
          continue;

        GlobalVector positiveP(trajPlus->momentum());
        GlobalVector negativeP(trajMins->momentum());
        GlobalVector totalP(positiveP + negativeP);

        // 2D pointing angle
        double dx = theVtx_KLM.x() - theBeamSpotPos.x();
        double dy = theVtx_KLM.y() - theBeamSpotPos.y();
        double px = totalP.x();
        double py = totalP.y();
        double angleXY = (dx * px + dy * py) / (sqrt(dx * dx + dy * dy) * sqrt(px * px + py * py));
        //if (angleXY < cosThetaXYCut_) continue;

        // 3D pointing angle
        double dz = theVtx_KLM.z() - theBeamSpotPos.z();
        double pz = totalP.z();
        double angleXYZ =
            (dx * px + dy * py + dz * pz) / (sqrt(dx * dx + dy * dy + dz * dz) * sqrt(px * px + py * py + pz * pz));
        //if (angleXYZ < cosThetaXYZCut_) continue;

        // step6: make kShort
        // calculate total energy of V0
        double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
        double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
        double kShortETot = piPlusE + piMinusE;
        const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);

        reco::Particle::Point vtx(theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z());
        const reco::Vertex::CovarianceMatrix vtxCov(theVtx_KLM.covariance());
        double vtxChi2(theVtx_KLM.chi2());
        double vtxNdof(theVtx_KLM.ndof());

        // Create the VertexCompositeCandidate object that will be stored in the Event
	reco::VertexCompositePtrCandidate theKshort_cand(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
	bool fromtrack1 = vectrk_idx[trkidx1].first;
        bool fromtrack2 = vectrk_idx[trkidx2].first;
        
        int posIdx, negIdx;
        bool fromPos, fromNeg;
        
        // Identify which one is positive and negative
        if (vectrk_ttrk[trkidx1].first.charge() > 0) {
            posIdx = vectrk_idx[trkidx1].second;
            negIdx = vectrk_idx[trkidx2].second;
            fromPos = fromtrack1;
            fromNeg = fromtrack2;
        } else {
            posIdx = vectrk_idx[trkidx2].second;
            negIdx = vectrk_idx[trkidx1].second;
            fromPos = fromtrack2;
            fromNeg = fromtrack1;
        }
        
        // Create edm::Ptr's
        edm::Ptr<pat::PackedCandidate> positivePtr = fromPos ?
            edm::Ptr<pat::PackedCandidate>(tracks, posIdx) :
            edm::Ptr<pat::PackedCandidate>(lostTracks, posIdx);
        
        edm::Ptr<pat::PackedCandidate> negativePtr = fromNeg ?
            edm::Ptr<pat::PackedCandidate>(tracks, negIdx) :
            edm::Ptr<pat::PackedCandidate>(lostTracks, negIdx);
        
        // Add daughters
        theKshort_cand.addDaughter(positivePtr);
        theKshort_cand.addDaughter(negativePtr);
        theKshort_cand.setPdgId(310);
        //if ( abs(theKshort_cand.mass() - kShortMass) > kShortMassCut_) continue;

	pat::CompositeCandidate V0cand;
        V0cand.setP4(theKshort_cand.p4());
	V0cand.setVertex( reco::Candidate::Point(theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z()));
        V0cand.addUserInt("Trk1_idx", posIdx);
	V0cand.addUserFloat("Trk1_pt",  Trk1.pt());
        V0cand.addUserFloat("Trk1_eta", Trk1.eta());
        V0cand.addUserFloat("Trk1_phi", Trk1.phi());
        V0cand.addUserFloat("Trk1_dca_beamspot", DCA_trk1_beamspot.first);
        V0cand.addUserFloat("Trk1_dcaErr_beamspot", DCA_trk1_beamspot.second);
        V0cand.addUserFloat("Trk1_PixelHits", Trk1.bestTrack()->hitPattern().numberOfValidPixelHits());
        V0cand.addUserFloat("Trk1_bestTrackHits", Trk1.bestTrack()->hitPattern().numberOfValidHits());
        V0cand.addUserFloat("Trk1_pseudoTrackHits", Trk1.pseudoTrack().numberOfValidHits());
	V0cand.addUserFloat("Trk1_ipsigXY", std::abs(Trk1.pseudoTrack().dxy(theBeamSpotPos) / Trk1.pseudoTrack().dxyError()));
        V0cand.addUserFloat("Trk1_ipsigZ", std::abs(Trk1.pseudoTrack().dz(theBeamSpotPos) / Trk1.pseudoTrack().dzError()));
        V0cand.addUserFloat("Trk1_normalizedChi2", Trk1.pseudoTrack().normalizedChi2());

        V0cand.addUserInt("Trk2_idx", negIdx);
	V0cand.addUserFloat("Trk2_pt",  Trk2.pt());
        V0cand.addUserFloat("Trk2_eta", Trk2.eta());
        V0cand.addUserFloat("Trk2_phi", Trk2.phi());
        V0cand.addUserFloat("Trk2_dca_beamspot", DCA_trk2_beamspot.first);
        V0cand.addUserFloat("Trk2_dcaErr_beamspot", DCA_trk2_beamspot.second);
        V0cand.addUserFloat("Trk2_PixelHits", Trk2.bestTrack()->hitPattern().numberOfValidPixelHits());
        V0cand.addUserFloat("Trk2_bestTrackHits", Trk2.bestTrack()->hitPattern().numberOfValidHits());
        V0cand.addUserFloat("Trk2_pseudoTrackHits", Trk2.pseudoTrack().numberOfValidHits());
        V0cand.addUserFloat("Trk2_ipsigXY", std::abs(Trk2.pseudoTrack().dxy(theBeamSpotPos) / Trk2.pseudoTrack().dxyError()));
        V0cand.addUserFloat("Trk2_ipsigZ", std::abs(Trk2.pseudoTrack().dz(theBeamSpotPos) / Trk2.pseudoTrack().dzError()));
        V0cand.addUserFloat("Trk2_normalizedChi2", Trk2.pseudoTrack().normalizedChi2());
	// closest xing point
        V0cand.addUserFloat("trk1_dot_trk2", trk1_dot_trk2);
        V0cand.addUserFloat("sigmaDistMagXY", sigmaDistMagXY);
        V0cand.addUserFloat("sigmaDistMagXYZ", sigmaDistMagXYZ);
        V0cand.addUserFloat("cxPtx", cxPt.x());
        V0cand.addUserFloat("cxPty", cxPt.y());
        V0cand.addUserFloat("cxPtz", cxPt.z());
        V0cand.addUserFloat("cxPtR2", cxPtR2);
        V0cand.addUserFloat("dca", dca);
	V0cand.addUserFloat("massSquared", massSquared);

	// KLM fit
	V0cand.addUserFloat("KLM_vtx_x", theVtx_KLM.x());
	V0cand.addUserFloat("KLM_vtx_y", theVtx_KLM.y());
	V0cand.addUserFloat("KLM_vtx_z", theVtx_KLM.z());
        V0cand.addUserFloat("KLM_chi2", theVtx_KLM.chi2());
        V0cand.addUserFloat("KLM_ndof", theVtx_KLM.ndof());
	V0cand.addUserFloat("KLM_normalizedChi2", theVtx_KLM.normalizedChi2());
	//V0cand.addUserFloat("KLM_prob", theVtx_KLM.prob());
        V0cand.addUserFloat("KLM_cos_theta_XY", angleXY);
        V0cand.addUserFloat("KLM_cos_theta_XYZ", angleXYZ);
	V0cand.addUserFloat("KLM_Trk1_pt",  tTrk1.track().pt());
	V0cand.addUserFloat("KLM_Trk1_eta", tTrk1.track().eta());
	V0cand.addUserFloat("KLM_Trk1_phi", tTrk1.track().phi());
        V0cand.addUserFloat("KLM_Trk2_pt",  tTrk2.track().pt());
        V0cand.addUserFloat("KLM_Trk2_eta", tTrk2.track().eta());
        V0cand.addUserFloat("KLM_Trk2_phi", tTrk2.track().phi());
        V0cand.addUserFloat("KLM_Ks0_pt",  theKshort_cand.pt());
        V0cand.addUserFloat("KLM_Ks0_eta", theKshort_cand.eta());
        V0cand.addUserFloat("KLM_Ks0_phi", theKshort_cand.phi());
        V0cand.addUserFloat("KLM_Ks0_mass",   theKshort_cand.mass());
        V0cand.addUserFloat("KLM_distMagXY",   distMagXY);
        V0cand.addUserFloat("KLM_sigmaDistMagXY",   sigmaDistMagXY);
	V0cand.addUserFloat("KLM_distMagXYZ",   distMagXYZ);
	V0cand.addUserFloat("KLM_sigmaDistMagXYZ",   sigmaDistMagXYZ);


	// Fit again with KinVtxFitter
        KinVtxFitter theKinVtxFitter(
        {tTrk1, tTrk2},
        {piMass, piMass},
        {K_SIGMA, K_SIGMA} );
        if (!theKinVtxFitter.success()) continue;
        if ( theKinVtxFitter.chi2()<0 || theKinVtxFitter.dof()<0) continue;
        V0cand.addUserFloat("Kin_vtx_x", theKinVtxFitter.fitted_vtx().x());
        V0cand.addUserFloat("Kin_vtx_y", theKinVtxFitter.fitted_vtx().y());
        V0cand.addUserFloat("Kin_vtx_z", theKinVtxFitter.fitted_vtx().z());
	V0cand.addUserFloat("Kin_chi2", theKinVtxFitter.chi2());
	//V0cand.addUserFloat("Kin_normalizedChi2", theKinVtxFitter.normalizedChi2());
	V0cand.addUserFloat("Kin_dof", theKinVtxFitter.dof());
        V0cand.addUserFloat("Kin_prob", theKinVtxFitter.prob());
        V0cand.addUserFloat("Kin_pt",  theKinVtxFitter.fitted_p4().pt());
        V0cand.addUserFloat("Kin_eta", theKinVtxFitter.fitted_p4().eta());
        V0cand.addUserFloat("Kin_phi", theKinVtxFitter.fitted_p4().phi());
	V0cand.addUserFloat("Kin_mass", theKinVtxFitter.fitted_candidate().mass());
        V0cand.addUserFloat("Kin_massErr", sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
        V0cand.addUserFloat("Kin_cos_theta_2D", cos_theta_2D(theKinVtxFitter, *theBeamSpot, V0cand.p4()));
        V0cand.addUserFloat("Kin_fitted_cos_theta_2D", cos_theta_2D(theKinVtxFitter, *theBeamSpot, theKinVtxFitter.fitted_p4()));
        auto lxy = l_xy(theKinVtxFitter, *theBeamSpot);
        V0cand.addUserFloat("Kin_l_xy", lxy.value());
        V0cand.addUserFloat("Kin_l_xy_unc", lxy.error());
        const auto& covMatrix = theKinVtxFitter.fitted_vtx_uncertainty();
        V0cand.addUserFloat("Kin_vtx_cxx", covMatrix.cxx());
        V0cand.addUserFloat("Kin_vtx_cyy", covMatrix.cyy());
        V0cand.addUserFloat("Kin_vtx_czz", covMatrix.czz());
        V0cand.addUserFloat("Kin_vtx_cyx", covMatrix.cyx());
        V0cand.addUserFloat("Kin_vtx_czx", covMatrix.czx());
        V0cand.addUserFloat("Kin_vtx_czy", covMatrix.czy());
        V0cand.addUserFloat("Kin_trk1_pt", theKinVtxFitter.daughter_p4(0).pt());
        V0cand.addUserFloat("Kin_trk1_eta", theKinVtxFitter.daughter_p4(0).eta());
        V0cand.addUserFloat("Kin_trk1_phi", theKinVtxFitter.daughter_p4(0).phi());
        V0cand.addUserFloat("Kin_trk2_pt", theKinVtxFitter.daughter_p4(1).pt());
        V0cand.addUserFloat("Kin_trk2_eta", theKinVtxFitter.daughter_p4(1).eta());
        V0cand.addUserFloat("Kin_trk2_phi", theKinVtxFitter.daughter_p4(1).phi());
        // save
        ret_val->push_back(V0cand);  // for save
        LooseKshorts_out->emplace_back(theKshort_cand); // for next sequence
        auto V0TT = theKinVtxFitter.fitted_candidate_ttrk();
        trans_out->emplace_back(V0TT);
	v0_idx+=1;

	// Get extra 2 pions
        //std::cout<<"main loop to find two pion from D0 "<< tracks->size()<<" "<< pu_tracks->size()<<" "<< ttracks->size() <<std::endl;
        for (unsigned int pionsidx = 0; pionsidx < vectrk_pions.size(); ++pionsidx) {
	      unsigned int trkidx3 = vectrk_pions[pionsidx].first;
              unsigned int trkidx4 = vectrk_pions[pionsidx].second;
              if(trkidx3 == trkidx1 || trkidx3 == trkidx2) continue;
              if(trkidx4 == trkidx1 || trkidx4 == trkidx2) continue;
              const auto Trk3 = vectrk_ttrk[trkidx3].first;
              const auto Trk4 = vectrk_ttrk[trkidx4].first;
              const auto& tTrk3 = vectrk_ttrk[trkidx3].second;
              const auto& tTrk4 = vectrk_ttrk[trkidx4].second;

              const auto& posImpact_pion = tTrk3.impactPointTSCP();
              const auto& negImpact_pion = tTrk4.impactPointTSCP();
              if (!posImpact_pion.isValid() || !negImpact_pion.isValid()) continue;
              FreeTrajectoryState const& posState_pion = posImpact_pion.theState();
              FreeTrajectoryState const& negState_pion = negImpact_pion.theState();
              ClosestApproachInRPhi cApp_pion;
              cApp_pion.calculate(posState_pion, negState_pion);
              if (!cApp_pion.status()) continue;
              //  the distance between the two trajectories at their closest approach in R-phi
              float dca_pion = std::abs(cApp_pion.distance());
              // the POCA should at least be in the sensitive volume
              GlobalPoint cxPt_pion = cApp_pion.crossingPoint();
              const double cxPtR2_pion = cxPt_pion.x() * cxPt_pion.x() + cxPt_pion.y() * cxPt_pion.y();
              // require the tracks point in the same quadrant
              TrajectoryStateClosestToPoint posTSCP_pion = tTrk3.trajectoryStateClosestToPoint(cxPt);
              TrajectoryStateClosestToPoint negTSCP_pion = tTrk4.trajectoryStateClosestToPoint(cxPt);
              if (!posTSCP_pion.isValid() || !negTSCP_pion.isValid())          continue;
              float trk3_dot_trk4 = posTSCP_pion.momentum().dot(negTSCP_pion.momentum());
              std::pair<double, double> DCA_trk3_beamspot = computeDCA(tTrk3, *theBeamSpot);
              std::pair<double, double> DCA_trk4_beamspot = computeDCA(tTrk4, *theBeamSpot);

  	      //D0
              KinVtxFitter D0_KinFitter(
                  {tTrk3, tTrk4, V0TT },
                  {piMass, piMass,  V0cand.userFloat("Kin_mass") },
                  {K_SIGMA, K_SIGMA, V0cand.userFloat("Kin_massErr")});
              if ( !D0_KinFitter.success() ) continue;
              if ( D0_KinFitter.chi2()<0 || D0_KinFitter.dof()<0) continue;
              //if ( D0_KinFitter.prob()<0.01) continue;
              auto D0_Kinfit_p4 = D0_KinFitter.fitted_p4();
              //if(abs(D0_Kinfit_p4.mass()-D0Mass)>0.2) continue;

              pat::CompositeCandidate D0;
              D0.addUserInt("Ks0_idx", v0_idx);
              D0.setP4(D0_KinFitter.fitted_p4());
              D0.setVertex( reco::Candidate::Point(
                  D0_KinFitter.fitted_vtx().x(),
                  D0_KinFitter.fitted_vtx().y(),
                  D0_KinFitter.fitted_vtx().z())
              );
	      // raw
              D0.addUserInt("Trk3_idx", trkidx3 );
              D0.addUserFloat("Trk3_pt",  Trk3.pt());
              D0.addUserFloat("Trk3_eta", Trk3.eta());
              D0.addUserFloat("Trk3_phi", Trk3.phi());
              D0.addUserFloat("Trk3_dca_beamspot", DCA_trk3_beamspot.first);
              D0.addUserFloat("Trk3_dcaErr_beamspot", DCA_trk3_beamspot.second);
              D0.addUserFloat("Trk3_PixelHits", Trk3.bestTrack()->hitPattern().numberOfValidPixelHits());
              D0.addUserFloat("Trk3_bestTrackHits", Trk3.bestTrack()->hitPattern().numberOfValidHits());
              D0.addUserFloat("Trk3_pseudoTrackHits", Trk3.pseudoTrack().numberOfValidHits());
              D0.addUserFloat("Trk3_ipsigXY", std::abs(Trk3.pseudoTrack().dxy(theBeamSpotPos) / Trk3.pseudoTrack().dxyError()));
              D0.addUserFloat("Trk3_ipsigZ", std::abs(Trk3.pseudoTrack().dz(theBeamSpotPos) / Trk3.pseudoTrack().dzError()));
              D0.addUserFloat("Trk3_normalizedChi2", Trk3.pseudoTrack().normalizedChi2());

              D0.addUserInt("Trk4_idx", trkidx4 );
	      D0.addUserFloat("Trk4_pt",  Trk4.pt());
              D0.addUserFloat("Trk4_eta", Trk4.eta());
              D0.addUserFloat("Trk4_phi", Trk4.phi());
              D0.addUserFloat("Trk4_dca_beamspot", DCA_trk4_beamspot.first);
              D0.addUserFloat("Trk4_dcaErr_beamspot", DCA_trk4_beamspot.second);
              D0.addUserFloat("Trk4_PixelHits", Trk4.bestTrack()->hitPattern().numberOfValidPixelHits());
              D0.addUserFloat("Trk4_bestTrackHits", Trk4.bestTrack()->hitPattern().numberOfValidHits());
              D0.addUserFloat("Trk4_pseudoTrackHits", Trk4.pseudoTrack().numberOfValidHits());
              D0.addUserFloat("Trk4_ipsigXY", std::abs(Trk4.pseudoTrack().dxy(theBeamSpotPos) / Trk4.pseudoTrack().dxyError()));
              D0.addUserFloat("Trk4_ipsigZ", std::abs(Trk4.pseudoTrack().dz(theBeamSpotPos) / Trk4.pseudoTrack().dzError()));
              D0.addUserFloat("Trk4_normalizedChi2", Trk4.pseudoTrack().normalizedChi2());
              // closest xing point
              D0.addUserFloat("trk3_dot_trk4", trk3_dot_trk4);
              D0.addUserFloat("cxPtx", cxPt_pion.x());
              D0.addUserFloat("cxPty", cxPt_pion.y());
              D0.addUserFloat("cxPtz", cxPt_pion.z());
              D0.addUserFloat("cxPtR2", cxPtR2_pion);
              D0.addUserFloat("dca", dca_pion);

	      // fitted
              D0.addUserFloat("Kin_vtx_x", D0_KinFitter.fitted_vtx().x());
              D0.addUserFloat("Kin_vtx_y", D0_KinFitter.fitted_vtx().y());
              D0.addUserFloat("Kin_vtx_z", D0_KinFitter.fitted_vtx().z());
              D0.addUserFloat("Kin_chi2", D0_KinFitter.chi2());
              D0.addUserFloat("Kin_dof", D0_KinFitter.dof());
              //D0.addUserFloat("Kin_normalizedChi2", D0_KinFitter.normalizedChi2());
	      D0.addUserFloat("Kin_prob", D0_KinFitter.prob());
              D0.addUserFloat("Kin_cos_theta_2D", cos_theta_2D(D0_KinFitter, *theBeamSpot, D0_Kinfit_p4));
              auto D0_lxy = l_xy(D0_KinFitter, *theBeamSpot);
              D0.addUserFloat("Kin_l_xy", D0_lxy.value());
              D0.addUserFloat("Kin_l_xy_unc", D0_lxy.error());
	      D0.addUserFloat("Kin_trk3_pt", D0_KinFitter.daughter_p4(0).pt());
              D0.addUserFloat("Kin_trk3_eta", D0_KinFitter.daughter_p4(0).eta());
              D0.addUserFloat("Kin_trk3_phi", D0_KinFitter.daughter_p4(0).phi());
              D0.addUserFloat("Kin_trk4_pt", D0_KinFitter.daughter_p4(1).pt());
              D0.addUserFloat("Kin_trk4_eta", D0_KinFitter.daughter_p4(1).eta());
              D0.addUserFloat("Kin_trk4_phi", D0_KinFitter.daughter_p4(1).phi());
              D0.addUserFloat("Kin_Ks0_pt", D0_KinFitter.daughter_p4(2).pt());
              D0.addUserFloat("Kin_Ks0_eta", D0_KinFitter.daughter_p4(2).eta());
              D0.addUserFloat("Kin_Ks0_phi", D0_KinFitter.daughter_p4(2).phi());
              D0.addUserFloat("Kin_D0_pt", D0_Kinfit_p4.pt() );
              D0.addUserFloat("Kin_D0_eta", D0_Kinfit_p4.eta() );
              D0.addUserFloat("Kin_D0_phi", D0_Kinfit_p4.phi() );
              D0.addUserFloat("Kin_D0_mass", D0_Kinfit_p4.mass() );

              reco::Candidate::Point Ks_decay_vtx(V0cand.userFloat("Kin_vtx_x"), V0cand.userFloat("Kin_vtx_y"), V0cand.userFloat("Kin_vtx_z") );
              reco::Candidate::Point Ks_production_vtx(D0.userFloat("Kin_vtx_x"), D0.userFloat("Kin_vtx_y"), D0.userFloat("Kin_vtx_z"));
	      math::XYZVector Ks_flight = Ks_decay_vtx - Ks_production_vtx;
              float Ks_flight_distance = Ks_flight.R();
              float Ks_flight_distance_2D = sqrt(Ks_flight.x() * Ks_flight.x() + Ks_flight.y() * Ks_flight.y());
              D0.addUserFloat("Ks_flight_distance", Ks_flight_distance);
              D0.addUserFloat("Ks_flight_distance_2D", Ks_flight_distance_2D);
              D0_out->push_back(D0);
        }
    }
  }

  if(isSignal)  iEvent.put(std::move(Genmatch_out), "Genmatch");
  iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");    
  iEvent.put(std::move(LooseKshorts_out), "LooseKshort");
  iEvent.put(std::move(ret_val), "SelectedV0Collection");
  iEvent.put(std::move(trans_out), "SelectedV0TransientCollection");
  iEvent.put(std::move(D0_out), "D0");
  return;

}

DEFINE_FWK_MODULE(BDPiFitter);

