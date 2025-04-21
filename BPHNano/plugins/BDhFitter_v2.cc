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
  const float  pion_sigma = piMass*1.e-6;
  const double BuMass = 5.27934;

}  // namespace

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

class BDhFitter_v2 : public edm::global::EDProducer<> {

public:
  explicit BDhFitter_v2(const edm::ParameterSet &theParameters):
      bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      token_beamSpot(consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"))),
      token_vertices(consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("tracks"))),
      lostTracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("lostTracks"))),
      genParticleToken_(consumes<reco::GenParticleCollection>(theParameters.getParameter<edm::InputTag>("genParticle"))),
      tkNHitsCut_( theParameters.getParameter<int>("tkNHitsCut")),
      tkPtCut_( theParameters.getParameter<double>("tkPtCut")),
      tkEtaCut_( theParameters.getParameter<double>("tkEtaCut")),
      tkChi2Cut_( theParameters.getParameter<double>("tkChi2Cut")),
      vtxChi2Cut_( theParameters.getParameter<double>("vtxChi2Cut")),
      vtxDecaySigXYCut_( theParameters.getParameter<double>("vtxDecaySigXYCut")),
      vtxDecaySigXYZCut_( theParameters.getParameter<double>("vtxDecaySigXYZCut")),
      cosThetaXYCut_( theParameters.getParameter<double>("cosThetaXYCut")),
      cosThetaXYZCut_( theParameters.getParameter<double>("cosThetaXYZCut")),
      mPiPiCut_( theParameters.getParameter<double>("mPiPiCut")),
      kShortMassCut_( theParameters.getParameter<double>("kShortMassCut")),
      D0MassCut_( theParameters.getParameter<double>("D0MassCut")),
      BMassCut_( theParameters.getParameter<double>("BMassCut")),
      savetrack_( theParameters.getParameter<bool>("savetrack")),
      verbose( theParameters.getParameter<int>("verbose"))
	{
          produces<pat::CompositeCandidateCollection>("Genmatch");
      	  produces<pat::CompositeCandidateCollection>("SelectedTracks");
      	  produces<pat::CompositeCandidateCollection>("SelectedDiTracks");
	  //produces<reco::VertexCompositePtrCandidateCollection>("LooseKshort");
          //produces<pat::CompositeCandidateCollection>("SelectedV0Collection");
          //produces<TransientTrackCollection>("SelectedV0TransientCollection");
	  produces<pat::CompositeCandidateCollection>("B");
  }

  ~BDhFitter_v2() override {}

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  const edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTracksToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;

  const int tkNHitsCut_;
  const double tkPtCut_;
  const double tkEtaCut_;
  const double tkChi2Cut_;

  const double vtxChi2Cut_;
  const double vtxDecaySigXYCut_;
  const double vtxDecaySigXYZCut_;
  const double cosThetaXYCut_;
  const double cosThetaXYZCut_;

  const double mPiPiCut_;
  const double kShortMassCut_;
  const double D0MassCut_;
  const double BMassCut_;
  const bool savetrack_;
  const int verbose;

  const bool isSignalMC=true;
  const bool onlyRecoMatchedPions=false;
  const bool onlyKeepGen=false;
};

int GetGlobalIndex(const reco::Candidate *cand, const edm::Handle<reco::GenParticleCollection> &handle) {
    for (size_t j = 0; j < handle->size(); j++) {
        if (&(*handle)[j] == cand) return j;
    }
    return -1;
}

std::tuple<int, float, float, float> MatchGenPart(
    const reco::Candidate* genCand,
    std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vec_trk_ttrk)
{
    int matchedIndex = -1;
    float thematch_pt = 0;
    float thematch_dr = 99;
    float thematch_iso = 0;
    if (genCand->pt()==0) return std::make_tuple(matchedIndex, thematch_pt, thematch_iso, thematch_dr);
    for (unsigned int iSort = 0; iSort < vec_trk_ttrk.size(); ++iSort) {
        const auto& track = vec_trk_ttrk[iSort].first;
        if(genCand->pdgId()*track.charge()<=0) continue;
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

std::vector<int> BuildBDecayChain(const edm::Handle<reco::GenParticleCollection> &genParticleHandle_,
                 std::vector<pat::CompositeCandidate> *Genmatch_out,
                 std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vec_trk_ttrk
		) {
    
    std::vector<int> genpion_index;
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
                idx_bpion = GetGlobalIndex(bDau, genParticleHandle_);
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
        Genmatch.addUserFloat("ks0_flight_distance", ks0_flight_distance);
        Genmatch.addUserFloat("ks0_flight_distance_2D", ks0_flight_distance_2D);
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
            std::tie(idx_kspip_thematch, idx_kspip_pt, idx_kspip_iso03, idx_kspip_dr) = MatchGenPart(genCand, vec_trk_ttrk);
        }
        Genmatch.addUserFloat("idx_kspip_pt", idx_kspip_pt);
        Genmatch.addUserFloat("idx_kspip_iso03", idx_kspip_iso03);
        Genmatch.addUserFloat("idx_kspip_dr", idx_kspip_dr);
        if (idx_kspim >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_kspim];
            std::tie(idx_kspim_thematch, idx_kspim_pt, idx_kspim_iso03, idx_kspim_dr) = MatchGenPart(genCand, vec_trk_ttrk);
        }
        Genmatch.addUserFloat("idx_kspim_pt", idx_kspim_pt);
        Genmatch.addUserFloat("idx_kspim_iso03", idx_kspim_iso03);
        Genmatch.addUserFloat("idx_kspim_dr", idx_kspim_dr);
	if (idx_pip >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_pip];
            std::tie(idx_pip_thematch, idx_pip_pt, idx_pip_iso03, idx_pip_dr) = MatchGenPart(genCand, vec_trk_ttrk);
        }
        Genmatch.addUserFloat("idx_pip_pt", idx_pip_pt);
        Genmatch.addUserFloat("idx_pip_iso03", idx_pip_iso03);
        Genmatch.addUserFloat("idx_pip_dr", idx_pip_dr);
        if (idx_pim >= 0) {
            const reco::Candidate* genCand = &(*genParticleHandle_)[idx_pim];
            std::tie(idx_pim_thematch, idx_pim_pt, idx_pim_iso03, idx_pim_dr) = MatchGenPart(genCand, vec_trk_ttrk);
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


void BDhFitter_v2::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {

  using std::vector;

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
  if (!theBeamSpotHandle.isValid()) {
    edm::LogError("BDhFitter_v2") << "No BeamSpot found!";
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
    edm::LogError("BDhFitter_v2") << "Track collections not found!";
    return;
  }

  edm::Handle<reco::GenParticleCollection> genParticleHandle_;
  iEvent.getByToken(genParticleToken_, genParticleHandle_);
  if (!genParticleHandle_.isValid()) {
    edm::LogError("BDhFitter_v2") << "GenParticles not found!";
    return;
  }

  // Filter track, track collection
  unsigned int nTracks = tracks->size();
  unsigned int totalTracks = nTracks + lostTracks->size();
  std::vector<std::pair<pat::PackedCandidate, reco::TransientTrack>> vec_trk_ttrk;   // track collection
  std::vector<std::pair<bool, int>> vectrk_idx;  // save track index
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
    if (fabs(bestTrack->eta()) > tkEtaCut_) continue;
    if (fabs(trk.pdgId()) != 211) continue; // pions only
    if (iTrk < nTracks && !trk.trackHighPurity()) continue;
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
    vec_trk_ttrk.emplace_back(trk, std::move(tmpTransient));
    if (iTrk < nTracks) {
        vectrk_idx.push_back(std::make_pair(true, iTrk));
    } else {
        vectrk_idx.push_back(std::make_pair(false, iTrk - nTracks));
    }
    if(savetrack_){
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
  }


  std::unique_ptr<pat::CompositeCandidateCollection> Genmatch_out(new pat::CompositeCandidateCollection());
  // Gen Study
  std::vector<int> genpion_index;
  if(isSignalMC){
      genpion_index = BuildBDecayChain(genParticleHandle_, Genmatch_out.get(), vec_trk_ttrk);
  }

  // Get general di-track collection 
  // could from Ks, could from D0
  std::vector<std::pair<unsigned int, unsigned int>> vec_ditrk;
  std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double>> vec_ditrk_properties;  // cxPtR2, cxPtZ, trk1_dot_trk2, ...

  std::unique_ptr<pat::CompositeCandidateCollection> ditracks_earlyout(new pat::CompositeCandidateCollection());
  for (unsigned int ditrkidx1 = 0; ditrkidx1 < vec_trk_ttrk.size(); ++ditrkidx1) {
    const auto& diTrk1 = vec_trk_ttrk[ditrkidx1].first;
    const auto& tdiTrk1 = vec_trk_ttrk[ditrkidx1].second;
    //if(onlyRecoMatchedPions && isSignalMC){
    //    if(int(ditrkidx1)!= genpion_index.at(0) && int(ditrkidx1)!= genpion_index.at(1) && int(ditrkidx1)!= genpion_index.at(2) && int(ditrkidx1)!= genpion_index.at(3)){
    //      continue;
    //    }
    //}
    for (unsigned int ditrkidx2 = ditrkidx1 + 1; ditrkidx2 < vec_trk_ttrk.size(); ++ditrkidx2) {
      const auto& diTrk2 = vec_trk_ttrk[ditrkidx2].first;
      const auto& tdiTrk2 = vec_trk_ttrk[ditrkidx2].second;
      //if(onlyRecoMatchedPions && isSignalMC){
      //    if(int(ditrkidx2)!= genpion_index.at(0) && int(ditrkidx2)!= genpion_index.at(1) && int(ditrkidx2)!= genpion_index.at(2) && int(ditrkidx2)!= genpion_index.at(3)){
      //      continue;
      //    }
      //}

      // Opposite charge requirement
      if (diTrk1.charge() * diTrk2.charge() >= 0) continue;
      // reorder, (tPosTrk, tNegTrk)
      const auto& tPosTrk = (diTrk1.charge() > 0) ? tdiTrk1 : tdiTrk2;
      const auto& tNegTrk = (diTrk1.charge() < 0) ? tdiTrk1 : tdiTrk2;
      unsigned int ditrkidx_pos = (diTrk1.charge() > 0) ? ditrkidx1 : ditrkidx2;
      unsigned int ditrkidx_neg = (diTrk1.charge() < 0) ? ditrkidx1 : ditrkidx2;

      // step1: DCA requirement before anything else
      // dca, cxPt, posTSCP, negTSCP
      const auto& posImpact = tPosTrk.impactPointTSCP();
      const auto& negImpact = tNegTrk.impactPointTSCP();
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
      if (cxPtR2 > 120. * 120. || std::abs(cxPt.z()) > 300.) continue;
      TrajectoryStateClosestToPoint posTSCP = tPosTrk.trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint negTSCP = tNegTrk.trajectoryStateClosestToPoint(cxPt);
      if (!posTSCP.isValid() || !negTSCP.isValid())          continue;
      float trk1_dot_trk2 = posTSCP.momentum().dot(negTSCP.momentum());
      if (trk1_dot_trk2 < 0.1 ) continue;
      //distance closest approach in x,y wrt beam spot
      std::pair<double, double> DCA_pos_beamspot = computeDCA(tPosTrk, *theBeamSpot);
      std::pair<double, double> DCA_neg_beamspot = computeDCA(tNegTrk, *theBeamSpot);

      // step2: Pre-fit mPiPi mass cut
      double p1Mag2 = posTSCP.momentum().mag2();
      double p2Mag2 = negTSCP.momentum().mag2();
      double totalE = std::sqrt(p1Mag2 + piMassSquared) + std::sqrt(p2Mag2 + piMassSquared);
      double totalESq = totalE * totalE;
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
      double massSquared = totalESq - totalPSq;

      // step3: Vertex Fit!
      // Create dummy covariance (error matrix) for the initial vertex position
      GlobalError dummyError(1.0e-3, 0.0, 1.0e-3, 0.0, 0.0, 1.0e-3);
      std::vector<reco::TransientTrack> transTracks{tPosTrk, tNegTrk};
      TransientVertex theRecoVertex(cxPt, dummyError, transTracks, 1.0e-3);
      KalmanVertexFitter theKalmanFitter(true);
      theRecoVertex = theKalmanFitter.vertex(transTracks);
      if (!theRecoVertex.isValid()) continue;
      reco::Vertex theVtx_KLM = reco::Vertex(theRecoVertex);
      if (theVtx_KLM.normalizedChi2() > vtxChi2Cut_) continue;
      // Get fitted vertex position
      GlobalPoint vtxPos(theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z());
      if(verbose>=1){
          std::cout << "A valid vertex found at ("
                << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z() << ")"
                << std::endl;
      }

      // step4: check the fit vertex
      // 2D decay significance
      SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx_KLM.covariance();
      SVector3 distVecXY(vtxPos.x() - theBeamSpotPos.x(), vtxPos.y() - theBeamSpotPos.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      double distMagXYErr = sqrt(ROOT::Math::Similarity(totalCov, distVecXY));
      // trajectory
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
            std::make_unique<TrajectoryStateClosestToPoint>(tdiTrk1.trajectoryStateClosestToPoint(vtxPos));
        trajMins =
            std::make_unique<TrajectoryStateClosestToPoint>(tdiTrk2.trajectoryStateClosestToPoint(vtxPos));
      }
      if (trajPlus.get() == nullptr || trajMins.get() == nullptr || !trajPlus->isValid() || !trajMins->isValid()) continue;
      GlobalVector positiveP(trajPlus->momentum());
      GlobalVector negativeP(trajMins->momentum());
      GlobalVector totalP(positiveP + negativeP);
      // 2D pointing angle
      double dx = theVtx_KLM.x() - theBeamSpotPos.x();
      double dy = theVtx_KLM.y() - theBeamSpotPos.y();
      double px = totalP.x();
      double py = totalP.y();
      double angleXY = (dx * px + dy * py) / (sqrt(dx * dx + dy * dy) * sqrt(px * px + py * py));
      if (angleXY < cosThetaXYCut_) continue;

      vec_ditrk.push_back(std::make_pair(ditrkidx_pos, ditrkidx_neg));
      vec_ditrk_properties.emplace_back(cxPtR2, cxPt.z(), trk1_dot_trk2, DCA_pos_beamspot.first, DCA_pos_beamspot.second, DCA_neg_beamspot.first, DCA_neg_beamspot.second, dca, massSquared, theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z(), theVtx_KLM.chi2(), theVtx_KLM.ndof(), theVtx_KLM.normalizedChi2(), distMagXY, distMagXYErr, angleXY);

      if(savetrack_){
          pat::CompositeCandidate pcand;
          pcand.addUserInt("leg1_idx",   int(ditrkidx_pos));
          pcand.addUserInt("leg2_idx",   int(ditrkidx_neg));
          ditracks_earlyout -> emplace_back(pcand);
      }

    }
  }
  if(verbose>=2) std::cout<<"vec_trk_ttrk.size() "<< vec_trk_ttrk.size()<<" vec_ditrk.size() "<< vec_ditrk.size()<<std::endl;


  // decide to save or not
  auto LooseKshorts_out = std::make_unique<reco::VertexCompositePtrCandidateCollection>();
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  std::unique_ptr<TransientTrackCollection> trans_out( new TransientTrackCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> Bu_out(new pat::CompositeCandidateCollection());

  // If early stop
  if( (isSignalMC && (genpion_index.size()<4 || (genpion_index.size()>=4 && genpion_index.at(3)==-1))) || (onlyKeepGen) ){
      iEvent.put(std::move(Genmatch_out), "Genmatch");
      if(savetrack_) iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");
      if(savetrack_) iEvent.put(std::move(ditracks_earlyout),       "SelectedDiTracks");
//      iEvent.put(std::move(LooseKshorts_out), "LooseKshort");
//      iEvent.put(std::move(ret_val), "SelectedV0Collection");
//      iEvent.put(std::move(trans_out), "SelectedV0TransientCollection");
      iEvent.put(std::move(Bu_out), "B");
      return;
  }




  // Two ditracks make one D0
  // vec_ditrk
  // Add one more track, make a B
  // vec_trk_ttrk
  for (unsigned int ditrkidx1 = 0; ditrkidx1 < vec_ditrk.size(); ++ditrkidx1) { 
      if(onlyRecoMatchedPions && isSignalMC){
          if(int(vec_ditrk[ditrkidx1].first)!= genpion_index.at(0) && int(vec_ditrk[ditrkidx1].first)!= genpion_index.at(1) && int(vec_ditrk[ditrkidx1].second)!= genpion_index.at(0) && int(vec_ditrk[ditrkidx1].second)!= genpion_index.at(1)){
            continue;
          }
      }
      if(std::get<4>(vec_ditrk_properties[ditrkidx1])==0 || std::get<6>(vec_ditrk_properties[ditrkidx1])==0) continue; 
      double cxPtR2_idx1               = std::get<0>(vec_ditrk_properties[ditrkidx1]);
      double cxPtR2trk_pos_dcasig_idx1 = std::get<3>(vec_ditrk_properties[ditrkidx1])/std::get<4>(vec_ditrk_properties[ditrkidx1]);
      double cxPtR2trk_neg_dcasig_idx1 = std::get<5>(vec_ditrk_properties[ditrkidx1])/std::get<6>(vec_ditrk_properties[ditrkidx1]);
      double massSquared_idx1          = std::get<8>(vec_ditrk_properties[ditrkidx1]);
      if(verbose>=2) std::cout<<"cxPtR2, cxPtR2trk_pos_dcasig, cxPtR2trk_neg_dcasig, massSquared: "<<cxPtR2_idx1<<" "<< cxPtR2trk_pos_dcasig_idx1<<" "<< cxPtR2trk_neg_dcasig_idx1<<" "<< massSquared_idx1<< std::endl;
      // reduce potential candidates to ~20, dedicated for Ks0
      if(massSquared_idx1 > mPiPiCut_*mPiPiCut_ ) continue;
      if (abs(cxPtR2trk_pos_dcasig_idx1)<2 || abs(cxPtR2trk_neg_dcasig_idx1)<2) continue;
      if(cxPtR2_idx1<0.02) continue;
      if(verbose>=2) std::cout<< "Pass Ks0 presel" << std::endl;

      // Kin fit for Ks0
      const auto& tTrk1 = vec_trk_ttrk[vec_ditrk[ditrkidx1].first].second;
      const auto& tTrk2 = vec_trk_ttrk[vec_ditrk[ditrkidx1].second].second;
      KinVtxFitter theKinVtxFitter(
      {tTrk1, tTrk2},
      {piMass, piMass},
      {K_SIGMA, K_SIGMA} );
      if ( !theKinVtxFitter.success() ) continue;
      if ( theKinVtxFitter.chi2()<0 || theKinVtxFitter.dof()<0) continue;
      if ( theKinVtxFitter.prob()<0.001) continue;
      auto theKinVtxFitter_p4 = theKinVtxFitter.fitted_p4();
      if(abs(theKinVtxFitter_p4.mass()-kShortMass)>kShortMassCut_) continue;
      if(abs(theKinVtxFitter.fitted_p4().eta())>2.4)  continue;
      if(verbose>=2) std::cout<< "Pass Kin fit" << std::endl;

      float Kin_mass = theKinVtxFitter.fitted_candidate().mass();
      float Kin_massErr = sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6));
      auto V0TT = theKinVtxFitter.fitted_candidate_ttrk();

      for (unsigned int ditrkidx2 = 0; ditrkidx2 < vec_ditrk.size(); ++ditrkidx2) {
	  //std::cout<<"Mass before fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
	  if(ditrkidx1==ditrkidx2) continue;
	  if(vec_ditrk[ditrkidx1].first == vec_ditrk[ditrkidx2].first || vec_ditrk[ditrkidx1].first == vec_ditrk[ditrkidx2].second || vec_ditrk[ditrkidx1].second == vec_ditrk[ditrkidx2].first || vec_ditrk[ditrkidx1].second == vec_ditrk[ditrkidx2].second) continue;

          if(std::get<4>(vec_ditrk_properties[ditrkidx2])==0 || std::get<6>(vec_ditrk_properties[ditrkidx2])==0) continue;
          double cxPtR2_idx2               = std::get<0>(vec_ditrk_properties[ditrkidx2]);
          double dca_idx2          = std::get<7>(vec_ditrk_properties[ditrkidx2]);
          // reduce potential candidates to ~20, dedicated for Ks0
          if(cxPtR2_idx2>6.8) continue; // layer 2 of BPIX
          if(dca_idx2>1) continue;

	  // A quick Mass fit
          const auto& tTrk3 = vec_trk_ttrk[vec_ditrk[ditrkidx2].first].second;
          const auto& tTrk4 = vec_trk_ttrk[vec_ditrk[ditrkidx2].second].second;
	  auto mom1 = vec_trk_ttrk[vec_ditrk[ditrkidx2].first].first.momentum();
	  auto mom2 = vec_trk_ttrk[vec_ditrk[ditrkidx2].second].first.momentum();
          reco::Candidate::LorentzVector p4_1(mom1.x(), mom1.y(), mom1.z(), sqrt(mom1.mag2() + piMass * piMass));
          reco::Candidate::LorentzVector p4_2(mom2.x(), mom2.y(), mom2.z(), sqrt(mom2.mag2() + piMass * piMass));
          reco::Candidate::LorentzVector total_p4 = p4_1 + p4_2 + theKinVtxFitter_p4;
          float D0_premass = total_p4.mass();
          if(abs(D0_premass-D0Mass)>D0MassCut_*2) continue;
          if(total_p4.pt()<0.3) continue;
          if(abs(total_p4.eta())>2.4) continue;

	  // Kin fit for D0
	  KinVtxFitter D0_KinFitter(
                  {tTrk3, tTrk4, V0TT },
                  {piMass, piMass,  Kin_mass},
                  {K_SIGMA, K_SIGMA, Kin_massErr}
		  );
          //std::cout<<"Mass after fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
	  if ( !D0_KinFitter.success() ) continue;
          if ( D0_KinFitter.chi2()<0 || D0_KinFitter.dof()<0) continue;
          if ( D0_KinFitter.prob()<0.00001) continue;
          auto D0_Kinfit_p4 = D0_KinFitter.fitted_p4();
          if(abs(D0_Kinfit_p4.mass()-D0Mass)>D0MassCut_) continue;    
	  if(verbose>=2) std::cout<<" D0 looks good "<<std::endl;

          auto DTT = D0_KinFitter.fitted_candidate_ttrk();
          float DKin_mass = D0_KinFitter.fitted_candidate().mass();
          float DKin_massErr = sqrt(D0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6));

/*
              // Save something
              pat::CompositeCandidate Bu;
              // idx
              Bu.addUserInt("DiTrack_idx1", int(ditrkidx1));
              Bu.addUserInt("DiTrack_idx2", int(ditrkidx2));
              Bu.addUserInt("Track_idx1", int(vec_ditrk[ditrkidx1].first));
              Bu.addUserInt("Track_idx2", int(vec_ditrk[ditrkidx1].second));
              Bu.addUserInt("Track_idx3", int(vec_ditrk[ditrkidx2].first));
              Bu.addUserInt("Track_idx4", int(vec_ditrk[ditrkidx2].second));
              // Ditrack variables
              Bu.addUserFloat("DiTrk1_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_dot",         std::get<2>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_trk1_dcaSig", std::get<3>(vec_ditrk_properties[ditrkidx1])/std::get<4>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_trk2_dcaSig", std::get<5>(vec_ditrk_properties[ditrkidx1])/std::get<6>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_dca",         std::get<7>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_x",          std::get<9>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_y",          std::get<10>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_distMagXY",      std::get<15>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_distMagXYErr", std::get<16>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk2_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dot",         std::get<2>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk1_dcaSig", std::get<3>(vec_ditrk_properties[ditrkidx2])/std::get<4>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk2_dcaSig", std::get<5>(vec_ditrk_properties[ditrkidx2])/std::get<6>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dca",         std::get<7>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_x",          std::get<9>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_y",          std::get<10>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_distMagXY",      std::get<15>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_distMagXYErr", std::get<16>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx2]));
              auto lxy = l_xy(theKinVtxFitter, *theBeamSpot);
              Bu.addUserFloat("Ks0_Kin_vtx_x", theKinVtxFitter.fitted_vtx().x());
              Bu.addUserFloat("Ks0_Kin_vtx_y", theKinVtxFitter.fitted_vtx().y());
              Bu.addUserFloat("Ks0_Kin_vtx_z", theKinVtxFitter.fitted_vtx().z());
              Bu.addUserFloat("Ks0_Kin_chi2",  theKinVtxFitter.chi2());
              Bu.addUserFloat("Ks0_Kin_dof",   theKinVtxFitter.dof()); // always 1 ?
              Bu.addUserFloat("Ks0_Kin_prob",  theKinVtxFitter.prob());
              Bu.addUserFloat("Ks0_Kin_pt",    theKinVtxFitter.fitted_p4().pt());
              Bu.addUserFloat("Ks0_Kin_eta",   theKinVtxFitter.fitted_p4().eta());
              Bu.addUserFloat("Ks0_Kin_phi",   theKinVtxFitter.fitted_p4().phi());
              Bu.addUserFloat("Ks0_Kin_mass",  theKinVtxFitter.fitted_candidate().mass());
              Bu.addUserFloat("Ks0_Kin_massErr", sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("Ks0_Kin_fitted_cos_theta_2D", cos_theta_2D(theKinVtxFitter, *theBeamSpot, theKinVtxFitter.fitted_p4()));
              Bu.addUserFloat("Ks0_Kin_l_xy", lxy.value());
              Bu.addUserFloat("Ks0_Kin_l_xy_unc", lxy.error());
              Bu.addUserFloat("Ks0_Kin_trk1_pt", theKinVtxFitter.daughter_p4(0).pt());
              Bu.addUserFloat("Ks0_Kin_trk1_eta", theKinVtxFitter.daughter_p4(0).eta());
              Bu.addUserFloat("Ks0_Kin_trk1_phi", theKinVtxFitter.daughter_p4(0).phi());
              Bu.addUserFloat("Ks0_Kin_trk2_pt", theKinVtxFitter.daughter_p4(1).pt());
              Bu.addUserFloat("Ks0_Kin_trk2_eta", theKinVtxFitter.daughter_p4(1).eta());
              Bu.addUserFloat("Ks0_Kin_trk2_phi", theKinVtxFitter.daughter_p4(1).phi());
              auto D0_lxy = l_xy(D0_KinFitter, *theBeamSpot);
              Bu.addUserFloat("D0_premass", D0_premass);
              Bu.addUserFloat("D0_Kinfitmass", D0_Kinfit_p4.mass());
              Bu.addUserFloat("D0_Kin_vtx_x", D0_KinFitter.fitted_vtx().x());
              Bu.addUserFloat("D0_Kin_vtx_y", D0_KinFitter.fitted_vtx().y());
              Bu.addUserFloat("D0_Kin_vtx_z", D0_KinFitter.fitted_vtx().z());
              Bu.addUserFloat("D0_Kin_chi2", D0_KinFitter.chi2());
              Bu.addUserFloat("D0_Kin_dof", D0_KinFitter.dof());
              Bu.addUserFloat("D0_Kin_prob", D0_KinFitter.prob());
              Bu.addUserFloat("D0_Kin_pt", D0_Kinfit_p4.pt() );
              Bu.addUserFloat("D0_Kin_eta", D0_Kinfit_p4.eta() );
              Bu.addUserFloat("D0_Kin_phi", D0_Kinfit_p4.phi() );
              Bu.addUserFloat("D0_Kin_mass", D0_Kinfit_p4.mass() );
              Bu.addUserFloat("D0_Kin_massErr", sqrt(D0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("D0_Kin_fitted_cos_theta_2D", cos_theta_2D(D0_KinFitter, *theBeamSpot, D0_KinFitter.fitted_p4()));
              Bu.addUserFloat("D0_Kin_l_xy", D0_lxy.value());
              Bu.addUserFloat("D0_Kin_l_xy_unc", D0_lxy.error());
              Bu.addUserFloat("D0_Kin_trk3_pt", D0_KinFitter.daughter_p4(0).pt());
              Bu.addUserFloat("D0_Kin_trk3_eta", D0_KinFitter.daughter_p4(0).eta());
              Bu.addUserFloat("D0_Kin_trk3_phi", D0_KinFitter.daughter_p4(0).phi());
              Bu.addUserFloat("D0_Kin_trk4_pt", D0_KinFitter.daughter_p4(1).pt());
              Bu.addUserFloat("D0_Kin_trk4_eta", D0_KinFitter.daughter_p4(1).eta());
              Bu.addUserFloat("D0_Kin_trk4_phi", D0_KinFitter.daughter_p4(1).phi());
              Bu.addUserFloat("D0_Kin_ks0_pt", D0_KinFitter.daughter_p4(2).pt());
              Bu.addUserFloat("D0_Kin_ks0_eta", D0_KinFitter.daughter_p4(2).eta());
              Bu.addUserFloat("D0_Kin_ks0_phi", D0_KinFitter.daughter_p4(2).phi());
              reco::Candidate::Point ks0_decay_vtx(Bu.userFloat("Ks0_Kin_vtx_x"), Bu.userFloat("Ks0_Kin_vtx_y"), Bu.userFloat("Ks0_Kin_vtx_z") );
              reco::Candidate::Point ks0_production_vtx(Bu.userFloat("D0_Kin_vtx_x"), Bu.userFloat("D0_Kin_vtx_y"), Bu.userFloat("D0_Kin_vtx_z"));
              math::XYZVector ks0_flight = ks0_decay_vtx - ks0_production_vtx;
              float ks0_flight_distance = ks0_flight.R();
              float ks0_flight_distance_2D = sqrt(ks0_flight.x() * ks0_flight.x() + ks0_flight.y() * ks0_flight.y());
              Bu.addUserFloat("ks0_flight_distance", ks0_flight_distance);
              Bu.addUserFloat("ks0_flight_distance_2D", ks0_flight_distance_2D);
              Bu_out->push_back(Bu);
*/

	  // Make a B
          for (unsigned int Btrack_idx = 0; Btrack_idx < vec_trk_ttrk.size(); ++Btrack_idx) {
              if(Btrack_idx == vec_ditrk[ditrkidx2].first || Btrack_idx == vec_ditrk[ditrkidx2].second || Btrack_idx == vec_ditrk[ditrkidx1].second || Btrack_idx == vec_ditrk[ditrkidx1].first) continue;
              const auto& Btrack = vec_trk_ttrk[Btrack_idx].first;
              const auto& tBtrack = vec_trk_ttrk[Btrack_idx].second;
              auto mom3 = vec_trk_ttrk[Btrack_idx].first.momentum();

              // A quick Mass fit
              reco::Candidate::LorentzVector p4_3(mom3.x(), mom3.y(), mom3.z(), sqrt(mom3.mag2() + piMass * piMass));
              reco::Candidate::LorentzVector B_pre_p4 = total_p4 + p4_3;
              float B_premass = B_pre_p4.mass();
              if(abs(B_premass-BuMass)>BMassCut_*2) continue;

	      // kin fit for B
              KinVtxFitter B_KinFitter(
                      {tBtrack, DTT },
                      {piMass, DKin_mass},
                      {K_SIGMA, DKin_massErr}
                      );
              if ( !B_KinFitter.success() ) continue;
              if ( B_KinFitter.chi2()<0 || B_KinFitter.dof()<0) continue;
              if ( B_KinFitter.prob()<0.00001) continue;
              auto B_Kinfit_p4 = B_KinFitter.fitted_p4();
              if(abs(B_Kinfit_p4.mass()-BuMass)>BMassCut_) continue;


              // Save something
              pat::CompositeCandidate Bu;
	      // idx
	      Bu.addUserInt("DiTrack_idx1", int(ditrkidx1));
              Bu.addUserInt("DiTrack_idx2", int(ditrkidx2));
              Bu.addUserInt("Track_idx1", int(vec_ditrk[ditrkidx1].first));
              Bu.addUserInt("Track_idx2", int(vec_ditrk[ditrkidx1].second));
              Bu.addUserInt("Track_idx3", int(vec_ditrk[ditrkidx2].first));
              Bu.addUserInt("Track_idx4", int(vec_ditrk[ditrkidx2].second));
              Bu.addUserInt("BTrack_idx", int(Btrack_idx));
	      // Ditrack variables
              Bu.addUserFloat("DiTrk1_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx1]));		      
	      Bu.addUserFloat("DiTrk1_dot",         std::get<2>(vec_ditrk_properties[ditrkidx1]));
	      Bu.addUserFloat("DiTrk1_trk1_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrkidx1]))/std::get<4>(vec_ditrk_properties[ditrkidx1]));
	      Bu.addUserFloat("DiTrk1_trk2_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrkidx1]))/std::get<6>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_dca",         std::get<7>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx1]));
	      Bu.addUserFloat("DiTrk1_KLM_vtx_x",          std::get<9>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_y",          std::get<10>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_distMagXY",      std::get<15>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_distMagXYErr", std::get<16>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk2_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dot",         std::get<2>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk1_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrkidx2]))/std::get<4>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk2_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrkidx2]))/std::get<6>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dca",         std::get<7>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_x",          std::get<9>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_y",          std::get<10>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_distMagXY",      std::get<15>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_distMagXYErr", std::get<16>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx2]));
	      // Ks0
	      auto lxy = l_xy(theKinVtxFitter, *theBeamSpot);
              Bu.addUserFloat("Ks0_Kin_vtx_x", theKinVtxFitter.fitted_vtx().x());
              Bu.addUserFloat("Ks0_Kin_vtx_y", theKinVtxFitter.fitted_vtx().y());
              Bu.addUserFloat("Ks0_Kin_vtx_z", theKinVtxFitter.fitted_vtx().z());
              Bu.addUserFloat("Ks0_Kin_chi2",  theKinVtxFitter.chi2());
              Bu.addUserFloat("Ks0_Kin_dof",   theKinVtxFitter.dof()); // always 1 ?
              Bu.addUserFloat("Ks0_Kin_prob",  theKinVtxFitter.prob());
              Bu.addUserFloat("Ks0_Kin_pt",    theKinVtxFitter.fitted_p4().pt());
              Bu.addUserFloat("Ks0_Kin_eta",   theKinVtxFitter.fitted_p4().eta());
              Bu.addUserFloat("Ks0_Kin_phi",   theKinVtxFitter.fitted_p4().phi());
              Bu.addUserFloat("Ks0_Kin_mass",  theKinVtxFitter.fitted_candidate().mass());
              Bu.addUserFloat("Ks0_Kin_massErr", sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("Ks0_Kin_fitted_cos_theta_2D", cos_theta_2D(theKinVtxFitter, *theBeamSpot, theKinVtxFitter.fitted_p4()));
              Bu.addUserFloat("Ks0_Kin_l_xy", lxy.value());
              Bu.addUserFloat("Ks0_Kin_l_xy_unc", lxy.error());
              Bu.addUserFloat("Ks0_Kin_trk1_pt", theKinVtxFitter.daughter_p4(0).pt());
              Bu.addUserFloat("Ks0_Kin_trk1_eta", theKinVtxFitter.daughter_p4(0).eta());
              Bu.addUserFloat("Ks0_Kin_trk1_phi", theKinVtxFitter.daughter_p4(0).phi());
              Bu.addUserFloat("Ks0_Kin_trk2_pt", theKinVtxFitter.daughter_p4(1).pt());
              Bu.addUserFloat("Ks0_Kin_trk2_eta", theKinVtxFitter.daughter_p4(1).eta());
              Bu.addUserFloat("Ks0_Kin_trk2_phi", theKinVtxFitter.daughter_p4(1).phi());

	      // D0
              auto D0_lxy = l_xy(D0_KinFitter, *theBeamSpot);
              Bu.addUserFloat("D0_premass", D0_premass);
              Bu.addUserFloat("D0_Kin_vtx_x", D0_KinFitter.fitted_vtx().x());
              Bu.addUserFloat("D0_Kin_vtx_y", D0_KinFitter.fitted_vtx().y());
              Bu.addUserFloat("D0_Kin_vtx_z", D0_KinFitter.fitted_vtx().z());
              Bu.addUserFloat("D0_Kin_chi2", D0_KinFitter.chi2());
              Bu.addUserFloat("D0_Kin_dof", D0_KinFitter.dof());
              Bu.addUserFloat("D0_Kin_prob", D0_KinFitter.prob());
              Bu.addUserFloat("D0_Kin_pt", D0_Kinfit_p4.pt() );
              Bu.addUserFloat("D0_Kin_eta", D0_Kinfit_p4.eta() );
              Bu.addUserFloat("D0_Kin_phi", D0_Kinfit_p4.phi() );
              Bu.addUserFloat("D0_Kin_mass", D0_Kinfit_p4.mass() );
              Bu.addUserFloat("D0_Kin_massErr", sqrt(D0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("D0_Kin_fitted_cos_theta_2D", cos_theta_2D(D0_KinFitter, *theBeamSpot, D0_KinFitter.fitted_p4()));
              Bu.addUserFloat("D0_Kin_l_xy", D0_lxy.value());
              Bu.addUserFloat("D0_Kin_l_xy_unc", D0_lxy.error());
	      Bu.addUserFloat("D0_Kin_trk3_pt", D0_KinFitter.daughter_p4(0).pt());
              Bu.addUserFloat("D0_Kin_trk3_eta", D0_KinFitter.daughter_p4(0).eta());
              Bu.addUserFloat("D0_Kin_trk3_phi", D0_KinFitter.daughter_p4(0).phi());
              Bu.addUserFloat("D0_Kin_trk4_pt", D0_KinFitter.daughter_p4(1).pt());
              Bu.addUserFloat("D0_Kin_trk4_eta", D0_KinFitter.daughter_p4(1).eta());
              Bu.addUserFloat("D0_Kin_trk4_phi", D0_KinFitter.daughter_p4(1).phi());
              Bu.addUserFloat("D0_Kin_ks0_pt", D0_KinFitter.daughter_p4(2).pt());
              Bu.addUserFloat("D0_Kin_ks0_eta", D0_KinFitter.daughter_p4(2).eta());
              Bu.addUserFloat("D0_Kin_ks0_phi", D0_KinFitter.daughter_p4(2).phi());
              reco::Candidate::Point ks0_decay_vtx(Bu.userFloat("Ks0_Kin_vtx_x"), Bu.userFloat("Ks0_Kin_vtx_y"), Bu.userFloat("Ks0_Kin_vtx_z") );
              reco::Candidate::Point ks0_production_vtx(Bu.userFloat("D0_Kin_vtx_x"), Bu.userFloat("D0_Kin_vtx_y"), Bu.userFloat("D0_Kin_vtx_z"));
              math::XYZVector ks0_flight = ks0_decay_vtx - ks0_production_vtx;
              float ks0_flight_distance = ks0_flight.R();
              float ks0_flight_distance_2D = sqrt(ks0_flight.x() * ks0_flight.x() + ks0_flight.y() * ks0_flight.y());
              Bu.addUserFloat("ks0_flight_distance", ks0_flight_distance);
              Bu.addUserFloat("ks0_flight_distance_2D", ks0_flight_distance_2D);	      
	      
	      // B
	      auto B_lxy = l_xy(B_KinFitter, *theBeamSpot);
	      Bu.addUserFloat("B_premass", B_premass);
              Bu.addUserFloat("B_Kin_vtx_x", B_KinFitter.fitted_vtx().x());
              Bu.addUserFloat("B_Kin_vtx_y", B_KinFitter.fitted_vtx().y());
              Bu.addUserFloat("B_Kin_vtx_z", B_KinFitter.fitted_vtx().z());
              Bu.addUserFloat("B_Kin_chi2", B_KinFitter.chi2());
              Bu.addUserFloat("B_Kin_dof", B_KinFitter.dof());
              Bu.addUserFloat("B_Kin_prob", B_KinFitter.prob());
              Bu.addUserFloat("B_Kin_pt", B_Kinfit_p4.pt() );
              Bu.addUserFloat("B_Kin_eta", B_Kinfit_p4.eta() );
              Bu.addUserFloat("B_Kin_phi", B_Kinfit_p4.phi() );
              Bu.addUserFloat("B_Kin_mass", B_Kinfit_p4.mass() );
              Bu.addUserFloat("B_Kin_massErr", sqrt(B_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("B_Kin_fitted_cos_theta_2D", cos_theta_2D(B_KinFitter, *theBeamSpot, B_KinFitter.fitted_p4()));
              Bu.addUserFloat("B_Kin_l_xy", B_lxy.value());
              Bu.addUserFloat("B_Kin_l_xy_unc", B_lxy.error());
              Bu.addUserFloat("B_Kin_trk_charge", Btrack.charge());
              Bu.addUserFloat("B_Kin_trk_pt", B_KinFitter.daughter_p4(0).pt());
              Bu.addUserFloat("B_Kin_trk_eta", B_KinFitter.daughter_p4(0).eta());
              Bu.addUserFloat("B_Kin_trk_phi", B_KinFitter.daughter_p4(0).phi());
              Bu.addUserFloat("B_Kin_D0_pt", B_KinFitter.daughter_p4(1).pt());
              Bu.addUserFloat("B_Kin_D0_eta", B_KinFitter.daughter_p4(1).eta());
              Bu.addUserFloat("B_Kin_D0_phi", B_KinFitter.daughter_p4(1).phi());
              reco::Candidate::Point D0_decay_vtx(Bu.userFloat("D0_Kin_vtx_x"), Bu.userFloat("D0_Kin_vtx_y"), Bu.userFloat("D0_Kin_vtx_z"));
              reco::Candidate::Point D0_production_vtx(Bu.userFloat("B_Kin_vtx_x"), Bu.userFloat("B_Kin_vtx_y"), Bu.userFloat("B_Kin_vtx_z"));
              math::XYZVector D0_flight = D0_decay_vtx - D0_production_vtx;
              float D0_flight_distance = D0_flight.R();
              float D0_flight_distance_2D = sqrt(D0_flight.x() * D0_flight.x() + D0_flight.y() * D0_flight.y());
              Bu.addUserFloat("D0_flight_distance", D0_flight_distance);
              Bu.addUserFloat("D0_flight_distance_2D", D0_flight_distance_2D);
	      Bu_out->push_back(Bu);
	  }
  
      }
  }

  if(isSignalMC)  iEvent.put(std::move(Genmatch_out), "Genmatch");
  if(savetrack_) iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");    
  if(savetrack_) iEvent.put(std::move(ditracks_earlyout),       "SelectedDiTracks");
//  iEvent.put(std::move(LooseKshorts_out), "LooseKshort");
//  iEvent.put(std::move(ret_val), "SelectedV0Collection");
//  iEvent.put(std::move(trans_out), "SelectedV0TransientCollection");
  iEvent.put(std::move(Bu_out), "B");
  return;

}

DEFINE_FWK_MODULE(BDhFitter_v2);

