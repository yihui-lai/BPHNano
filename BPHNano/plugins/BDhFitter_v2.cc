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

#include "PhysicsTools/BPHNano/plugins/XGBooster.h"

// pdg mass constants
namespace {
  const double piMass = 0.13957018;
  const double piMassSquared = piMass * piMass;
  const double kShortMass = 0.497614;
  const double kShortMassSquared = kShortMass * kShortMass;
  const double kplusMass = 0.493677;
  const double D0Mass = 1.86484;
  const float  pion_sigma = piMass*1.e-6;
  const double BuMass = 5.27934;

}  // namespace

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

//class BDhFitter_v2 : public edm::global::EDProducer<> {
class BDhFitter_v2 : public edm::stream::EDProducer<> {

public:
  explicit BDhFitter_v2(const edm::ParameterSet &theParameters):
      bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      xgboost_models_( theParameters.getParameter<std::vector<std::string>>( "xgboost_models" ) ),
      xgboost_variable_names_( theParameters.getParameter<std::vector<std::string>>( "xgboost_variable_names" ) ),
      token_beamSpot(consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"))),
      token_vertices(consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("tracks"))),
      lostTracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("lostTracks"))),
      tkNHitsCut_( theParameters.getParameter<int>("tkNHitsCut")),
      tkPtCut_( theParameters.getParameter<double>("tkPtCut")),
      DtkPtCut_( theParameters.getParameter<double>("DtkPtCut")),
      BtkPtCut_( theParameters.getParameter<double>("BtkPtCut")),
      tkEtaCut_( theParameters.getParameter<double>("tkEtaCut")),
      tkChi2Cut_( theParameters.getParameter<double>("tkChi2Cut")),
      vtxChi2Cut_( theParameters.getParameter<double>("vtxChi2Cut")),
      vtxDecaySigXYCut_( theParameters.getParameter<double>("vtxDecaySigXYCut")),
      TrkSigXYCut_( theParameters.getParameter<double>("TrkSigXYCut")),
      vtxDecaySigXYZCut_( theParameters.getParameter<double>("vtxDecaySigXYZCut")),
      cosThetaXYCut_( theParameters.getParameter<double>("cosThetaXYCut")),
      cosThetaXYZCut_( theParameters.getParameter<double>("cosThetaXYZCut")),
      diTrack2_dca_( theParameters.getParameter<double>("diTrack2_dca")),
      mPiPiCut_( theParameters.getParameter<double>("mPiPiCut")),
      Ks0_l_xyzSigCut_( theParameters.getParameter<double>("Ks0_l_xyzSigCut")),
      D0_PtCut_( theParameters.getParameter<double>("D0_PtCut")),
      D0vtxDecaySigXYCut_( theParameters.getParameter<double>("D0vtxDecaySigXYCut")),
      B_PtCut_( theParameters.getParameter<double>("B_PtCut")),
      Btrk_dcaSigCut_( theParameters.getParameter<double>("Btrk_dcaSigCut")),
      kShortMassCut_( theParameters.getParameter<double>("kShortMassCut")),
      D0MassCut_( theParameters.getParameter<double>("D0MassCut")),
      BMassCut_( theParameters.getParameter<double>("BMassCut")),
      savetrack_( theParameters.getParameter<bool>("savetrack")),
      verbose( theParameters.getParameter<int>("verbose"))
	{
      	  produces<pat::CompositeCandidateCollection>("SelectedTracks");
      	  produces<pat::CompositeCandidateCollection>("SelectedDiTracks");
	  //produces<reco::VertexCompositePtrCandidateCollection>("LooseKshort");
          //produces<pat::CompositeCandidateCollection>("SelectedV0Collection");
          //produces<TransientTrackCollection>("SelectedV0TransientCollection");
	  produces<pat::CompositeCandidateCollection>("B");
          for (auto model: xgboost_models_)
              BMva_.push_back(XGBooster(edm::FileInPath("PhysicsTools/BPHNano/data/BDh_mva_May4/" + model + ".model").fullPath(),
                                        edm::FileInPath("PhysicsTools/BPHNano/data/BDh_mva_May4/" + model + ".features").fullPath()));

  }

  ~BDhFitter_v2() override {}
  
  //void produce(edm::StreamID, edm::Event&, const edm::EventSetup&);
  //void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  //float Bmva_estimator(pat::CompositeCandidate& Bu, unsigned int i) const;
  float Bmva_estimator(pat::CompositeCandidate& Bu, unsigned int i);

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {}

private:
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;

  std::vector<std::string> features_;
  std::vector<std::string> xgboost_models_;
  std::vector<std::string> xgboost_variable_names_;
  std::vector<XGBooster> BMva_;

  const edm::EDGetTokenT<reco::BeamSpot> token_beamSpot;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTracksToken_;
  
  const int tkNHitsCut_;
  const double tkPtCut_;
  const double DtkPtCut_;
  const double BtkPtCut_;
  const double tkEtaCut_;
  const double tkChi2Cut_;

  const double vtxChi2Cut_;
  const double vtxDecaySigXYCut_;
  const double TrkSigXYCut_;
  const double vtxDecaySigXYZCut_;
  const double cosThetaXYCut_;
  const double cosThetaXYZCut_;
  const double diTrack2_dca_;
  
  const double mPiPiCut_;
  const double Ks0_l_xyzSigCut_;
  const double D0_PtCut_;
  const double D0vtxDecaySigXYCut_;
  const double B_PtCut_;
  const double Btrk_dcaSigCut_;
  const double kShortMassCut_;
  const double D0MassCut_;
  const double BMassCut_;
  const bool savetrack_;
  const int verbose;

  const bool isSignalMC=true;
  const bool onlyRecoMatchedPions=false;
  const bool onlyKeepGen=false;
};

float BDhFitter_v2::Bmva_estimator(pat::CompositeCandidate& Bu, unsigned int i){
       BMva_.at(i).set("B_DiTrk1_dca",            Bu.userFloat("DiTrk1_dca"));
       BMva_.at(i).set("B_Ks0_Kin_vtx_r",         Bu.userFloat("Ks0_Kin_vtx_r"));
       BMva_.at(i).set("B_Ks0_Kin_d0_alpha_2D",   Bu.userFloat("Ks0_Kin_d0_alpha_2D"));
       BMva_.at(i).set("B_Ks0_Kin_d0_l_xy",       Bu.userFloat("Ks0_Kin_d0_l_xy"));
       BMva_.at(i).set("B_Ks0_Kin_d0_l_xySig",    Bu.userFloat("Ks0_Kin_d0_l_xySig"));
       BMva_.at(i).set("B_Ks0_Kin_d0_dca",        Bu.userFloat("Ks0_Kin_d0_dca"));
       BMva_.at(i).set("B_Ks0_Kin_d0_dcaSig",     Bu.userFloat("Ks0_Kin_d0_dcaSig"));
       BMva_.at(i).set("B_D0_Kin_vtx_r",          Bu.userFloat("D0_Kin_vtx_r"));
       BMva_.at(i).set("B_D0_Kin_prob",           Bu.userFloat("D0_Kin_prob"));
       BMva_.at(i).set("B_D0_Kin_b_alpha_2D",     Bu.userFloat("D0_Kin_b_alpha_2D"));
       BMva_.at(i).set("B_D0_Kin_b_l_xy",         Bu.userFloat("D0_Kin_b_l_xy"));
       BMva_.at(i).set("B_D0_Kin_b_l_xyz",        Bu.userFloat("D0_Kin_b_l_xyz"));
       BMva_.at(i).set("B_D0_Kin_b_l_xySig",      Bu.userFloat("D0_Kin_b_l_xySig"));
       BMva_.at(i).set("B_D0_Kin_b_l_xyzSig",     Bu.userFloat("D0_Kin_b_l_xyzSig"));
       BMva_.at(i).set("B_D0_Kin_b_dcaSig",       Bu.userFloat("D0_Kin_b_dcaSig"));
       BMva_.at(i).set("B_B_Kin_vtx_r",           Bu.userFloat("B_Kin_vtx_r"));
       BMva_.at(i).set("B_B_Kin_prob",            Bu.userFloat("B_Kin_prob"));
       BMva_.at(i).set("B_B_Kin_pv_alpha_2D",     Bu.userFloat("B_Kin_pv_alpha_2D"));
       BMva_.at(i).set("B_B_Kin_pv_l_xy",         Bu.userFloat("B_Kin_pv_l_xy"));
       BMva_.at(i).set("B_B_Kin_pv_l_xySig",      Bu.userFloat("B_Kin_pv_l_xySig"));
       BMva_.at(i).set("B_B_Kin_trk_b_dca",       Bu.userFloat("B_Kin_trk_b_dca"));
       BMva_.at(i).set("B_B_Kin_trk_b_dcaSig",    Bu.userFloat("B_Kin_trk_b_dcaSig"));
       BMva_.at(i).set("B_B_Kin_pt",              Bu.userFloat("B_Kin_pt"));
       BMva_.at(i).set("B_B_Kin_D0_pt",           Bu.userFloat("B_Kin_D0_pt"));
       BMva_.at(i).set("B_B_Kin_trk_pt",          Bu.userFloat("B_Kin_trk_pt"));
       return BMva_.at(i).predict();
}

//void BDhFitter_v2::produce(edm::StreamID, edm::Event &iEvent, edm::EventSetup const &iSetup) const {
void BDhFitter_v2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
  const reco::Vertex* referenceVtx = &(*vertices)[0];
  math::XYZPoint referencePos(referenceVtx->position());

  edm::Handle<edm::View<pat::PackedCandidate>> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  edm::Handle<edm::View<pat::PackedCandidate>> lostTracks;
  iEvent.getByToken(lostTracksToken_, lostTracks);
  if (!tracks.isValid() || !lostTracks.isValid()) {
    edm::LogError("BDhFitter_v2") << "Track collections not found!";
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
    //if (fabs(trk.pdgId()) != 211) continue; // pions only
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
        pcand.addUserFloat("dxySig", trk.dxy() / trk.dxyError()); // cut 2 in cmssw
        pcand.addUserFloat("dz",   trk.dz());
        pcand.addUserFloat("dzSig",  trk.dz() / trk.dzError()); // not cut
        pcand.addUserFloat("bt_pt",         bestTrack->pt());
        if (iTrk < nTracks && trk.trackHighPurity()){
            pcand.addUserFloat("trackHighPurity", 1);
        }else if (iTrk < nTracks){
            pcand.addUserFloat("trackHighPurity", 0);
        }else{
            pcand.addUserFloat("trackHighPurity", -1);
        }
        pcand.addUserFloat("dxy_bs",       ipXY_bs);
        pcand.addUserFloat("dz_bs",        ipZ_bs);
        pcand.addUserFloat("dxySig_bs",       ipsigXY_bs);
        pcand.addUserFloat("dzSig_bs",        ipsigZ_bs);
	pcand.addUserFloat("dxy_pv",       ipXY_pv);
        pcand.addUserFloat("dz_pv",        ipZ_pv);
        pcand.addUserFloat("dxySig_pv",       ipsigXY_pv);
        pcand.addUserFloat("dzSig_pv",        ipsigZ_pv);
        pcand.addUserFloat("ptErr",         bestTrack->ptError());
        pcand.addUserFloat("normChi2",      bestTrack->normalizedChi2());
        pcand.addUserInt("nValidPixelHits", bestTrack->hitPattern().numberOfValidPixelHits());
        pcand.addUserInt("nValidHits",      bestTrack->hitPattern().numberOfValidHits());
        tracks_earlyout -> emplace_back(pcand);
    }
  }

  // Get general di-track collection 
  // could from Ks, could from D0
  std::vector<std::pair<unsigned int, unsigned int>> vec_ditrk;
  std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double>> vec_ditrk_properties;

  std::unique_ptr<pat::CompositeCandidateCollection> ditracks_earlyout(new pat::CompositeCandidateCollection());
  for (unsigned int ditrkidx1 = 0; ditrkidx1 < vec_trk_ttrk.size(); ++ditrkidx1) {
    const auto& diTrk1 = vec_trk_ttrk[ditrkidx1].first;
    const auto& tdiTrk1 = vec_trk_ttrk[ditrkidx1].second;
    for (unsigned int ditrkidx2 = ditrkidx1 + 1; ditrkidx2 < vec_trk_ttrk.size(); ++ditrkidx2) {
      const auto& diTrk2 = vec_trk_ttrk[ditrkidx2].first;
      const auto& tdiTrk2 = vec_trk_ttrk[ditrkidx2].second;
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
      //the tracks should at least point in the same quadrant
      if (trk1_dot_trk2 < 0 ) continue;
      //distance closest approach in x,y wrt beam spot
      std::pair<double, double> DCA_pos_beamspot = computeDCA(tPosTrk, *theBeamSpot);
      std::pair<double, double> DCA_neg_beamspot = computeDCA(tNegTrk, *theBeamSpot);
      std::pair<double, double> DCA_pos_pv = computeDCA(tPosTrk, referencePos.x(), referencePos.y(), referencePos.z());
      std::pair<double, double> DCA_neg_pv = computeDCA(tNegTrk, referencePos.x(), referencePos.y(), referencePos.z());

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
      // Get fitted vertex position
      GlobalPoint vtxPos(theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z());
      if(verbose>=10){
          std::cout << "A valid vertex found at ("
                << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z() << ")"
                << std::endl;
      }

      // step4: check the fit vertex
      // 2D decay significance
      // relative to beamspot
      SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx_KLM.covariance();
      SVector3 distVecXY(vtxPos.x() - theBeamSpotPos.x(), vtxPos.y() - theBeamSpotPos.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      double distMagXYErr = sqrt(ROOT::Math::Similarity(totalCov, distVecXY));

      // vertex relative to PV
      SMatrixSym3D totalCov_pv = referenceVtx->covariance() + theVtx_KLM.covariance();
      SVector3 distVecXY_pv(vtxPos.x() - referencePos.x(), vtxPos.y() - referencePos.y(), 0.);
      double distMagXY_pv = ROOT::Math::Mag(distVecXY_pv);
      double distMagXYErr_pv = sqrt(ROOT::Math::Similarity(totalCov_pv, distVecXY_pv));
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
        if (thePositiveRefTrack == nullptr || theNegativeRefTrack == nullptr) continue;
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
      // 2D pointing angle relative to BeamSpot
      double dx = theVtx_KLM.x() - theBeamSpotPos.x();
      double dy = theVtx_KLM.y() - theBeamSpotPos.y();
      // 2D pointing angle relative to PV
      double dx_pv = theVtx_KLM.x() - referencePos.x();
      double dy_pv = theVtx_KLM.y() - referencePos.y();
      double px = totalP.x();
      double py = totalP.y();
      double angleXY = ((dx * px + dy * py) / (sqrt(dx * dx + dy * dy) * sqrt(px * px + py * py)));
      double angleXY_pv = ((dx_pv * px + dy_pv * py) / (sqrt(dx_pv * dx_pv + dy_pv * dy_pv) * sqrt(px * px + py * py)));

      vec_ditrk.push_back(std::make_pair(ditrkidx_pos, ditrkidx_neg));
      vec_ditrk_properties.emplace_back(cxPtR2, cxPt.z(), trk1_dot_trk2, DCA_pos_beamspot.first, DCA_pos_beamspot.second, DCA_neg_beamspot.first, DCA_neg_beamspot.second, dca, massSquared, theVtx_KLM.x(), theVtx_KLM.y(), theVtx_KLM.z(), theVtx_KLM.chi2(), theVtx_KLM.ndof(), theVtx_KLM.normalizedChi2(), distMagXY, distMagXYErr, angleXY, DCA_pos_pv.first, DCA_pos_pv.second, DCA_neg_pv.first, DCA_neg_pv.second, distMagXY_pv, distMagXYErr_pv, angleXY_pv);
      if(savetrack_){
          pat::CompositeCandidate pcand;
          pcand.addUserInt("leg1_idx",   int(ditrkidx_pos));
          pcand.addUserInt("leg2_idx",   int(ditrkidx_neg));
          ditracks_earlyout -> emplace_back(pcand);
      }

    }
  }
  if(verbose>=5) std::cout<<"vec_trk_ttrk.size() "<< vec_trk_ttrk.size()<<" vec_ditrk.size() "<< vec_ditrk.size()<<std::endl;


  // decide to save or not
  auto LooseKshorts_out = std::make_unique<reco::VertexCompositePtrCandidateCollection>();
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  std::unique_ptr<TransientTrackCollection> trans_out( new TransientTrackCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> Bu_out(new pat::CompositeCandidateCollection());

  // Two ditracks make one D0
  // vec_ditrk
  // Add one more track, make a B
  // vec_trk_ttrk
  for (unsigned int ditrkidx1 = 0; ditrkidx1 < vec_ditrk.size(); ++ditrkidx1) { 
      if(std::get<4>(vec_ditrk_properties[ditrkidx1])==0 || std::get<6>(vec_ditrk_properties[ditrkidx1])==0) continue; 
      
      if(std::get<23>(vec_ditrk_properties[ditrkidx1])==0) continue;
      if(std::get<21>(vec_ditrk_properties[ditrkidx1])==0) continue;
      if(std::get<19>(vec_ditrk_properties[ditrkidx1])==0) continue;
      double massSquared_idx1      = std::get<8>(vec_ditrk_properties[ditrkidx1]);
      double normchi2_idx1         = std::get<14>(vec_ditrk_properties[ditrkidx1]); 
      double vtxDecaySigXY_pv_idx1 = abs(std::get<22>(vec_ditrk_properties[ditrkidx1])/std::get<23>(vec_ditrk_properties[ditrkidx1])); 
      double angleXY_pv_idx1       = std::get<24>(vec_ditrk_properties[ditrkidx1]);
      double Trk1SigXY_pv_idx1     = abs(std::get<18>(vec_ditrk_properties[ditrkidx1]))/std::get<19>(vec_ditrk_properties[ditrkidx1]);
      double Trk2SigXY_pv_idx1     = abs(std::get<20>(vec_ditrk_properties[ditrkidx1]))/std::get<21>(vec_ditrk_properties[ditrkidx1]);

      if(massSquared_idx1 > mPiPiCut_*mPiPiCut_ ) continue;
      if(normchi2_idx1>vtxChi2Cut_) continue;
      if(vtxDecaySigXY_pv_idx1<vtxDecaySigXYCut_) continue;
      if(Trk1SigXY_pv_idx1<TrkSigXYCut_ || Trk2SigXY_pv_idx1<TrkSigXYCut_) continue;
      //if(vtxDecaySigXYZ_pv_idx1<vtxDecaySigXYZCut) continue;
      if (angleXY_pv_idx1 < cosThetaXYCut_) continue;
      //if (angleXYZ_pv_idx1 < cosThetaXYZCut_) continue;
      if(verbose>=3) std::cout<< "Pass Ks0 presel" << std::endl;
      
      // Kin fit for Ks0
      const auto& tTrk1 = vec_trk_ttrk[vec_ditrk[ditrkidx1].first].second;
      const auto& tTrk2 = vec_trk_ttrk[vec_ditrk[ditrkidx1].second].second;
      KinVtxFitter theKinVtxFitter(
      {tTrk1, tTrk2},
      {piMass, piMass},
      {PI_SIGMA, PI_SIGMA} );
      if ( !theKinVtxFitter.success() ) continue;
      if ( theKinVtxFitter.chi2()<0 || theKinVtxFitter.dof()<0) continue;
      if ( theKinVtxFitter.prob()<0.001) continue;
      auto theKinVtxFitter_p4 = theKinVtxFitter.fitted_p4();
      if(abs(theKinVtxFitter_p4.mass()-kShortMass)>kShortMassCut_) continue;
      if(abs(theKinVtxFitter.fitted_p4().eta())>2.4)  continue;
      if(verbose>=2) std::cout<< "Pass Ks0 Kin fit" << std::endl;

      float Kin_mass = theKinVtxFitter.fitted_candidate().mass();
      float Kin_massErr = sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6));
      auto V0TT = theKinVtxFitter.fitted_candidate_ttrk();

      for (unsigned int ditrkidx2 = 0; ditrkidx2 < vec_ditrk.size(); ++ditrkidx2) {
	  //std::cout<<"Mass before fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
	  if(ditrkidx1==ditrkidx2) continue;
	  if(vec_ditrk[ditrkidx1].first == vec_ditrk[ditrkidx2].first || vec_ditrk[ditrkidx1].first == vec_ditrk[ditrkidx2].second || vec_ditrk[ditrkidx1].second == vec_ditrk[ditrkidx2].first || vec_ditrk[ditrkidx1].second == vec_ditrk[ditrkidx2].second) continue;

	  if(std::get<4>(vec_ditrk_properties[ditrkidx2])==0 || std::get<6>(vec_ditrk_properties[ditrkidx2])==0) continue;
          double dca_idx2          = std::get<7>(vec_ditrk_properties[ditrkidx2]);
          if(dca_idx2 > diTrack2_dca_) continue;
          if (vec_trk_ttrk[vec_ditrk[ditrkidx2].first].first.pt() < DtkPtCut_) continue;
          if (vec_trk_ttrk[vec_ditrk[ditrkidx2].second].first.pt() < DtkPtCut_) continue;

	  // A quick Mass fit
          const auto& tTrk3 = vec_trk_ttrk[vec_ditrk[ditrkidx2].first].second;
          const auto& tTrk4 = vec_trk_ttrk[vec_ditrk[ditrkidx2].second].second;
	  auto mom1 = vec_trk_ttrk[vec_ditrk[ditrkidx2].first].first.momentum();
	  auto mom2 = vec_trk_ttrk[vec_ditrk[ditrkidx2].second].first.momentum();
          reco::Candidate::LorentzVector p4_1(mom1.x(), mom1.y(), mom1.z(), sqrt(mom1.mag2() + piMass * piMass));
          reco::Candidate::LorentzVector p4_2(mom2.x(), mom2.y(), mom2.z(), sqrt(mom2.mag2() + piMass * piMass));
          reco::Candidate::LorentzVector total_p4 = p4_1 + p4_2 + theKinVtxFitter_p4;
          float D0_premass = total_p4.mass();
          if(abs(D0_premass-D0Mass)>D0MassCut_*1.5) continue;
          if(total_p4.pt()<D0_PtCut_) continue;
          if(abs(total_p4.eta())>2.4) continue;
          if(verbose>=3) std::cout<< "Pass D0 presel" << std::endl;

	  // Kin fit for D0
	  KinVtxFitter D0_KinFitter(
                  {tTrk3, tTrk4, V0TT },
                  {piMass, piMass,  Kin_mass},
                  {PI_SIGMA, PI_SIGMA, Kin_massErr}
		  );
          //std::cout<<"Mass after fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
	  if ( !D0_KinFitter.success() ) continue;
          if ( D0_KinFitter.chi2()<0 || D0_KinFitter.dof()<0) continue;
          if ( D0_KinFitter.prob()<0.001) continue;
          auto D0_Kinfit_p4 = D0_KinFitter.fitted_p4();
          if(abs(D0_Kinfit_p4.mass()-D0Mass)>D0MassCut_) continue;    
	  if(D0_Kinfit_p4.pt() < D0_PtCut_) continue;
          auto D0_pv_lxy  =  l_xy(D0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
          if (D0_pv_lxy.error()!=0 && abs(D0_pv_lxy.value()/D0_pv_lxy.error()) < D0vtxDecaySigXYCut_ ) continue;

	  if(verbose>=2) std::cout<<" D0 looks good "<<std::endl;
	  if((cos_theta_3D(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4()))<cosThetaXYCut_) continue;

          auto DTT = D0_KinFitter.fitted_candidate_ttrk();
          float DKin_mass = D0_KinFitter.fitted_candidate().mass();
          float DKin_massErr = sqrt(D0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6));

          auto Ks0_d0_lxyz =  l_xyz(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z());
          double Ks0_Kin_d0_l_xyzSig = Ks0_d0_lxyz.error()==0? -99: abs(Ks0_d0_lxyz.value()/Ks0_d0_lxyz.error());
	  if(Ks0_Kin_d0_l_xyzSig < Ks0_l_xyzSigCut_) continue;

	  // Make a B
          for (unsigned int Btrack_idx = 0; Btrack_idx < vec_trk_ttrk.size(); ++Btrack_idx) {
              if(Btrack_idx == vec_ditrk[ditrkidx2].first || Btrack_idx == vec_ditrk[ditrkidx2].second || Btrack_idx == vec_ditrk[ditrkidx1].second || Btrack_idx == vec_ditrk[ditrkidx1].first) continue;
              const auto& Btrack = vec_trk_ttrk[Btrack_idx].first;
              if(Btrack.pt() < BtkPtCut_) continue;
              const auto& tBtrack = vec_trk_ttrk[Btrack_idx].second;
              auto mom3 = Btrack.momentum();
              std::pair<double, double> DCA_Btrack_pv = computeDCA(tBtrack, referencePos.x(), referencePos.y(), referencePos.z());
              double B_Kin_trk_pv_dcaSig = DCA_Btrack_pv.second==0? -99: abs(DCA_Btrack_pv.first/DCA_Btrack_pv.second);
              if( B_Kin_trk_pv_dcaSig < Btrk_dcaSigCut_) continue;

              // A quick Mass fit, assuming pi
              reco::Candidate::LorentzVector p4_3(mom3.x(), mom3.y(), mom3.z(), sqrt(mom3.mag2() + piMass * piMass));
              reco::Candidate::LorentzVector B_pre_p4 = total_p4 + p4_3;
              float B_premass = B_pre_p4.mass();
	      if (B_pre_p4.pt() < B_PtCut_) continue;
	      // A quick Mass fit, assuming kaon
	      reco::Candidate::LorentzVector p4kaon_3(mom3.x(), mom3.y(), mom3.z(), sqrt(mom3.mag2() + kplusMass * kplusMass));
              reco::Candidate::LorentzVector B_pre_p4kaon = total_p4 + p4kaon_3;
              float B_premass_kaon = B_pre_p4kaon.mass();
              
	      //if(abs(B_premass-BuMass)>BMassCut_*1.5 && abs(B_premass_kaon-BuMass)>BMassCut_*1.5) continue;
	      if(  B_premass < BuMass - 0.12 - BMassCut_*1.5  || B_premass > BuMass + BMassCut_*1.5 ) continue;
              if(verbose>=3) std::cout<< "Pass B presel" << std::endl;

	      // kin fit for B, assuming pi
	      KinVtxFitter B_KinFitter(
                      {tBtrack, DTT },
                      {piMass, DKin_mass},
                      {PI_SIGMA, DKin_massErr}
                      );
              // kin fit for B, assuming kaon
              KinVtxFitter B_KinFitter_kaon(
                      {tBtrack, DTT },
                      {kplusMass, DKin_mass},
                      {K_SIGMA, DKin_massErr}
                      );

              if ( !B_KinFitter.success() && !B_KinFitter_kaon.success() ) continue;
              if ( (B_KinFitter.chi2()<0 || B_KinFitter.dof()<0) && (B_KinFitter_kaon.chi2()<0 || B_KinFitter_kaon.dof()<0) ) continue;
	                    
	      auto B_Kinfit_p4 = B_KinFitter.fitted_p4();
              auto B_Kinfit_p4kaon = B_KinFitter_kaon.fitted_p4();
              if (B_Kinfit_p4.pt() < B_PtCut_) continue;

	      if(verbose>=3) std::cout<< "B prob, dmass, pt, eta: "<<B_KinFitter.prob()<<" "<<abs(B_Kinfit_p4.mass()-BuMass)<<" "<<B_Kinfit_p4.pt()<<" "<<abs(B_Kinfit_p4.eta())<<" "<<cos_theta_3D(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z(),  D0_KinFitter.fitted_p4())<<" "<<cos_theta_3D(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter.fitted_p4()) << std::endl;

	      if ( B_KinFitter.prob()<=0 && B_KinFitter_kaon.prob()<=0) continue;
	      //if(abs(B_Kinfit_p4.mass()-BuMass)>BMassCut_ && abs(B_Kinfit_p4kaon.mass()-BuMass)>BMassCut_) continue;
              if(  B_premass < BuMass - 0.12 - BMassCut_  || B_premass > BuMass + BMassCut_ ) continue;
	      if(B_Kinfit_p4.pt()<1 && B_Kinfit_p4kaon.pt()<1) continue;
              if(abs(B_Kinfit_p4.eta())>2.4 && abs(B_Kinfit_p4kaon.eta())>2.4) continue;

	      //if((cos_theta_3D(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z(),  D0_KinFitter.fitted_p4()))<cosThetaXYCut_) continue;
	      //if((cos_theta_3D(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter.fitted_p4()))<cosThetaXYCut_) continue;
              if(verbose>=2) std::cout<< "B looks good" << std::endl;

              // Save something
              pat::CompositeCandidate Bu;
	      // idx
	      Bu.addUserInt("DiTrack_idx1", int(ditrkidx1));
              Bu.addUserInt("DiTrack_idx2", int(ditrkidx2));
              Bu.addUserInt("Track_idx1",   int(vec_ditrk[ditrkidx1].first));
              Bu.addUserInt("Track_idx2",   int(vec_ditrk[ditrkidx1].second));
              Bu.addUserInt("Track_idx3",   int(vec_ditrk[ditrkidx2].first));
              Bu.addUserInt("Track_idx4",   int(vec_ditrk[ditrkidx2].second));
	      // Ditrack variables
	      // dca relative to BS or PV
              Bu.addUserFloat("DiTrk1_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx1]));		      
	      Bu.addUserFloat("DiTrk1_dot",         std::get<2>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_dca",         abs(std::get<7>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_vtx_r",   sqrt(std::get<9>(vec_ditrk_properties[ditrkidx1])*std::get<9>(vec_ditrk_properties[ditrkidx1])+std::get<10>(vec_ditrk_properties[ditrkidx1])*std::get<10>(vec_ditrk_properties[ditrkidx1]))   );
              Bu.addUserFloat("DiTrk1_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx1]));
	      // dca relative to bs and pv
              Bu.addUserFloat("DiTrk1_trk1_bs_dca", abs(std::get<3>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_trk2_bs_dca", abs(std::get<5>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_trk1_pv_dca", abs(std::get<18>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_trk2_pv_dca", abs(std::get<20>(vec_ditrk_properties[ditrkidx1])));
	      Bu.addUserFloat("DiTrk1_trk1_bs_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrkidx1]))/std::get<4>(vec_ditrk_properties[ditrkidx1]));
	      Bu.addUserFloat("DiTrk1_trk2_bs_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrkidx1]))/std::get<6>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_trk1_pv_dcaSig", abs(std::get<18>(vec_ditrk_properties[ditrkidx1]))/std::get<19>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_trk2_pv_dcaSig", abs(std::get<20>(vec_ditrk_properties[ditrkidx1]))/std::get<21>(vec_ditrk_properties[ditrkidx1]));
	      // displacement relative to bs and pv
              Bu.addUserFloat("DiTrk1_KLM_bs_lxy",            abs(std::get<15>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_KLM_bs_lxyErr",         std::get<16>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_bs_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_pv_lxy",            abs(std::get<22>(vec_ditrk_properties[ditrkidx1])));
              Bu.addUserFloat("DiTrk1_KLM_pv_lxyErr",         std::get<23>(vec_ditrk_properties[ditrkidx1]));
              Bu.addUserFloat("DiTrk1_KLM_pv_cos_theta_XY",   std::get<24>(vec_ditrk_properties[ditrkidx1]));

              Bu.addUserFloat("DiTrk2_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_cxPtz",       std::get<1>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dot",         std::get<2>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_dca",         abs(std::get<7>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_massSquared", std::get<8>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_vtx_r",          sqrt(std::get<9>(vec_ditrk_properties[ditrkidx2])*std::get<9>(vec_ditrk_properties[ditrkidx2])+std::get<10>(vec_ditrk_properties[ditrkidx2])*std::get<10>(vec_ditrk_properties[ditrkidx2]))   );
              Bu.addUserFloat("DiTrk2_KLM_vtx_z",          std::get<11>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_chi2",           std::get<12>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_ndof",           std::get<13>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_normalizedChi2", std::get<14>(vec_ditrk_properties[ditrkidx2]));
              // dca relative to bs and pv
              Bu.addUserFloat("DiTrk2_trk1_bs_dca", abs(std::get<3>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_trk2_bs_dca", abs(std::get<5>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_trk1_pv_dca", abs(std::get<18>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_trk2_pv_dca", abs(std::get<20>(vec_ditrk_properties[ditrkidx2])));
	      Bu.addUserFloat("DiTrk2_trk1_bs_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrkidx2]))/std::get<4>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk2_bs_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrkidx2]))/std::get<6>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk1_pv_dcaSig", abs(std::get<18>(vec_ditrk_properties[ditrkidx2]))/std::get<19>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_trk2_pv_dcaSig", abs(std::get<20>(vec_ditrk_properties[ditrkidx2]))/std::get<21>(vec_ditrk_properties[ditrkidx2]));
              // displacement relative to bs and pv
              Bu.addUserFloat("DiTrk2_KLM_bs_lxy",            abs(std::get<15>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_KLM_bs_lxyErr",         std::get<16>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_bs_cos_theta_XY",   std::get<17>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_pv_lxy",            abs(std::get<22>(vec_ditrk_properties[ditrkidx2])));
              Bu.addUserFloat("DiTrk2_KLM_pv_lxyErr",         std::get<23>(vec_ditrk_properties[ditrkidx2]));
              Bu.addUserFloat("DiTrk2_KLM_pv_cos_theta_XY",   std::get<24>(vec_ditrk_properties[ditrkidx2]));

	      // Ks0
	      // Basic
              Bu.addUserFloat("Ks0_Kin_vtx_x", theKinVtxFitter.fitted_vtx().x());
              Bu.addUserFloat("Ks0_Kin_vtx_y", theKinVtxFitter.fitted_vtx().y());
              Bu.addUserFloat("Ks0_Kin_vtx_r", sqrt( theKinVtxFitter.fitted_vtx().x()*theKinVtxFitter.fitted_vtx().x() + theKinVtxFitter.fitted_vtx().y()*theKinVtxFitter.fitted_vtx().y() ) );
	      Bu.addUserFloat("Ks0_Kin_vtx_z", theKinVtxFitter.fitted_vtx().z());
              Bu.addUserFloat("Ks0_Kin_chi2",  theKinVtxFitter.chi2());
              Bu.addUserFloat("Ks0_Kin_dof",   theKinVtxFitter.dof()); // always 1 ?
              Bu.addUserFloat("Ks0_Kin_prob",  theKinVtxFitter.prob());
              Bu.addUserFloat("Ks0_Kin_pt",    theKinVtxFitter.fitted_p4().pt());
              Bu.addUserFloat("Ks0_Kin_eta",   theKinVtxFitter.fitted_p4().eta());
              Bu.addUserFloat("Ks0_Kin_phi",   theKinVtxFitter.fitted_p4().phi());
              Bu.addUserFloat("Ks0_Kin_mass",  theKinVtxFitter.fitted_candidate().mass());
              Bu.addUserFloat("Ks0_Kin_massErr", sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("Ks0_Kin_trk1_pt", theKinVtxFitter.daughter_p4(0).pt());
              Bu.addUserFloat("Ks0_Kin_trk1_eta", theKinVtxFitter.daughter_p4(0).eta());
              Bu.addUserFloat("Ks0_Kin_trk1_phi", theKinVtxFitter.daughter_p4(0).phi());
              Bu.addUserFloat("Ks0_Kin_trk2_pt", theKinVtxFitter.daughter_p4(1).pt());
              Bu.addUserFloat("Ks0_Kin_trk2_eta", theKinVtxFitter.daughter_p4(1).eta());
              Bu.addUserFloat("Ks0_Kin_trk2_phi", theKinVtxFitter.daughter_p4(1).phi());
	      // pointing angle relative to bs, pv, D0
              Bu.addUserFloat("Ks0_Kin_bs_alpha_2D", (cos_theta_2D(theKinVtxFitter, *theBeamSpot, theKinVtxFitter.fitted_p4())));
              Bu.addUserFloat("Ks0_Kin_pv_alpha_2D", (cos_theta_2D(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z(), theKinVtxFitter.fitted_p4())));
              Bu.addUserFloat("Ks0_Kin_pv_alpha_3D", (cos_theta_3D(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z(), theKinVtxFitter.fitted_p4())));
              Bu.addUserFloat("Ks0_Kin_d0_alpha_2D", (cos_theta_2D(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4())));
              Bu.addUserFloat("Ks0_Kin_d0_alpha_3D", (cos_theta_3D(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4())));
	      // displacement relative to bs, pv, D0
              // dca relative to bs and pv
              auto Ks0_bs_lxy  =  l_xy(theKinVtxFitter, *theBeamSpot);
              auto Ks0_pv_lxy  =  l_xy(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto Ks0_pv_lxyz =  l_xyz(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto Ks0_d0_lxy  =  l_xy(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z());
              auto Ks0_d0_lxyz =  l_xyz(theKinVtxFitter, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z());
              double Ks0_Kin_bs_l_xySig  = Ks0_bs_lxy.error()==0?  -99: abs(Ks0_bs_lxy.value()/Ks0_bs_lxy.error());
              double Ks0_Kin_pv_l_xySig  = Ks0_pv_lxy.error()==0?  -99: abs(Ks0_pv_lxy.value()/Ks0_pv_lxy.error());
              double Ks0_Kin_pv_l_xyzSig = Ks0_pv_lxyz.error()==0? -99: abs(Ks0_pv_lxyz.value()/Ks0_pv_lxyz.error());
              double Ks0_Kin_d0_l_xySig  = Ks0_d0_lxy.error()==0?  -99: abs(Ks0_d0_lxy.value()/Ks0_d0_lxy.error());
              double Ks0_Kin_d0_l_xyzSig = Ks0_d0_lxyz.error()==0? -99: abs(Ks0_d0_lxyz.value()/Ks0_d0_lxyz.error());
              Bu.addUserFloat("Ks0_Kin_bs_l_xy",     Ks0_bs_lxy.value());
              Bu.addUserFloat("Ks0_Kin_bs_l_xySig",  Ks0_Kin_bs_l_xySig);
              Bu.addUserFloat("Ks0_Kin_pv_l_xy",     Ks0_pv_lxy.value());
              Bu.addUserFloat("Ks0_Kin_pv_l_xySig",  Ks0_Kin_pv_l_xySig);
              Bu.addUserFloat("Ks0_Kin_pv_l_xyz",    Ks0_pv_lxyz.value());
              Bu.addUserFloat("Ks0_Kin_pv_l_xyzSig", Ks0_Kin_pv_l_xyzSig);
              Bu.addUserFloat("Ks0_Kin_d0_l_xy",     Ks0_d0_lxy.value());
              Bu.addUserFloat("Ks0_Kin_d0_l_xySig",  Ks0_Kin_d0_l_xySig);
              Bu.addUserFloat("Ks0_Kin_d0_l_xyz",    Ks0_d0_lxyz.value());
              Bu.addUserFloat("Ks0_Kin_d0_l_xyzSig", Ks0_Kin_d0_l_xyzSig);
              // dca relative to bs and pv
              std::pair<double, double> DCA_ks0 = computeDCA(V0TT, D0_KinFitter.fitted_vtx().x(), D0_KinFitter.fitted_vtx().y(), D0_KinFitter.fitted_vtx().z());
              double Ks0_Kin_b_dcaSig = DCA_ks0.second==0? -99: abs(DCA_ks0.first/DCA_ks0.second);
              Bu.addUserFloat("Ks0_Kin_d0_dca", abs(DCA_ks0.first));
              Bu.addUserFloat("Ks0_Kin_d0_dcaSig", Ks0_Kin_b_dcaSig);


	      // D0
              // Basic
              Bu.addUserFloat("D0_premass", D0_premass);
              Bu.addUserFloat("D0_Kin_vtx_x", D0_KinFitter.fitted_vtx().x());
              Bu.addUserFloat("D0_Kin_vtx_y", D0_KinFitter.fitted_vtx().y());
              Bu.addUserFloat("D0_Kin_vtx_r", sqrt( D0_KinFitter.fitted_vtx().x()*D0_KinFitter.fitted_vtx().x() + D0_KinFitter.fitted_vtx().y()*D0_KinFitter.fitted_vtx().y() ) );
	      Bu.addUserFloat("D0_Kin_vtx_z", D0_KinFitter.fitted_vtx().z());
              Bu.addUserFloat("D0_Kin_chi2", D0_KinFitter.chi2());
              Bu.addUserFloat("D0_Kin_dof", D0_KinFitter.dof());
              Bu.addUserFloat("D0_Kin_prob", D0_KinFitter.prob());
              Bu.addUserFloat("D0_Kin_pt", D0_Kinfit_p4.pt() );
              Bu.addUserFloat("D0_Kin_eta", D0_Kinfit_p4.eta() );
              Bu.addUserFloat("D0_Kin_phi", D0_Kinfit_p4.phi() );
              Bu.addUserFloat("D0_Kin_mass", D0_Kinfit_p4.mass() );
              Bu.addUserFloat("D0_Kin_massErr", sqrt(D0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("D0_Kin_trk3_pt", D0_KinFitter.daughter_p4(0).pt());
              Bu.addUserFloat("D0_Kin_trk3_eta", D0_KinFitter.daughter_p4(0).eta());
              Bu.addUserFloat("D0_Kin_trk3_phi", D0_KinFitter.daughter_p4(0).phi());
              Bu.addUserFloat("D0_Kin_trk4_pt", D0_KinFitter.daughter_p4(1).pt());
              Bu.addUserFloat("D0_Kin_trk4_eta", D0_KinFitter.daughter_p4(1).eta());
              Bu.addUserFloat("D0_Kin_trk4_phi", D0_KinFitter.daughter_p4(1).phi());
              Bu.addUserFloat("D0_Kin_ks0_pt", D0_KinFitter.daughter_p4(2).pt());
              Bu.addUserFloat("D0_Kin_ks0_eta", D0_KinFitter.daughter_p4(2).eta());
              Bu.addUserFloat("D0_Kin_ks0_phi", D0_KinFitter.daughter_p4(2).phi());
              // pointing angle relative to bs, pv, D0
              Bu.addUserFloat("D0_Kin_bs_alpha_2D", (cos_theta_2D(D0_KinFitter, *theBeamSpot, D0_KinFitter.fitted_p4())));
              Bu.addUserFloat("D0_Kin_pv_alpha_2D", (cos_theta_2D(D0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), D0_KinFitter.fitted_p4())));
              Bu.addUserFloat("D0_Kin_pv_alpha_3D", (cos_theta_3D(D0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), D0_KinFitter.fitted_p4())));
              Bu.addUserFloat("D0_Kin_b_alpha_2D",  (cos_theta_2D(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z(),  D0_KinFitter.fitted_p4())));
              Bu.addUserFloat("D0_Kin_b_alpha_3D",  (cos_theta_3D(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z(),  D0_KinFitter.fitted_p4())));
              // displacement relative to bs, pv, D0
              auto D0_bs_lxy  =  l_xy(D0_KinFitter, *theBeamSpot);
              auto D0_pv_lxy  =  l_xy(D0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto D0_pv_lxyz =  l_xyz(D0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto D0_b_lxy  =  l_xy(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z());
              auto D0_b_lxyz =  l_xyz(D0_KinFitter, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z());
              double D0_Kin_bs_l_xySig  = D0_bs_lxy.error()==0?  -99: abs(D0_bs_lxy.value()/D0_bs_lxy.error());
              double D0_Kin_pv_l_xySig  = D0_pv_lxy.error()==0?  -99: abs(D0_pv_lxy.value()/D0_pv_lxy.error());
              double D0_Kin_pv_l_xyzSig = D0_pv_lxyz.error()==0? -99: abs(D0_pv_lxyz.value()/D0_pv_lxyz.error());
              double D0_Kin_b_l_xySig   = D0_b_lxy.error()==0?   -99: abs(D0_b_lxy.value()/D0_b_lxy.error());
              double D0_Kin_b_l_xyzSig  = D0_b_lxyz.error()==0?  -99: abs(D0_b_lxyz.value()/D0_b_lxyz.error());
              Bu.addUserFloat("D0_Kin_bs_l_xy",     D0_bs_lxy.value());
              Bu.addUserFloat("D0_Kin_bs_l_xySig",  D0_Kin_bs_l_xySig);
              Bu.addUserFloat("D0_Kin_pv_l_xy",     D0_pv_lxy.value());
              Bu.addUserFloat("D0_Kin_pv_l_xySig",  D0_Kin_pv_l_xySig);
              Bu.addUserFloat("D0_Kin_pv_l_xyz",    D0_pv_lxyz.value());
              Bu.addUserFloat("D0_Kin_pv_l_xyzSig", D0_Kin_pv_l_xyzSig);
              Bu.addUserFloat("D0_Kin_b_l_xy",      D0_b_lxy.value());
              Bu.addUserFloat("D0_Kin_b_l_xySig",   D0_Kin_b_l_xySig);
              Bu.addUserFloat("D0_Kin_b_l_xyz",     D0_b_lxyz.value());
              Bu.addUserFloat("D0_Kin_b_l_xyzSig",  D0_Kin_b_l_xyzSig);
              // dca relative to bs and pv
              std::pair<double, double> DCA_d0 = computeDCA(DTT, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z());
              double D0_Kin_b_dcaSig = DCA_d0.second==0? -99: abs(DCA_d0.first/DCA_d0.second);
	      Bu.addUserFloat("D0_Kin_b_dca", abs(DCA_d0.first));
              Bu.addUserFloat("D0_Kin_b_dcaSig", D0_Kin_b_dcaSig);

	      // B
              // D0pi
    	      Bu.addUserInt("BTrack_idx", int(Btrack_idx));
	      Bu.addUserFloat("B_premass", B_premass);
              Bu.addUserFloat("B_Kin_vtx_x", B_KinFitter.fitted_vtx().x());
              Bu.addUserFloat("B_Kin_vtx_y", B_KinFitter.fitted_vtx().y());
              Bu.addUserFloat("B_Kin_vtx_r", sqrt(B_KinFitter.fitted_vtx().x()*B_KinFitter.fitted_vtx().x() + B_KinFitter.fitted_vtx().y()*B_KinFitter.fitted_vtx().y() ));
	      Bu.addUserFloat("B_Kin_vtx_z", B_KinFitter.fitted_vtx().z());
              Bu.addUserFloat("B_Kin_chi2", B_KinFitter.chi2());
              Bu.addUserFloat("B_Kin_dof", B_KinFitter.dof());
              Bu.addUserFloat("B_Kin_prob", B_KinFitter.prob());
              Bu.addUserFloat("B_Kin_pt", B_Kinfit_p4.pt() );
              Bu.addUserFloat("B_Kin_eta", B_Kinfit_p4.eta() );
              Bu.addUserFloat("B_Kin_phi", B_Kinfit_p4.phi() );
              Bu.addUserFloat("B_Kin_mass", B_Kinfit_p4.mass() );
              Bu.addUserFloat("B_Kin_massErr", sqrt(B_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("B_Kin_trk_charge", Btrack.charge());
              Bu.addUserFloat("B_Kin_trk_pt", B_KinFitter.daughter_p4(0).pt());
              Bu.addUserFloat("B_Kin_trk_eta", B_KinFitter.daughter_p4(0).eta());
              Bu.addUserFloat("B_Kin_trk_phi", B_KinFitter.daughter_p4(0).phi());
              Bu.addUserFloat("B_Kin_D0_pt", B_KinFitter.daughter_p4(1).pt());
              Bu.addUserFloat("B_Kin_D0_eta", B_KinFitter.daughter_p4(1).eta());
              Bu.addUserFloat("B_Kin_D0_phi", B_KinFitter.daughter_p4(1).phi());

	      // pointing angle relative to bs, pv, D0
              Bu.addUserFloat("B_Kin_bs_alpha_2D", (cos_theta_2D(B_KinFitter, *theBeamSpot, B_KinFitter.fitted_p4())));
              Bu.addUserFloat("B_Kin_pv_alpha_2D", (cos_theta_2D(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter.fitted_p4())));
              Bu.addUserFloat("B_Kin_pv_alpha_3D", (cos_theta_3D(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter.fitted_p4())));
             // displacement relative to bs, pv, D0
              auto B_bs_lxy  =  l_xy(B_KinFitter, *theBeamSpot);
              auto B_pv_lxy  =  l_xy(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto B_pv_lxyz =  l_xyz(B_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              double B_Kin_bs_l_xySig  = B_bs_lxy.error()==0?  -99: abs(B_bs_lxy.value()/B_bs_lxy.error());
              double B_Kin_pv_l_xySig  = B_pv_lxy.error()==0?  -99: abs(B_pv_lxy.value()/B_pv_lxy.error());
              double B_Kin_pv_l_xyzSig = B_pv_lxyz.error()==0? -99: abs(B_pv_lxyz.value()/B_pv_lxyz.error());
              Bu.addUserFloat("B_Kin_bs_l_xy",     B_bs_lxy.value());
              Bu.addUserFloat("B_Kin_bs_l_xySig",  B_Kin_bs_l_xySig);
              Bu.addUserFloat("B_Kin_pv_l_xy",     B_pv_lxy.value());
              Bu.addUserFloat("B_Kin_pv_l_xySig",  B_Kin_pv_l_xySig);
              Bu.addUserFloat("B_Kin_pv_l_xyz",    B_pv_lxyz.value());
              Bu.addUserFloat("B_Kin_pv_l_xyzSig", B_Kin_pv_l_xyzSig);
              std::pair<double, double> DCA_Btrack_beamspot = computeDCA(tBtrack, *theBeamSpot);
	      double B_Kin_trk_bs_dcaSig = DCA_Btrack_beamspot.second==0? -99: abs(DCA_Btrack_beamspot.first/DCA_Btrack_beamspot.second);
	      Bu.addUserFloat("B_Kin_trk_bs_dca", abs(DCA_Btrack_beamspot.first));
	      Bu.addUserFloat("B_Kin_trk_bs_dcaSig", B_Kin_trk_bs_dcaSig);
              Bu.addUserFloat("B_Kin_trk_pv_dca", abs(DCA_Btrack_pv.first));
	      Bu.addUserFloat("B_Kin_trk_pv_dcaSig", B_Kin_trk_pv_dcaSig);
              std::pair<double, double> DCA_Btrack_b = computeDCA(tBtrack, B_KinFitter.fitted_vtx().x(), B_KinFitter.fitted_vtx().y(), B_KinFitter.fitted_vtx().z());
              double B_Kin_trk_b_dcaSig = DCA_Btrack_b.second==0? -99: abs(DCA_Btrack_b.first/DCA_Btrack_b.second);
	      Bu.addUserFloat("B_Kin_trk_b_dca", abs(DCA_Btrack_b.first));
              Bu.addUserFloat("B_Kin_trk_b_dcaSig", B_Kin_trk_b_dcaSig);
	      Bu.setCharge(Btrack.charge());
	      Bu.setP4(B_Kinfit_p4);
	      //std::cout<<"B_Kin_trk_b_dcaSig "<<DCA_Btrack_b.first<<" "<<DCA_Btrack_b.second<<" "<< Bu.userFloat("B_Kin_trk_b_dcaSig")<<std::endl;

              //if(abs(B_Kinfit_p4.mass()-BuMass) < abs(B_Kinfit_p4kaon.mass()-BuMass)){
	      //        Bu.setP4(B_Kinfit_p4);
	      //}else{
	      //        Bu.setP4(B_Kinfit_p4kaon);
	      //}

              // D0K
              Bu.addUserFloat("B_kaon_premass", B_premass_kaon);
              Bu.addUserFloat("B_kaon_Kin_vtx_x", B_KinFitter_kaon.fitted_vtx().x());
              Bu.addUserFloat("B_kaon_Kin_vtx_y", B_KinFitter_kaon.fitted_vtx().y());
              Bu.addUserFloat("B_kaon_Kin_vtx_r", sqrt(B_KinFitter_kaon.fitted_vtx().x()*B_KinFitter_kaon.fitted_vtx().x() + B_KinFitter_kaon.fitted_vtx().y()*B_KinFitter_kaon.fitted_vtx().y() ));
              Bu.addUserFloat("B_kaon_Kin_vtx_z", B_KinFitter_kaon.fitted_vtx().z());
              Bu.addUserFloat("B_kaon_Kin_chi2", B_KinFitter_kaon.chi2());
              Bu.addUserFloat("B_kaon_Kin_dof", B_KinFitter_kaon.dof());
              Bu.addUserFloat("B_kaon_Kin_prob", B_KinFitter_kaon.prob());
              Bu.addUserFloat("B_kaon_Kin_pt", B_Kinfit_p4kaon.pt() );
              Bu.addUserFloat("B_kaon_Kin_eta", B_Kinfit_p4kaon.eta() );
              Bu.addUserFloat("B_kaon_Kin_phi", B_Kinfit_p4kaon.phi() );
              Bu.addUserFloat("B_kaon_Kin_mass", B_Kinfit_p4kaon.mass() );
              Bu.addUserFloat("B_kaon_Kin_massErr", sqrt(B_KinFitter_kaon.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Bu.addUserFloat("B_kaon_Kin_trk_charge", Btrack.charge());
              Bu.addUserFloat("B_kaon_Kin_trk_pt", B_KinFitter_kaon.daughter_p4(0).pt());
              Bu.addUserFloat("B_kaon_Kin_trk_eta", B_KinFitter_kaon.daughter_p4(0).eta());
              Bu.addUserFloat("B_kaon_Kin_trk_phi", B_KinFitter_kaon.daughter_p4(0).phi());
              Bu.addUserFloat("B_kaon_Kin_D0_pt", B_KinFitter_kaon.daughter_p4(1).pt());
              Bu.addUserFloat("B_kaon_Kin_D0_eta", B_KinFitter_kaon.daughter_p4(1).eta());
              Bu.addUserFloat("B_kaon_Kin_D0_phi", B_KinFitter_kaon.daughter_p4(1).phi());

              // pointing angle relative to bs, pv, D0
              Bu.addUserFloat("B_kaon_Kin_bs_alpha_2D", (cos_theta_2D(B_KinFitter_kaon, *theBeamSpot, B_KinFitter_kaon.fitted_p4())));
              Bu.addUserFloat("B_kaon_Kin_pv_alpha_2D", (cos_theta_2D(B_KinFitter_kaon, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter_kaon.fitted_p4())));
              Bu.addUserFloat("B_kaon_Kin_pv_alpha_3D", (cos_theta_3D(B_KinFitter_kaon, referencePos.x(), referencePos.y(), referencePos.z(), B_KinFitter_kaon.fitted_p4())));
              // displacement relative to bs, pv, D0
              auto   B_kaon_bs_lxy  =  l_xy(B_KinFitter_kaon, *theBeamSpot);
              auto   B_kaon_pv_lxy  =  l_xy(B_KinFitter_kaon, referencePos.x(), referencePos.y(), referencePos.z());
              auto   B_kaon_pv_lxyz =  l_xyz(B_KinFitter_kaon, referencePos.x(), referencePos.y(), referencePos.z());
              double B_kaon_Kin_bs_l_xySig  = B_kaon_bs_lxy.error()==0?  -99: abs(B_kaon_bs_lxy.value() /B_kaon_bs_lxy.error());
              double B_kaon_Kin_pv_l_xySig  = B_kaon_pv_lxy.error()==0?  -99: abs(B_kaon_pv_lxy.value() /B_kaon_pv_lxy.error());
              double B_kaon_Kin_pv_l_xyzSig = B_kaon_pv_lxyz.error()==0? -99: abs(B_kaon_pv_lxyz.value()/B_kaon_pv_lxyz.error());
              Bu.addUserFloat("B_kaon_Kin_bs_l_xy",     B_kaon_bs_lxy.value());
              Bu.addUserFloat("B_kaon_Kin_bs_l_xySig",  B_kaon_Kin_bs_l_xySig);
              Bu.addUserFloat("B_kaon_Kin_pv_l_xy",     B_kaon_pv_lxy.value());
              Bu.addUserFloat("B_kaon_Kin_pv_l_xySig",  B_kaon_Kin_pv_l_xySig);
              Bu.addUserFloat("B_kaon_Kin_pv_l_xyz",    B_kaon_pv_lxyz.value());
              Bu.addUserFloat("B_kaon_Kin_pv_l_xyzSig", B_kaon_Kin_pv_l_xyzSig);
              //for (unsigned int i=0; i < BMva_.size(); ++i) {
              //  BMva_.at(i).set("B_DiTrk1_dca",            Bu.userFloat("DiTrk1_dca"));
              //  BMva_.at(i).set("B_Ks0_Kin_vtx_r",         Bu.userFloat("Ks0_Kin_vtx_r"));
              //  BMva_.at(i).set("B_Ks0_Kin_d0_alpha_2D",   Bu.userFloat("Ks0_Kin_d0_alpha_2D"));
              //  BMva_.at(i).set("B_Ks0_Kin_d0_l_xy",       Bu.userFloat("Ks0_Kin_d0_l_xy"));
              //  BMva_.at(i).set("B_Ks0_Kin_d0_l_xySig",    Bu.userFloat("Ks0_Kin_d0_l_xySig"));
              //  BMva_.at(i).set("B_Ks0_Kin_d0_dca",        Bu.userFloat("Ks0_Kin_d0_dca"));
              //  BMva_.at(i).set("B_Ks0_Kin_d0_dcaSig",     Bu.userFloat("Ks0_Kin_d0_dcaSig"));
              //  BMva_.at(i).set("B_D0_Kin_vtx_r",          Bu.userFloat("D0_Kin_vtx_r"));
              //  BMva_.at(i).set("B_D0_Kin_prob",           Bu.userFloat("D0_Kin_prob"));
              //  BMva_.at(i).set("B_D0_Kin_b_alpha_2D",     Bu.userFloat("D0_Kin_b_alpha_2D"));
              //  BMva_.at(i).set("B_D0_Kin_b_l_xy",         Bu.userFloat("D0_Kin_b_l_xy"));
              //  BMva_.at(i).set("B_D0_Kin_b_l_xyz",        Bu.userFloat("D0_Kin_b_l_xyz"));
              //  BMva_.at(i).set("B_D0_Kin_b_l_xySig",      Bu.userFloat("D0_Kin_b_l_xySig"));
              //  BMva_.at(i).set("B_D0_Kin_b_l_xyzSig",     Bu.userFloat("D0_Kin_b_l_xyzSig"));
              //  BMva_.at(i).set("B_D0_Kin_b_dcaSig",       Bu.userFloat("D0_Kin_b_dcaSig"));
              //  BMva_.at(i).set("B_B_Kin_vtx_r",           Bu.userFloat("B_Kin_vtx_r"));
              //  BMva_.at(i).set("B_B_Kin_prob",            Bu.userFloat("B_Kin_prob"));
              //  BMva_.at(i).set("B_B_Kin_pv_alpha_2D",     Bu.userFloat("B_Kin_pv_alpha_2D"));
              //  BMva_.at(i).set("B_B_Kin_pv_l_xy",         Bu.userFloat("B_Kin_pv_l_xy"));
              //  BMva_.at(i).set("B_B_Kin_pv_l_xySig",      Bu.userFloat("B_Kin_pv_l_xySig"));
              //  BMva_.at(i).set("B_B_Kin_trk_b_dca",       Bu.userFloat("B_Kin_trk_b_dca"));
              //  BMva_.at(i).set("B_B_Kin_trk_b_dcaSig",    Bu.userFloat("B_Kin_trk_b_dcaSig"));
              //  BMva_.at(i).set("B_B_Kin_pt",              Bu.userFloat("B_Kin_pt"));
              //  BMva_.at(i).set("B_B_Kin_D0_pt",           Bu.userFloat("B_Kin_D0_pt"));
              //  BMva_.at(i).set("B_B_Kin_trk_pt",          Bu.userFloat("B_Kin_trk_pt"));
              //  Bu.addUserFloat("xgb_" + xgboost_variable_names_.at(i), BMva_.at(i).predict());
              //}
	                    
	      if(std::isnan(abs(Bu.userFloat("Ks0_Kin_d0_l_xySig")))) continue;
	      if(std::isnan(abs(Bu.userFloat("Ks0_Kin_d0_dcaSig")))) continue;
	      if(std::isnan(abs(Bu.userFloat("D0_Kin_b_l_xySig")))) continue;
	      if(std::isnan(abs(Bu.userFloat("D0_Kin_b_l_xyzSig")))) continue;
              if(std::isnan(abs(Bu.userFloat("D0_Kin_b_dcaSig")))) continue;
              if(std::isnan(abs(Bu.userFloat("B_Kin_pv_l_xySig")))) continue;
              if(std::isnan(abs(Bu.userFloat("B_Kin_trk_b_dcaSig")))) continue;
	      for (unsigned int i=0; i < BMva_.size(); ++i) {
                  Bu.addUserFloat("xgb_" + xgboost_variable_names_.at(i), Bmva_estimator(Bu, i));
	      }
	      if(Bu.userFloat("xgb_" + xgboost_variable_names_.at(0))<0.2) continue;
	      Bu_out->push_back(Bu);
	  }
  
      }
  }

  if(savetrack_) iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");    
  if(savetrack_) iEvent.put(std::move(ditracks_earlyout),       "SelectedDiTracks");
//  iEvent.put(std::move(LooseKshorts_out), "LooseKshort");
//  iEvent.put(std::move(ret_val), "SelectedV0Collection");
//  iEvent.put(std::move(trans_out), "SelectedV0TransientCollection");
  iEvent.put(std::move(Bu_out), "B");
  return;

}

DEFINE_FWK_MODULE(BDhFitter_v2);

