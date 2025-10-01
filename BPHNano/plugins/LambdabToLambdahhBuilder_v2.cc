/// original authors: Yihui Lai (Princeton)
// lambdab0 -> lambda0 + hh


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
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include <tuple>  // For std::tuple

#include "PhysicsTools/BPHNano/plugins/XGBooster.h"

#include <iostream>
#include <chrono>

// pdg mass constants
namespace {
  const double protonMass = 0.938;
  
  const double piMass = 0.13957018;
  const float  pion_sigma = piMass*1.e-6;
  const double piMassSquared = piMass * piMass;

  const double ld0_Mass = 1.115683;
  const double ld0_MassSquared = ld0_Mass * ld0_Mass;
  const double kplusMass = 0.493677;
  const double Ldb0Mass = 5.61957;
}  // namespace

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

class LambdabToLambdahhBuilder_v2 : public edm::stream::EDProducer<> {

public:
  explicit LambdabToLambdahhBuilder_v2(const edm::ParameterSet &theParameters):
      bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      xgboost_models_( theParameters.getParameter<std::vector<std::string>>( "xgboost_models" ) ),
      xgboost_variable_names_( theParameters.getParameter<std::vector<std::string>>( "xgboost_variable_names" ) ),
      token_beamSpot(consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"))),
      token_vertices(consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"))),
      tracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("tracks"))),
      lostTracksToken_(consumes<edm::View<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("lostTracks"))),
      tkNHitsCut_( theParameters.getParameter<int>("tkNHitsCut")),
      tkPtCut_( theParameters.getParameter<double>("tkPtCut")),
      Lambdab_tkPtCut_( theParameters.getParameter<double>("Lambdab_tkPtCut")),
      tkEtaCut_( theParameters.getParameter<double>("tkEtaCut")),
      tkChi2Cut_( theParameters.getParameter<double>("tkChi2Cut")),
      tkIPSigXYCut_( theParameters.getParameter<double>("tkIPSigXYCut")),
      vtxChi2Cut_( theParameters.getParameter<double>("vtxChi2Cut")),
      vtxDecaySigXYCut_( theParameters.getParameter<double>("vtxDecaySigXYCut")),
      TrkSigXYCut_( theParameters.getParameter<double>("TrkSigXYCut")),
      vtxDecaySigXYZCut_( theParameters.getParameter<double>("vtxDecaySigXYZCut")),
      cosThetaXYCut_( theParameters.getParameter<double>("cosThetaXYCut")),
      cosThetaXYZCut_( theParameters.getParameter<double>("cosThetaXYZCut")),
      diTrack2_dca_( theParameters.getParameter<double>("diTrack2_dca")),
      Trk34SigXYCut_( theParameters.getParameter<double>("Trk34SigXYCut")),
      ld0_l_xyzSigCut_( theParameters.getParameter<double>("ld0_l_xyzSigCut")),
      ldb0_PtCut_( theParameters.getParameter<double>("ldb0_PtCut")),
      ldb0_vtxDecaySigXYCut_( theParameters.getParameter<double>("ldb0_vtxDecaySigXYCut")),
      ld0_MassCut_( theParameters.getParameter<double>("ld0_MassCut")),
      ldb0_MassCut_( theParameters.getParameter<double>("ldb0_MassCut")),
      savetrack_( theParameters.getParameter<bool>("savetrack")),
      verbose( theParameters.getParameter<int>("verbose"))
	{
      	  produces<pat::CompositeCandidateCollection>("SelectedTracks");
      	  produces<pat::CompositeCandidateCollection>("SelectedDiTracks");
	  produces<pat::CompositeCandidateCollection>("Ldb0");
          for (auto model: xgboost_models_)
              BMva_.push_back(XGBooster(edm::FileInPath("PhysicsTools/BPHNano/data/BDh_mva_May4/" + model + ".model").fullPath(),
                                        edm::FileInPath("PhysicsTools/BPHNano/data/BDh_mva_May4/" + model + ".features").fullPath()));

  }

  ~LambdabToLambdahhBuilder_v2() override {}
  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  float Bmva_estimator(pat::CompositeCandidate& Ldb0, unsigned int i);

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
  const double Lambdab_tkPtCut_;
  const double tkEtaCut_;
  const double tkChi2Cut_;
  const double tkIPSigXYCut_;

  const double vtxChi2Cut_;
  const double vtxDecaySigXYCut_;
  const double TrkSigXYCut_;
  const double vtxDecaySigXYZCut_;
  const double cosThetaXYCut_;
  const double cosThetaXYZCut_;
  const double diTrack2_dca_;
  const double Trk34SigXYCut_;
  
  const double ld0_l_xyzSigCut_;
  const double ldb0_PtCut_;
  const double ldb0_vtxDecaySigXYCut_;
  const double ld0_MassCut_;
  const double ldb0_MassCut_;
  const bool savetrack_;
  const int verbose;

  const bool isSignalMC=true;
  const bool onlyRecoMatchedPions=false;
  const bool onlyKeepGen=false;
};

float LambdabToLambdahhBuilder_v2::Bmva_estimator(pat::CompositeCandidate& Ldb0, unsigned int i){
       BMva_.at(i).set("Ldb0DiTrk1_dca",            Ldb0.userFloat("DiTrk1_dca"));
       return BMva_.at(i).predict();
}

void LambdabToLambdahhBuilder_v2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;

  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
  if (!theBeamSpotHandle.isValid()) {
    edm::LogError("LambdabToLambdahhBuilder_v2") << "No BeamSpot found!";
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
    edm::LogError("LambdabToLambdahhBuilder_v2") << "Track collections not found!";
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
    double ipsigXY_pv = std::abs(tmpTrack.dxy(referencePos) / tmpTrack.dxyError());
    if ( ipsigXY_pv < tkIPSigXYCut_) continue;
    double ipXY_bs = tmpTrack.dxy(theBeamSpotPos);
    double ipZ_bs  = tmpTrack.dz(theBeamSpotPos);
    double ipXY_pv = tmpTrack.dxy(referencePos);
    double ipZ_pv  = tmpTrack.dz(referencePos);
    double ipsigXY_bs = std::abs(tmpTrack.dxy(theBeamSpotPos) / tmpTrack.dxyError());
    double ipsigZ_bs  = std::abs(tmpTrack.dz(theBeamSpotPos) / tmpTrack.dzError());
    double ipsigZ_pv  = std::abs(tmpTrack.dz(referencePos) / tmpTrack.dzError());

    const reco::TransientTrack tmpTransient( (*trk.bestTrack()) , &theMagneticField);
    //distance closest approach in x,y wrt beam spot
    std::pair<double, double> DCA = computeDCA(tmpTransient, *theBeamSpot);
    float DCABS = DCA.first;
    float DCABSErr = DCA.second;
    float DCASig = (DCABSErr != 0 && float(DCABSErr) == DCABSErr) ? fabs(DCABS / DCABSErr) : -1;

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
        pcand.addUserInt("isPacked", (iTrk < nTracks));
        pcand.addUserInt("isLostTrk", (iTrk < nTracks) ? 0 : 1);
        pcand.addUserFloat("dxy",  trk.dxy());
        pcand.addUserFloat("dxySig", trk.dxy() / trk.dxyError()); // cut 2 in cmssw
        pcand.addUserFloat("dz",   trk.dz());
        pcand.addUserFloat("dzSig",  trk.dz() / trk.dzError()); // not cut
        pcand.addUserFloat("bt_pt",         bestTrack->pt());
        pcand.addUserFloat("DCASig", DCASig);

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

  //auto timer_trk = std::chrono::high_resolution_clock::now();

  // Get general di-track collection 
  // from Ld0, point to vec_ditrk
  std::vector<unsigned int> vec_ditrk_lambda;
  // from Ldb0, point to vec_ditrk
  std::vector<unsigned int> vec_ditrk_lambdab;

  // both Ld0, Ldb0
  std::vector<std::pair<unsigned int, unsigned int>> vec_ditrk;
  std::vector<std::tuple<double, float, float, double, double, double, double, float, double, double, double, double, double>> vec_ditrk_properties;

  std::unique_ptr<pat::CompositeCandidateCollection> ditracks_earlyout(new pat::CompositeCandidateCollection());

  for (unsigned int ditrk_leg1 = 0; ditrk_leg1 < vec_trk_ttrk.size(); ++ditrk_leg1) {

    const auto& diTrk1 = vec_trk_ttrk[ditrk_leg1].first;
    const auto& tdiTrk1 = vec_trk_ttrk[ditrk_leg1].second;
    
    for (unsigned int ditrk_leg2 = ditrk_leg1 + 1; ditrk_leg2 < vec_trk_ttrk.size(); ++ditrk_leg2) {
      
      const auto& diTrk2 = vec_trk_ttrk[ditrk_leg2].first;
      const auto& tdiTrk2 = vec_trk_ttrk[ditrk_leg2].second;

      // Opposite charge requirement
      if (diTrk1.charge() + diTrk2.charge() != 0) continue;
      // reorder, (tPosTrk, tNegTrk)
      const auto& tPosTrk = (diTrk1.charge() > 0) ? tdiTrk1 : tdiTrk2;
      const auto& tNegTrk = (diTrk1.charge() < 0) ? tdiTrk1 : tdiTrk2;
      unsigned int ditrkidx_pos = (diTrk1.charge() > 0) ? ditrk_leg1 : ditrk_leg2;
      unsigned int ditrkidx_neg = (diTrk1.charge() < 0) ? ditrk_leg1 : ditrk_leg2;

      //auto timer_1 = std::chrono::high_resolution_clock::now();

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
      double cxPtR2 = cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y();
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

      //auto timer_2 = std::chrono::high_resolution_clock::now();

      // step2: Pre-fit mPiPi mass cut
      double p1Mag2 = posTSCP.momentum().mag2();
      double p2Mag2 = negTSCP.momentum().mag2();
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();

      double totalE1 = std::sqrt(p1Mag2 + protonMass*protonMass) + std::sqrt(p2Mag2 + piMassSquared);
      double totalESq1 = totalE1 * totalE1;
      double massSquared1 = totalESq1 - totalPSq;

      double totalE2 = std::sqrt(p1Mag2 + piMassSquared) + std::sqrt(p2Mag2 + protonMass*protonMass);
      double totalESq2 = totalE2 * totalE2;
      double massSquared2 = totalESq2 - totalPSq;

      if(DCA_pos_beamspot.second==0 || DCA_neg_beamspot.second==0 || DCA_pos_pv.second==0 || DCA_neg_pv.second==0) continue;


      bool pass_ditrk_lambda = false;
      bool pass_ditrk_lambdab = false;
      // Good for ld0?
      if( ( abs(sqrt(massSquared1)-ld0_Mass) < 1.5*ld0_MassCut_ || abs(sqrt(massSquared2)-ld0_Mass) < 1.5*ld0_MassCut_ ) &&
          abs(DCA_pos_beamspot.first/DCA_pos_beamspot.second) > TrkSigXYCut_ &&
          abs(DCA_neg_beamspot.first/DCA_neg_beamspot.second) > TrkSigXYCut_) pass_ditrk_lambda = true;

      // Good for ldb0 ?
      if( dca < diTrack2_dca_ &&
          diTrk1.pt() > Lambdab_tkPtCut_ &&
	  diTrk2.pt() > Lambdab_tkPtCut_ && 	  
	  abs(DCA_pos_beamspot.first/DCA_pos_beamspot.second) > Trk34SigXYCut_ &&
          abs(DCA_neg_beamspot.first/DCA_neg_beamspot.second) > Trk34SigXYCut_
	  ) pass_ditrk_lambdab = true;

      if (pass_ditrk_lambda || pass_ditrk_lambdab){
          vec_ditrk.push_back(std::make_pair(ditrkidx_pos, ditrkidx_neg));
          vec_ditrk_properties.emplace_back(cxPtR2, cxPt.z(), trk1_dot_trk2, DCA_pos_beamspot.first, DCA_pos_beamspot.second, DCA_neg_beamspot.first, DCA_neg_beamspot.second, dca, massSquared1, DCA_pos_pv.first, DCA_pos_pv.second, DCA_neg_pv.first, DCA_neg_pv.second);
          if(savetrack_){
              pat::CompositeCandidate pcand;
              pcand.addUserInt("leg1_idx",   int(ditrkidx_pos));
              pcand.addUserInt("leg2_idx",   int(ditrkidx_neg));
	      pcand.addUserFloat("cxPtR2",         cxPtR2);
              pcand.addUserFloat("cxPtx",          cxPt.x());
              pcand.addUserFloat("cxPty",          cxPt.y());
              pcand.addUserFloat("cxPtz",          cxPt.z());
              pcand.addUserFloat("dot",            trk1_dot_trk2);
              pcand.addUserFloat("dca",            dca);
              pcand.addUserFloat("massSquared",    massSquared1);
              pcand.addUserFloat("trk1_bs_dca",    DCA_pos_beamspot.first);
              pcand.addUserFloat("trk2_bs_dca",    DCA_neg_beamspot.first);
              pcand.addUserFloat("trk1_pv_dca",    DCA_pos_pv.first);
              pcand.addUserFloat("trk2_pv_dca",    DCA_neg_pv.first);
              pcand.addUserFloat("trk1_bs_dcaErr", DCA_pos_beamspot.second);
              pcand.addUserFloat("trk2_bs_dcaErr", DCA_neg_beamspot.second);
              pcand.addUserFloat("trk1_pv_dcaErr", DCA_pos_pv.second);
              pcand.addUserFloat("trk2_pv_dcaErr", DCA_neg_pv.second);
	      ditracks_earlyout -> emplace_back(pcand);
          }
	  if(pass_ditrk_lambda)	  vec_ditrk_lambda.push_back( vec_ditrk.size()-1 );
          if(pass_ditrk_lambdab)       vec_ditrk_lambdab.push_back( vec_ditrk.size()-1 );
      }
      //std::chrono::duration<double> duration_1 = timer_2 - timer_1;
      //std::chrono::duration<double> duration_2 = timer_3 - timer_2;
      //std::chrono::duration<double> duration_3 = timer_4 - timer_3;
      //std::chrono::duration<double> duration_4 = timer_5 - timer_4;
      //std::cout << "duration_ditrk steps " << duration_1.count() << " "<< duration_2.count() << " " << duration_3.count() << " " << duration_4.count() << " " << std::endl;

    }
  }
  if(verbose>=5) std::cout<<"vec_trk_ttrk.size() "<< vec_trk_ttrk.size()<<" vec_ditrk.size() "<< vec_ditrk.size()<<std::endl;
  if(verbose>=5) 
	  std::cout<<"vec_ditrk_lambda.size() "<< vec_ditrk_lambda.size()<<" vec_ditrk_lambdab.size() "<< vec_ditrk_lambdab.size()<<std::endl;


//          auto timer_ditrk = std::chrono::high_resolution_clock::now();
//  if(timing){
//          std::chrono::duration<double> duration_trk = timer_trk - start;
//          std::chrono::duration<double> duration_ditrk = timer_ditrk - timer_trk;
//          std::cout << "duration_trk took " << duration_trk.count() << " seconds" << std::endl;
//          std::cout << "duration_ditrk took " << duration_ditrk.count() << " seconds" << std::endl;
//  }

  auto LooseLd0_out = std::make_unique<reco::VertexCompositePtrCandidateCollection>();
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  std::unique_ptr<TransientTrackCollection> trans_out( new TransientTrackCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> Ldb0_out(new pat::CompositeCandidateCollection());


  // Two ditracks make one ldb0
  for (unsigned int ditrk_idx1 = 0; ditrk_idx1 < vec_ditrk_lambda.size(); ++ditrk_idx1) { 
      
      // Kin fit for Ld0
      unsigned int ditrk_trueidx1 = vec_ditrk_lambda[ditrk_idx1];
      const auto& tTrk1 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx1].first].second;
      const auto& tTrk2 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx1].second].second;
      int mass_hyp = -1;
      for ( std::pair<double, double> masses : { std::pair<double, double>(protonMass, piMass), std::pair<double, double>(piMass, protonMass) } ) {
	  mass_hyp += 1;
          KinVtxFitter theKinVtxFitter(
          {tTrk1, tTrk2},
          {masses.first, masses.second},
          {0.00001, 0.00001} );
          if ( !theKinVtxFitter.success() ) continue;
          if ( theKinVtxFitter.chi2()<0 || theKinVtxFitter.dof()<0) continue;
          if ( theKinVtxFitter.prob()<0.001) continue;
          auto theKinVtxFitter_p4 = theKinVtxFitter.fitted_p4();
          if(abs(theKinVtxFitter_p4.mass()-ld0_Mass)>ld0_MassCut_) continue;
          if(abs(theKinVtxFitter.fitted_p4().eta())>2.4)  continue;
          if(verbose>=2) std::cout<< "Pass ld0 Kin fit" << std::endl;
          float Kin_mass = theKinVtxFitter.fitted_candidate().mass();
          float Kin_massErr = sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6));
          auto V0TT = theKinVtxFitter.fitted_candidate_ttrk();

          for (unsigned int ditrk_idx2 = 0; ditrk_idx2 < vec_ditrk_lambdab.size(); ++ditrk_idx2) {
              //std::cout<<"Mass before fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
              unsigned int ditrk_trueidx2 = vec_ditrk_lambdab[ditrk_idx2];
              if( ditrk_trueidx1 == ditrk_trueidx2 ) continue;
              if( vec_ditrk[ditrk_trueidx1].first == vec_ditrk[ditrk_trueidx2].first || vec_ditrk[ditrk_trueidx1].first == vec_ditrk[ditrk_trueidx2].second || vec_ditrk[ditrk_trueidx1].second == vec_ditrk[ditrk_trueidx2].first || vec_ditrk[ditrk_trueidx1].second == vec_ditrk[ditrk_trueidx2].second) continue;

              // A quick Mass fit
              const auto& tTrk3 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx2].first].second;
              const auto& tTrk4 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx2].second].second;
              auto mom1 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx2].first].first.momentum();
              auto mom2 = vec_trk_ttrk[vec_ditrk[ditrk_trueidx2].second].first.momentum();
              reco::Candidate::LorentzVector p4_1(mom1.x(), mom1.y(), mom1.z(), sqrt(mom1.mag2() + kplusMass * kplusMass));
              reco::Candidate::LorentzVector p4_2(mom2.x(), mom2.y(), mom2.z(), sqrt(mom2.mag2() + kplusMass * kplusMass));
              reco::Candidate::LorentzVector total_p4 = p4_1 + p4_2 + theKinVtxFitter_p4;
              float Ldb0_premass = total_p4.mass();
              if(abs(Ldb0_premass-Ldb0Mass)>ldb0_MassCut_*1.5) continue;
              if(total_p4.pt()<ldb0_PtCut_) continue;
              if(abs(total_p4.eta())>2.4) continue;
              if(verbose>=3) std::cout<< "Pass ldb0 presel" << std::endl;

              // Kin fit for Ldb0
              KinVtxFitter Ldb0_KinFitter(
                      {tTrk3, tTrk4, V0TT },
		      {kplusMass, kplusMass,  Kin_mass},
                      {K_SIGMA, K_SIGMA, Kin_massErr}
            	  );
              //std::cout<<"Mass after fit "<< theKinVtxFitter.fitted_candidate().mass()<<std::endl;
              if ( !Ldb0_KinFitter.success() ) continue;
              if ( Ldb0_KinFitter.chi2()<0 || Ldb0_KinFitter.dof()<0) continue;
              if ( Ldb0_KinFitter.prob()<0.001) continue;
              auto Ldb0_Kinfit_p4 = Ldb0_KinFitter.fitted_p4();
              if(abs(Ldb0_Kinfit_p4.mass()-Ldb0Mass)>ldb0_MassCut_) continue;    
              if(Ldb0_Kinfit_p4.pt() < ldb0_PtCut_) continue;
              auto Ldb0_pv_lxy  =  l_xy(Ldb0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              if (Ldb0_pv_lxy.error()!=0 && abs(Ldb0_pv_lxy.value()/Ldb0_pv_lxy.error()) < ldb0_vtxDecaySigXYCut_ ) continue;
              if((cos_theta_3D(theKinVtxFitter, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4()))<cosThetaXYCut_) continue;
              if(verbose>=2) std::cout<<" ldb0 looks good "<<std::endl;

              auto Ld0_ldb0_lxyz =  l_xyz(theKinVtxFitter, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z());
              double Ld0_Kin_ldb0_l_xyzSig = Ld0_ldb0_lxyz.error()==0? -99: abs(Ld0_ldb0_lxyz.value()/Ld0_ldb0_lxyz.error());
              if(Ld0_Kin_ldb0_l_xyzSig < ld0_l_xyzSigCut_) continue;
              // Save something
              pat::CompositeCandidate Ldb0;
              // idx
              Ldb0.addUserInt("DiTrack_idx1", int(ditrk_trueidx1));
              Ldb0.addUserInt("DiTrack_idx2", int(ditrk_trueidx2));
              Ldb0.addUserInt("Track_idx1",   int(vec_ditrk[ditrk_trueidx1].first));
              Ldb0.addUserInt("Track_idx2",   int(vec_ditrk[ditrk_trueidx1].second));
              Ldb0.addUserInt("Track_idx3",   int(vec_ditrk[ditrk_trueidx2].first));
              Ldb0.addUserInt("Track_idx4",   int(vec_ditrk[ditrk_trueidx2].second));
              // Ditrack variables
              // dca relative to BS or PV
              Ldb0.addUserFloat("DiTrk1_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_cxPtz",       std::get<1>(vec_ditrk_properties[ditrk_trueidx1]));		      
              Ldb0.addUserFloat("DiTrk1_dot",         std::get<2>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_dca",         abs(std::get<7>(vec_ditrk_properties[ditrk_trueidx1])));
              Ldb0.addUserFloat("DiTrk1_massSquared", std::get<8>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk1_mass", masses.first);
              Ldb0.addUserFloat("DiTrk1_trk2_mass", masses.second);

              // dca relative to bs and pv
              Ldb0.addUserFloat("DiTrk1_trk1_bs_dca", std::get<3>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk2_bs_dca", std::get<5>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk1_pv_dca", std::get<9>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk2_pv_dca", std::get<11>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk1_bs_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrk_trueidx1]))/std::get<4>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk2_bs_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrk_trueidx1]))/std::get<6>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk1_pv_dcaSig", abs(std::get<9>(vec_ditrk_properties[ditrk_trueidx1]))/std::get<10>(vec_ditrk_properties[ditrk_trueidx1]));
              Ldb0.addUserFloat("DiTrk1_trk2_pv_dcaSig", abs(std::get<11>(vec_ditrk_properties[ditrk_trueidx1]))/std::get<12>(vec_ditrk_properties[ditrk_trueidx1]));
              // displacement relative to bs and pv

              Ldb0.addUserFloat("DiTrk2_cxPtR2",      std::get<0>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_cxPtz",       std::get<1>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_dot",         std::get<2>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_dca",         abs(std::get<7>(vec_ditrk_properties[ditrk_trueidx2])));
              Ldb0.addUserFloat("DiTrk2_massSquared", std::get<8>(vec_ditrk_properties[ditrk_trueidx2]));
              // dca relative to bs and pv
              Ldb0.addUserFloat("DiTrk2_trk1_bs_dca", std::get<3>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk2_bs_dca", std::get<5>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk1_pv_dca", std::get<9>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk2_pv_dca", std::get<11>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk1_bs_dcaSig", abs(std::get<3>(vec_ditrk_properties[ditrk_trueidx2]))/std::get<4>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk2_bs_dcaSig", abs(std::get<5>(vec_ditrk_properties[ditrk_trueidx2]))/std::get<6>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk1_pv_dcaSig", abs(std::get<9>(vec_ditrk_properties[ditrk_trueidx2]))/std::get<10>(vec_ditrk_properties[ditrk_trueidx2]));
              Ldb0.addUserFloat("DiTrk2_trk2_pv_dcaSig", abs(std::get<11>(vec_ditrk_properties[ditrk_trueidx2]))/std::get<12>(vec_ditrk_properties[ditrk_trueidx2]));

              // Ld0
              // Basic
              Ldb0.addUserFloat("Ld0_Kin_vtx_x", theKinVtxFitter.fitted_vtx().x());
              Ldb0.addUserFloat("Ld0_Kin_vtx_y", theKinVtxFitter.fitted_vtx().y());
              Ldb0.addUserFloat("Ld0_Kin_vtx_r", sqrt( theKinVtxFitter.fitted_vtx().x()*theKinVtxFitter.fitted_vtx().x() + theKinVtxFitter.fitted_vtx().y()*theKinVtxFitter.fitted_vtx().y() ) );
              Ldb0.addUserFloat("Ld0_Kin_vtx_z", theKinVtxFitter.fitted_vtx().z());
              Ldb0.addUserFloat("Ld0_Kin_chi2",  theKinVtxFitter.chi2());
              Ldb0.addUserFloat("Ld0_Kin_dof",   theKinVtxFitter.dof()); // always 1 ?
              Ldb0.addUserFloat("Ld0_Kin_prob",  theKinVtxFitter.prob());
              Ldb0.addUserFloat("Ld0_Kin_pt",    theKinVtxFitter.fitted_p4().pt());
              Ldb0.addUserFloat("Ld0_Kin_eta",   theKinVtxFitter.fitted_p4().eta());
              Ldb0.addUserFloat("Ld0_Kin_phi",   theKinVtxFitter.fitted_p4().phi());
              Ldb0.addUserFloat("Ld0_Kin_mass",  theKinVtxFitter.fitted_candidate().mass());
              Ldb0.addUserInt("Ld0_Kin_mass_hyp",  mass_hyp);
              Ldb0.addUserFloat("Ld0_Kin_massErr", sqrt(theKinVtxFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Ldb0.addUserFloat("Ld0_Kin_trk1_pt", theKinVtxFitter.daughter_p4(0).pt());
              Ldb0.addUserFloat("Ld0_Kin_trk1_p", theKinVtxFitter.daughter_p4(0).mag());
              Ldb0.addUserFloat("Ld0_Kin_trk1_eta", theKinVtxFitter.daughter_p4(0).eta());
              Ldb0.addUserFloat("Ld0_Kin_trk1_phi", theKinVtxFitter.daughter_p4(0).phi());
              Ldb0.addUserFloat("Ld0_Kin_trk1_mass", theKinVtxFitter.daughter_p4(0).mass());
              Ldb0.addUserFloat("Ld0_Kin_trk2_pt", theKinVtxFitter.daughter_p4(1).pt());
              Ldb0.addUserFloat("Ld0_Kin_trk2_p", theKinVtxFitter.daughter_p4(1).mag());
              Ldb0.addUserFloat("Ld0_Kin_trk2_eta", theKinVtxFitter.daughter_p4(1).eta());
              Ldb0.addUserFloat("Ld0_Kin_trk2_phi", theKinVtxFitter.daughter_p4(1).phi());
              Ldb0.addUserFloat("Ld0_Kin_trk2_mass", theKinVtxFitter.daughter_p4(1).mass());
              // pointing angle relative to bs, pv, Ldb0
              Ldb0.addUserFloat("Ld0_Kin_bs_alpha_2D", (cos_theta_2D(theKinVtxFitter, *theBeamSpot, theKinVtxFitter.fitted_p4())));
              Ldb0.addUserFloat("Ld0_Kin_pv_alpha_2D", (cos_theta_2D(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z(), theKinVtxFitter.fitted_p4())));
              Ldb0.addUserFloat("Ld0_Kin_pv_alpha_3D", (cos_theta_3D(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z(), theKinVtxFitter.fitted_p4())));
              Ldb0.addUserFloat("Ld0_Kin_ldb0_alpha_2D", (cos_theta_2D(theKinVtxFitter, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4())));
              Ldb0.addUserFloat("Ld0_Kin_ldb0_alpha_3D", (cos_theta_3D(theKinVtxFitter, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z(), theKinVtxFitter.fitted_p4())));
              // displacement relative to bs, pv, Ldb0
              // dca relative to bs and pv
              auto Ld0_bs_lxy  =  l_xy(theKinVtxFitter, *theBeamSpot);
              auto Ld0_pv_lxy  =  l_xy(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto Ld0_pv_lxyz =  l_xyz(theKinVtxFitter, referencePos.x(), referencePos.y(), referencePos.z());
              auto Ld0_ldb0_lxy  =  l_xy(theKinVtxFitter, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z());
              double Ld0_Kin_bs_l_xySig  = Ld0_bs_lxy.error()==0?  -99: abs(Ld0_bs_lxy.value()/Ld0_bs_lxy.error());
              double Ld0_Kin_pv_l_xySig  = Ld0_pv_lxy.error()==0?  -99: abs(Ld0_pv_lxy.value()/Ld0_pv_lxy.error());
              double Ld0_Kin_pv_l_xyzSig = Ld0_pv_lxyz.error()==0? -99: abs(Ld0_pv_lxyz.value()/Ld0_pv_lxyz.error());
              double Ld0_Kin_ldb0_l_xySig  = Ld0_ldb0_lxy.error()==0?  -99: abs(Ld0_ldb0_lxy.value()/Ld0_ldb0_lxy.error());
              Ldb0.addUserFloat("Ld0_Kin_bs_l_xy",     Ld0_bs_lxy.value());
              Ldb0.addUserFloat("Ld0_Kin_bs_l_xySig",  Ld0_Kin_bs_l_xySig);
              Ldb0.addUserFloat("Ld0_Kin_pv_l_xy",     Ld0_pv_lxy.value());
              Ldb0.addUserFloat("Ld0_Kin_pv_l_xySig",  Ld0_Kin_pv_l_xySig);
              Ldb0.addUserFloat("Ld0_Kin_pv_l_xyz",    Ld0_pv_lxyz.value());
              Ldb0.addUserFloat("Ld0_Kin_pv_l_xyzSig", Ld0_Kin_pv_l_xyzSig);
              Ldb0.addUserFloat("Ld0_Kin_ldb0_l_xy",     Ld0_ldb0_lxy.value());
              Ldb0.addUserFloat("Ld0_Kin_ldb0_l_xySig",  Ld0_Kin_ldb0_l_xySig);
              Ldb0.addUserFloat("Ld0_Kin_ldb0_l_xyz",    Ld0_ldb0_lxyz.value());
              Ldb0.addUserFloat("Ld0_Kin_ldb0_l_xyzSig", Ld0_Kin_ldb0_l_xyzSig);
              // dca relative to bs and pv
              std::pair<double, double> DCA_ld0 = computeDCA(V0TT, Ldb0_KinFitter.fitted_vtx().x(), Ldb0_KinFitter.fitted_vtx().y(), Ldb0_KinFitter.fitted_vtx().z());
              double Ld0_Kin_ldb0_dcaSig = DCA_ld0.second==0? -99: abs(DCA_ld0.first/DCA_ld0.second);
              Ldb0.addUserFloat("Ld0_Kin_ldb0_dca", DCA_ld0.first);
              Ldb0.addUserFloat("Ld0_Kin_ldb0_dcaSig", Ld0_Kin_ldb0_dcaSig);


              // Ldb0
              // Basic
              Ldb0.addUserFloat("Ldb0_premass", Ldb0_premass);
              Ldb0.addUserFloat("Ldb0_Kin_vtx_x", Ldb0_KinFitter.fitted_vtx().x());
              Ldb0.addUserFloat("Ldb0_Kin_vtx_y", Ldb0_KinFitter.fitted_vtx().y());
              Ldb0.addUserFloat("Ldb0_Kin_vtx_r", sqrt( Ldb0_KinFitter.fitted_vtx().x()*Ldb0_KinFitter.fitted_vtx().x() + Ldb0_KinFitter.fitted_vtx().y()*Ldb0_KinFitter.fitted_vtx().y() ) );
              Ldb0.addUserFloat("Ldb0_Kin_vtx_z", Ldb0_KinFitter.fitted_vtx().z());
              Ldb0.addUserFloat("Ldb0_Kin_chi2", Ldb0_KinFitter.chi2());
              Ldb0.addUserFloat("Ldb0_Kin_dof", Ldb0_KinFitter.dof());
              Ldb0.addUserFloat("Ldb0_Kin_prob", Ldb0_KinFitter.prob());
              Ldb0.addUserFloat("Ldb0_Kin_pt", Ldb0_Kinfit_p4.pt() );
              Ldb0.addUserFloat("Ldb0_Kin_eta", Ldb0_Kinfit_p4.eta() );
              Ldb0.addUserFloat("Ldb0_Kin_phi", Ldb0_Kinfit_p4.phi() );
              Ldb0.addUserFloat("Ldb0_Kin_mass", Ldb0_Kinfit_p4.mass() );
              Ldb0.addUserFloat("Ldb0_Kin_massErr", sqrt(Ldb0_KinFitter.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
              Ldb0.addUserFloat("Ldb0_Kin_trk3_pt", Ldb0_KinFitter.daughter_p4(0).pt());
              Ldb0.addUserFloat("Ldb0_Kin_trk3_eta", Ldb0_KinFitter.daughter_p4(0).eta());
              Ldb0.addUserFloat("Ldb0_Kin_trk3_phi", Ldb0_KinFitter.daughter_p4(0).phi());
              Ldb0.addUserFloat("Ldb0_Kin_trk4_pt", Ldb0_KinFitter.daughter_p4(1).pt());
              Ldb0.addUserFloat("Ldb0_Kin_trk4_eta", Ldb0_KinFitter.daughter_p4(1).eta());
              Ldb0.addUserFloat("Ldb0_Kin_trk4_phi", Ldb0_KinFitter.daughter_p4(1).phi());
              Ldb0.addUserFloat("Ldb0_Kin_ld0_pt", Ldb0_KinFitter.daughter_p4(2).pt());
              Ldb0.addUserFloat("Ldb0_Kin_ld0_eta", Ldb0_KinFitter.daughter_p4(2).eta());
              Ldb0.addUserFloat("Ldb0_Kin_ld0_phi", Ldb0_KinFitter.daughter_p4(2).phi());
              // pointing angle relative to bs, pv, Ldb0
              Ldb0.addUserFloat("Ldb0_Kin_bs_alpha_2D", (cos_theta_2D(Ldb0_KinFitter, *theBeamSpot, Ldb0_KinFitter.fitted_p4())));
              Ldb0.addUserFloat("Ldb0_Kin_pv_alpha_2D", (cos_theta_2D(Ldb0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), Ldb0_KinFitter.fitted_p4())));
              Ldb0.addUserFloat("Ldb0_Kin_pv_alpha_3D", (cos_theta_3D(Ldb0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z(), Ldb0_KinFitter.fitted_p4())));
              // displacement relative to bs, pv, Ldb0
              auto Ldb0_bs_lxy  =  l_xy(Ldb0_KinFitter, *theBeamSpot);
              auto Ldb0_pv_lxyz =  l_xyz(Ldb0_KinFitter, referencePos.x(), referencePos.y(), referencePos.z());
              double Ldb0_Kin_bs_l_xySig  = Ldb0_bs_lxy.error()==0?  -99: abs(Ldb0_bs_lxy.value()/Ldb0_bs_lxy.error());
              double Ldb0_Kin_pv_l_xySig  = Ldb0_pv_lxy.error()==0?  -99: abs(Ldb0_pv_lxy.value()/Ldb0_pv_lxy.error());
              double Ldb0_Kin_pv_l_xyzSig = Ldb0_pv_lxyz.error()==0? -99: abs(Ldb0_pv_lxyz.value()/Ldb0_pv_lxyz.error());
              Ldb0.addUserFloat("Ldb0_Kin_bs_l_xy",     Ldb0_bs_lxy.value());
              Ldb0.addUserFloat("Ldb0_Kin_bs_l_xySig",  Ldb0_Kin_bs_l_xySig);
              Ldb0.addUserFloat("Ldb0_Kin_pv_l_xy",     Ldb0_pv_lxy.value());
              Ldb0.addUserFloat("Ldb0_Kin_pv_l_xySig",  Ldb0_Kin_pv_l_xySig);
              Ldb0.addUserFloat("Ldb0_Kin_pv_l_xyz",    Ldb0_pv_lxyz.value());
              Ldb0.addUserFloat("Ldb0_Kin_pv_l_xyzSig", Ldb0_Kin_pv_l_xyzSig);
              Ldb0.setP4(Ldb0_KinFitter.fitted_p4());

              // dca relative to bs and pv
              if(std::isnan(abs(Ldb0.userFloat("Ld0_Kin_ldb0_l_xySig")))) continue;
              if(std::isnan(abs(Ldb0.userFloat("Ld0_Kin_ldb0_dcaSig")))) continue;


              // Kin fit for Ldb0 with lambda0 set to world average
              KinVtxFitter Ldb0_KinFitter_fixmass(
                      {tTrk3, tTrk4, V0TT },
                      {kplusMass, kplusMass,  ld0_Mass},
                      {K_SIGMA, K_SIGMA, K_SIGMA}
                  );
	      if (Ldb0_KinFitter_fixmass.success()){
		      Ldb0.addUserFloat("Ldb0_Kin_mass_fixlb0", Ldb0_KinFitter_fixmass.fitted_p4().mass());
                      Ldb0.addUserFloat("Ldb0_Kin_massErr_fixlb0", sqrt(Ldb0_KinFitter_fixmass.fitted_candidate().kinematicParametersError().matrix()(6, 6)));
	      }else{
		      Ldb0.addUserFloat("Ldb0_Kin_mass_fixlb0", -1);
                      Ldb0.addUserFloat("Ldb0_Kin_massErr_fixlb0", -1);
	      }
  
	      //compute isolation
              std::vector<std::string> dnames{ "trk3", "trk4", "ld0" };
	      std::vector<float> iso(dnames.size(), 0);
              unsigned int nTracks = tracks->size();
              unsigned int totalTracks = nTracks + lostTracks->size();
              for (unsigned int iTrk = 0; iTrk < totalTracks; ++iTrk) {
                const pat::PackedCandidate& trk = (iTrk < nTracks) ? (*tracks)[iTrk] : (*lostTracks)[iTrk - nTracks];
                if (!trk.hasTrackDetails() || !trk.bestTrack()) continue;
                if (fabs(trk.pdgId()) != 211) continue; //do we want also to keep muons?
                if (iTrk < nTracks && !trk.trackHighPurity()) continue;
		const reco::Track* bestTrack = trk.bestTrack();
                for ( size_t iname = 0; iname < dnames.size(); ++iname) {
                  float dr = deltaR( Ldb0.userFloat("Ldb0_Kin_" + dnames[iname] + "_eta"), Ldb0.userFloat("Ldb0_Kin_" + dnames[iname] + "_phi"), bestTrack->eta(), bestTrack->phi());
                  if (dr > 0 && dr < 0.4)
                    iso[iname] += bestTrack->pt();
                }
              }
              for (size_t idaughter = 0; idaughter < dnames.size(); idaughter++) {
                Ldb0.addUserFloat("Ldb0_Kin_" + dnames[idaughter] + "_iso04", iso[idaughter]);
              }

              for (unsigned int i=0; i < BMva_.size(); ++i) {
                  Ldb0.addUserFloat("xgb_" + xgboost_variable_names_.at(i), Bmva_estimator(Ldb0, i));
              }
              //if(Ldb0.userFloat("xgb_" + xgboost_variable_names_.at(0))<0.2) continue;
              Ldb0_out->push_back(Ldb0);
	      if(verbose>=5) std::cout<<"Save a good candidate \n"<<std::endl;
          }
      }
  }

//        auto end = std::chrono::high_resolution_clock::now();
//  if(timing){
//      std::chrono::duration<double> duration_loop = end - timer_ditrk;
//      std::cout << "recoB took " << duration_loop.count() << " seconds" << std::endl;
//  }

  if(savetrack_) iEvent.put(std::move(tracks_earlyout),       "SelectedTracks");    
  if(savetrack_) iEvent.put(std::move(ditracks_earlyout),       "SelectedDiTracks");
  iEvent.put(std::move(Ldb0_out), "Ldb0");
  return;

}

DEFINE_FWK_MODULE(LambdabToLambdahhBuilder_v2);

