// original author: RK18 team
#include "KinVtxFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" // MIGHT be useful for Phi->KK?
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas):
  n_particles_{masses.size()} {
  
  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;
  for(size_t i = 0; i < tracks.size(); ++i) {
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  KinematicParticleVertexFitter kcv_fitter;    
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles);

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) { 
    success_=false; 
    return;
  }
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
}



KinVtxFitter::KinVtxFitter(const std::vector<reco::TransientTrack> tracks, 
                           const std::vector<double> masses, 
                           std::vector<float> sigmas, ParticleMass dilep_mass):
  n_particles_{masses.size()} {
  
  KinematicParticleFactoryFromTransientTrack factory;
  std::vector<RefCountedKinematicParticle> particles;
  for(size_t i = 0; i < tracks.size(); ++i) {
    particles.emplace_back(
      factory.particle(
        tracks.at(i), masses.at(i), kin_chi2_, 
        kin_ndof_, sigmas[i]
        )
      );
  }

  MultiTrackKinematicConstraint *  dilep_const;
 
  if(tracks.size()==2){
      MultiTrackKinematicConstraint * dilep_const = new TwoTrackMassKinematicConstraint(dilep_mass);
  }else{
      //dilep_const = new MultiTrackKinematicConstraint(dilep_mass);
      KinematicConstraint * dilep_const = new MassKinematicConstraint(dilep_mass, 1e-6);
  }
  KinematicConstrainedVertexFitter kcv_fitter;    
  RefCountedKinematicTree vtx_tree = kcv_fitter.fit(particles,dilep_const);

  if (vtx_tree->isEmpty() || !vtx_tree->isValid() || !vtx_tree->isConsistent()) {
    success_ = false; 
    return;
  }

  vtx_tree->movePointerToTheTop(); 
  fitted_particle_ = vtx_tree->currentParticle();
  fitted_vtx_ = vtx_tree->currentDecayVertex();
  if (!fitted_particle_->currentState().isValid() || !fitted_vtx_->vertexIsValid()){ 
    success_ = false; 
    return;
  }
  fitted_state_ = fitted_particle_->currentState();
  fitted_children_ = vtx_tree->finalStateParticles();
  if(fitted_children_.size() != n_particles_) { 
    success_=false; 
    return;
  }
  fitted_track_ = fitted_particle_->refittedTransientTrack();
  success_ = true;
}

