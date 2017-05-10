#include "particles_boundary.hpp"

// cycle through all particles
int BC_Particle::computeParticleBCs(std::vector<Particle> *pl) {
  long size = pl->size();
  for(long i = 0; i < size; i++){
    (*pl)[i].isGhost = (*pl)[i].isGhost || 
      particle_BC(&(*pl)[i]);
  }
  return completeBC(pl);
}
