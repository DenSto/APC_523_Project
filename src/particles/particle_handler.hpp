#ifndef PARTICLE_HANDLER_HPP
#define PARTICLE_HANDLER_HPP

//! Class that handles all particle-relevant operations.
/*!
        Particle handler handles all the particle operations. This includes deposition,
        boundary conditions, particle pushing, and communication between MPI nodes if
        needed
*/
#include <vector>
#include <stdio.h>
#include "../IO/input.hpp"
#include "../domain/domain.hpp"
#include "../grid/grid.hpp"
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../boundaries/particles_boundary.hpp"

class Particle_Handler {
public:
  Particle_Handler(); // list of np particles and their fields
  ~Particle_Handler();

  void Load(Input_Info_t *input_info, Domain* domain, Grid* grid); // Initialize particles
  void Push(double dt);   // Push all particles
  long nParticles();

  void incrementNParticles(int inc);

  void SortParticles(Particle_Compare comp); // quicksort particle list
  void CountingSortParticles(Grid* grid); // quicksort particle list

  void setPusher(Pusher* pusher) {pusher_=pusher;};
  void clearGhosts();  // remove all ghost particles in particle list

  void InterpolateEB(Grid* grid);
  void depositRhoJ(Grid *grid); // deposit current and charge density from particles to grid

  double computeCFLTimestep(Grid* grid); // return timestep computed from max velocity and grid size


  void setParticleBoundaries(BC_Particle** bc){boundaries_=bc;}
  void executeParticleBoundaryConditions(Grid* grid);

  void outputParticles(const char* basename, long nstep, Input_Info_t *input_info); //should be in its own class.

#ifdef PART_IN_CELL
  void updateCellLists(Grid* grid, short justGhosts);
  void insertParticle(Particle p,long id);
#endif

private:
  BC_Particle** boundaries_; /* Particle Boundary Conditions */
  Pusher* pusher_;
  int* pa_;
  int* pa_save_;
  int ncell_;


   // Output parameters (should be in its own class)
  double dT_, nextT_; // variables for time cadencing
  long dstep_, nextStep_; //variables for step cadencing
  long outputCount_; // How many particles OF EACH SPECIES per core to output
//  int lz1,lz2;  // leading zeros for output filename

  long np_;                      /* total number of particles */
#ifdef PART_IN_CELL
  std::vector<Particle>* parts_;   /* array of particle vectors, size of NCell */
#else
  std::vector<Particle> parts_;    /* Vector of particles */
#endif

};

#endif
