#include "../globals.hpp"
#include "particle_handler.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../pusher/boris.hpp"
#include "../grid/grid.hpp"
#include "../utils/RNG.hpp"
#include <string.h>
#include <iostream>
#include <limits>
#define _USE_MATH_DEFINES
#if USE_MPI
#include "mpi.h"
#endif

#define CALC_ID -10000

Particle_Handler::Particle_Handler(){
  np_=0;
  dT_=-1;
  dstep_=-1;
  nextT_=0;
  nextStep_=0;
}

Particle_Handler::~Particle_Handler(){
}

//! Load and initialize the particle handler. Should be called at the beginning of the run.
void Particle_Handler::Load(Input_Info_t *input_info, Domain* domain, Grid* grid){

  int restart = input_info->restart;
  double* L = domain->getLxyz();
  double* x0 = domain->getxyz0();
  double idx[3]; grid->getInvSpacing(idx);

  int* n= domain->getnxyz();

  int cpt = input_info->cellPerTile;

  for(int i = 0; i < 3; i++){
    cpt_[i] = cpt;
    idx_[i] = idx[i] / cpt;
    x0_[i] = x0[i];
    xMax_[i] = x0_[i] + L[i];
    nT_[i] = (int)ceil(((double)n[i])/((double)cpt));
  }
  nTiles_ = nT_[0]*nT_[1]*nT_[2];

  parts_ = new std::vector<Particle>[nTiles_+1];

  pa_= (int*) malloc(sizeof(int) *(cpt*cpt*cpt  +1));
  pa_save_ = (int*) malloc(sizeof(int) *( cpt*cpt*cpt + 1));

  int npart      = input_info->nparticles_tot*n[0]*n[1]*n[2];
  double *mass   = input_info->mass_ratio;
  double *charge = input_info->charge_ratio;
  double *dens   = input_info->dens_frac; 


  if(restart==0){// initial run
    if(rank_MPI==0)printf("    Loading random particles from distribution...\n");
    Random_Number_Generator *rng = new Random_Number_Generator(-1);
    int ispec = 0; // temporaty counter  
    double cden = dens[0]; // cummulative density fraction
    double vth;
    for(long ip=0; ip < npart; ip++){
      Particle p = new_particle();
      if(ip >= cden*npart){
        ispec += 1;
        cden  += dens[ispec];
      }
      p.q = charge[ispec]; // super particle charge
      p.m = mass[ispec];   // super particle mass
      p.type = ispec;

      vth=UNIT_VTH*sqrt(input_info->temp[ispec]/p.m);

      p.x[0]=rng->getUniform()*L[0]+x0[0];
      p.x[1]=rng->getUniform()*L[1]+x0[1];
      p.x[2]=rng->getUniform()*L[2]+x0[2];

      p.v[0]=rng->getGaussian(0.0,vth);
      p.v[1]=rng->getGaussian(0.0,vth);
      p.v[2]=rng->getGaussian(0.0,vth);

      p.my_id=ip;
      p.initRank=rank_MPI;

      int id = getTileID(p.x[0],p.x[1],p.x[2]);
      parts_[id].push_back(p);
      np_++;
    } 
  } else {
    if(rank_MPI==0)printf("    Loading particles from file...\n");
    Particle p = new_particle();
    p.q = 1.0;
    p.m = 1.0;
    p.x[0]=L[0]/2+x0[0];
    p.x[1]=L[1]/2+x0[1];
    p.x[2]=L[2]/2+x0[2];
    p.v[0]=0.01;
    p.v[1]=0.001;
    p.v[2]=0.001;
    int id = getTileID(p.x[0],p.x[1],p.x[2]);
    parts_[id].push_back(p);
    np_++;
  }
  printf("N Particles %ld\n",np_);
}

void Particle_Handler::Push(double dt){
  for(long i = 0; i < nTiles_; i++){
    long size = (long) parts_[i].size();
#pragma simd
    for(long ip=0;ip < size; ip++){
      double qom = parts_[i][ip].q/parts_[i][ip].m;
      double ret[6];
      pusher_->Step(parts_[i][ip].x,parts_[i][ip].v,qom, parts_[i][ip].field,dt,ret);

      parts_[i][ip].x[0] = ret[0];
      parts_[i][ip].x[1] = ret[1];
      parts_[i][ip].x[2] = ret[2];
      parts_[i][ip].v[0] = ret[3];
      parts_[i][ip].v[1] = ret[4];
      parts_[i][ip].v[2] = ret[5];
    }
  }
}

long Particle_Handler::nParticles(){
  return -1;
}

void Particle_Handler::incrementNParticles(int inc){
#ifndef PART_IN_CELL
  np_+=inc;
#endif
}

void Particle_Handler::InterpolateEB(Grid* grid){
  double pos[3]; //Vector of position of particle.
  double icell[3]; //Vector of inverse lengths of unit cell.
  double L0[3]; //Vector of lengths of unit cell.
  double weight[3][3][3];


  int is,js,ks;
  int ic,jc,kc;

  //Get lengths of grid cells.
  grid->getInvSpacing(icell);
  grid->getL0(L0);

  //int count =0;

  for (long j = 0; j < nTiles_; j++){
    long size = (long) parts_[j].size();
#pragma vector
    for (long i=0; i < size; i++) {
      //Get position of particle.
      pos[0] = parts_[j][i].x[0];
      pos[1] = parts_[j][i].x[1];
      pos[2] = parts_[j][i].x[2];

      getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);
      Couple cp;
      grid->getCouple(is,js,ks,&cp);
      interpolateFields(&cp,weight,&parts_[j][i].field);
    }
  }
}

void Particle_Handler::depositRhoJ(Grid *grid){
  double pos[3]; //Vector of position of particle.
  double icell[3]; //Vector of inverse lengths of unit cell.
  double L0[3]; //Vector of lengths of unit cell.
  double weight[3][3][3];

  if(useVecDeposition){
    depositRho(grid,parts_);
    return;
  }


  double ***rho = grid->getRho();
  int nghost = grid->getGhost();

  //Get lengths of grid cells.
  grid->getInvSpacing(icell);
  grid->getL0(L0);

  for (long n=0; n < nTiles_ ; n++) {
    long size = (long) parts_[n].size();
    for (long ip=0; ip< size; ip++) {
      int is,js,ks;
      int ic,jc,kc;
 //     Particle *p = &parts_[n][ip];
  
      //Get position of particle.
      pos[0] = parts_[n][ip].x[0];
      pos[1] = parts_[n][ip].x[1];
      pos[2] = parts_[n][ip].x[2];
      double q = parts_[n][ip].q;

      getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);

      for( int i = 0; i < 3; i++){
        for( int j = 0; j < 3; j++){
#pragma simd
          for( int k = 0; k < 3; k++){
            rho[is +nghost + i][js + nghost + j][ks + nghost + k] += q*weight[i][j][k];	    
          }
        }
      }
    }
  }
}

//! Sort particles based on grid location. 
/*
 * Sorts particles based on grid location using std::sort which is O(n log n).
 * This should be called often to ensure cache hits.  Should be emperically determined.
 */
void Particle_Handler::SortParticles(Particle_Compare comp){
  if(debug>0 && rank_MPI == 0) fprintf(stderr,"Sorting Particles.\n");
  for(int i = 0; i < nTiles_; i++)
    std::sort(parts_[i].begin(),parts_[i].end(),comp);
}

void Particle_Handler::CountingSortParticles(Grid* grid){
  long i,j,k,se;
  int cellid;
  double pos[3];

  Particle *in, *out, *stop, *tmp;
  Particle tmp1;
  tmp = &tmp1;

  for(i = 0; i < ncell_; i++){
    pa_[i] = 0;
  }
  for (i= 0; i<np_; i++) {
    pos[0] = parts_[i].x[0];
    pos[1] = parts_[i].x[1];
    pos[2] = parts_[i].x[2];
    cellid = grid->getCellID(pos[0],pos[1],pos[2]);
    pa_[cellid]++;  
  }
  for (i=k=0; i<=ncell_; i++) { 
    j = pa_[i]; 
    pa_save_[i] = pa_[i] = k; 
    k+=j;
  }

  i = 0;

  while( i < ncell_){
    if(pa_[i] >= pa_save_[i+1]) {
      i++;
    } else {
      in = stop = &parts_[pa_[i]];
      do{
        pos[0] = in->x[0];
        pos[1] = in->x[1];
        pos[2] = in->x[2];
        cellid = grid->getCellID(pos[0],pos[1],pos[2]);
        out = &parts_[pa_[cellid]];
        pa_[cellid]++;

        if( out != stop){
          tmp1 = *out;
          *out = *in;
          *in=tmp1;
        } else if( out != in) {
          *out=*in;
        }
      } while( out != stop);
    }
  }
}


double Particle_Handler::computeCFLTimestep(Grid* grid){
  double dx[3];

  double mindt = 1e20;

  double v;
  grid->getSpacing(dx);
  for (int i= 0; i<nTiles_; i++) {
    long size = (long) parts_[i].size();
    for (int ip= 0; ip < size; ip++) {
      for(int j = 0; j < 3; j++){
        v = parts_[i][ip].v[j];
        if(v*v > 0){
          if(fabs(dx[j]/v) < mindt) mindt = fabs(dx[j]/v);
        }
      } 
    }
  }
  return mindt;

#if USE_MPI
//  double mindtall[3];
//  int ierr = MPI_Allreduce(mindt,mindtall,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  //return *std::min_element(mindtall,mindtall+3);
  return *std::min_element(mindt,mindt+3);
#endif
}

//! Clear all ghost particles. Uses a swap-to-back and pop-last-element for speed.
void Particle_Handler::clearGhosts(){
  long size = 0;
  for(int i = 0; i < nTiles_; i++) size += (long) parts_[i].size();
  parts_[nTiles_].clear();
  assert(size == np_);
}

void Particle_Handler::updateTiles(Grid *grid,short justGhosts){
  long i = justGhosts ? nTiles_ : 0;
  for(; i < nTiles_ + 1; i++){
    for(std::vector<Particle>::iterator iter = parts_[i].begin(); iter != parts_[i].end();){
      long id = getTileID_wGhost(iter->x[0],iter->x[1],iter->x[2]);
      id = id < 0 ? nTiles_ : id;
      if(id != i){ // particle is not in it's cell anymore
        // swap pop
        Particle pt = *iter;
        std::swap(*iter, parts_[i].back());
        parts_[i].pop_back();

        insertParticle(pt,id);
      } else {
        iter++;
      }
    }
  }
}

void Particle_Handler::insertParticle(Particle p, long id){
  long ind = id < 0 ? nTiles_ : id;
  parts_[ind].push_back(p);
}


void Particle_Handler::executeParticleBoundaryConditions(Grid* grid){
  for(int i = 0; i < 6; i++){
    // determine whether particles are ghost
        // place ghost particles
        // change the number of particles np_ in each domain
    boundaries_[i]->computeParticleBCs(&parts_[nTiles_]);
    updateTiles(grid,1);
  }
}

//! Output particles
/*!
 *  Output particles. Currently outputs time, position and velocity.
 *
 *  Can work with either cadencing on time (output every dT) or 
 *  cadencing on steps (output ever dsteps), or both. Either are 
 *  optional parameters in the input file and will default to -1.
 *
 *  Particles are written to the same directory as the program input file and 
 *  named with initial rank and id.
 */
void Particle_Handler::outputParticles(const char* basename,long step, Input_Info_t *input_info){

  if(debug>2 && rank_MPI==0) printf("    ti=%ld: writing particle tracks...\n",step);

  double t = time_phys;
  dstep_ = input_info->nstep_parts;
  outputCount_ = input_info->output_pCount;

  static bool init = true;
  bool needsOutput=false;
  if(init){
    init=false;
    //mkdir("./tracks", 0775); //create the particle directory
    for(int i = 0; i < nTiles_; i++){
      for(std::vector<Particle>::iterator iter = parts_[i].begin();iter!=parts_[i].end();++iter){
        if(iter->my_id < outputCount_){
          char fname[100];
          sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);  
          FILE *pout=fopen(fname,"w");
          fprintf(pout,"[1] time [2] x [3] y [4] z  [5] vx   [6] vy   [7] vz\n");
          fclose(pout);
        }
      }
    }
  }

  // cadence on i or t
  if(dT_ > 0 && t  >= nextT_){
    needsOutput=true;
    nextT_ += dT_;
  }
  
  if(dstep_ > 0 && step >= nextStep_){
    needsOutput=true;
    nextStep_+=dstep_;
  }

  if(!needsOutput)
    return;

 if(debug){
   fprintf(stderr,"rank=%d: writing tracks for %ld particles...\n",
           rank_MPI,outputCount_);
  }
  char fname[100];
  FILE *pout;
  for(int i =0; i<nTiles_;i++){
    long size = (long) parts_[i].size();
    for (int ip= 0; ip < size; ip++) {
      Particle *iter = &parts_[i][ip];
      if(iter->my_id < outputCount_){
        sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);  
                        if(debug>1)fprintf(stderr,"    track file name %s\n",fname);
        pout=fopen(fname,"a");
        //assert(pout != NULL);
        fprintf(pout,"%e %.15e %.15e %.15e %.15e %.15e %.15e\n",t,
          iter->x[0],iter->x[1],iter->x[2],
          iter->v[0],iter->v[1],iter->v[2]);
        fclose(pout);
      }
    }
  }
  if(debug)fprintf(stderr,"rank=%d: finish writing particle tracks!\n",rank_MPI);
}
#undef _USE_MATH_DEFINES
