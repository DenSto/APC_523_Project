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

  int* n= domain->getnxyz();
  ncell_ = n[0]*n[1]*n[2];

#ifdef PART_IN_CELL
  parts_ = new std::vector<Particle>[ncell_+1];
#else
  pa_= (int*) malloc(sizeof(int) *(ncell_+1));
  pa_save_ = (int*) malloc(sizeof(int) *( ncell_ +1));
#endif

  int nspec      = input_info->nspecies;
  int npart      = input_info->nparticles_tot*ncell_;
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

#ifdef PART_IN_CELL
      int id = grid->getCellID(p.x[0],p.x[1],p.x[2]);
      parts_[id].push_back(p);
#else
      parts_.push_back(p);
#endif
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
#ifdef PART_IN_CELL
    int id = grid->getCellID(p.x[0],p.x[1],p.x[2]);
    parts_[id].push_back(p);
#else
    parts_.push_back(p);
#endif
    np_++;
  }
}

void Particle_Handler::Push(double dt){
#ifdef PART_IN_CELL
  for(long i = 0; i < ncell_; i++){
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
#else
  for(long ip=0;ip<np_;ip++){
    double qom = parts_[ip].q/parts_[ip].m;
    double ret[6];
    pusher_->Step(parts_[ip].x,parts_[ip].v,qom, parts_[ip].field,dt,ret);

    parts_[ip].x[0] = ret[0];
    parts_[ip].x[1] = ret[1];
    parts_[ip].x[2] = ret[2];
    parts_[ip].v[0] = ret[3];
    parts_[ip].v[1] = ret[4];
    parts_[ip].v[2] = ret[5];
  }
#endif
}

long Particle_Handler::nParticles(){
#ifdef PART_IN_CELL
  return -1;
#else
  return np_;
#endif
}

void Particle_Handler::incrementNParticles(int inc){
#ifndef PART_IN_CELL
  np_+=inc;
#endif
}

void Particle_Handler::InterpolateEB(Grid* grid){
  int nxyz[3];
  double pos[3]; //Vector of position of particle.
  double icell[3]; //Vector of inverse lengths of unit cell.
  double L0[3]; //Vector of lengths of unit cell.
  double weight[3][3][3];

  int pID, cellID = -1;

  Couple cp;

  int is,js,ks;
  int ic,jc,kc;

  //Get lengths of grid cells.
  grid->getnxyz(nxyz);
  grid->getInvSpacing(icell);
  grid->getL0(L0);

  //int count =0;

#ifdef PART_IN_CELL
  for (long j = 0; j < ncell_; j++){
    grid->getCoupleID(j,&cp);
    long size = (long) parts_[j].size();
#pragma vector
    for (long i=0; i < size; i++) {
      //Get position of particle.
      pos[0] = parts_[j][i].x[0];
      pos[1] = parts_[j][i].x[1];
      pos[2] = parts_[j][i].x[2];

      getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);
      interpolateFields(&cp,weight,&parts_[j][i].field);
    }
  }
#else
  for (long i=0; i<np_; i++) {
    assert(!parts_[i].isGhost);
  
    //Get position of particle.
    pos[0] = parts_[i].x[0];
    pos[1] = parts_[i].x[1];
    pos[2] = parts_[i].x[2];

    getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);

    pID = nxyz[1]*nxyz[2]*ic + nxyz[2]*jc + kc;
    if(pID != cellID){
      cellID=pID;
      grid->getCouple(is,js,ks,&cp);
    }
    interpolateFields(&cp,weight,&parts_[i].field);
  }
#endif
}

void Particle_Handler::depositRhoJ(Grid *grid){
  int nxyz[3];
  double pos[3]; //Vector of position of particle.
  double icell[3]; //Vector of inverse lengths of unit cell.
  double L0[3]; //Vector of lengths of unit cell.
  double weight[3][3][3];

  if(useVecDeposition){
#ifdef PART_IN_CELL
    depositRho(grid,parts_);
#else
    depositRho(grid,&parts_);
#endif
    return;
  }

  int pID, cellID = -1;

  double ***rho = grid->getRho();
  int nghost = grid->getGhost();

  Couple cp;


  //Get lengths of grid cells.
  grid->getInvSpacing(icell);
  grid->getL0(L0);

  PartProp cur;

#ifdef PART_IN_CELL
  for (long n=0; n < ncell_ ; n++) {
    long size = (long) parts_[n].size();
    for (long ip=0; ip< size; ip++) {
      int is,js,ks;
      int ic,jc,kc;
      Particle *p = &parts_[n][ip];
  
      //Get position of particle.
      pos[0] = p->x[0];
      pos[1] = p->x[1];
      pos[2] = p->x[2];

      getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);

      double q = p->q;
      for( int i = 0; i < 3; i++){
        for( int j = 0; j < 3; j++){
#pragma simd
          for( int k = 0; k < 3; k++){
            rho[is +nghost + i][js + nghost + j][ks + nghost + k] += q*weight[i][j][k];	    
  //        cur.d = p->m*weight[i][j][k];
  //        cur.jx = p->q*p->v[0]*weight[i][j][k];
  //        cur.jy = p->q*p->v[1]*weight[i][j][k];
  //        cur.jz = p->q*p->v[2]*weight[i][j][k];
  //        grid->addPartProp(is + i, js + j, ks + k, cur);
          }
        }
      }
    }
  }
#else
  for (long ip=0; ip<np_; ip++) {
    int is,js,ks;
    int ic,jc,kc;
    assert(!parts_[ip].isGhost);
  
    //Get position of particle.
    pos[0] = parts_[ip].x[0];
    pos[1] = parts_[ip].x[1];
    pos[2] = parts_[ip].x[2];

    getwei_TSC(L0,icell,pos[0],pos[1],pos[2],weight,&is,&js,&ks,&ic,&jc,&kc);

    for( int i = 0; i < 3; i++){
      for( int j = 0; j < 3; j++){
        for( int k = 0; k < 3; k++){
          cur.d = parts_[ip].m*weight[i][j][k];
          cur.jx = parts_[ip].q*parts_[ip].v[0]*weight[i][j][k];
          cur.jy = parts_[ip].q*parts_[ip].v[1]*weight[i][j][k];
          cur.jz = parts_[ip].q*parts_[ip].v[2]*weight[i][j][k];
          grid->addPartProp(is + i, js + j, ks + k, cur);
        }
      }
    }
  }
#endif
}

//! Sort particles based on grid location. 
/*
 * Sorts particles based on grid location using std::sort which is O(n log n).
 * This should be called often to ensure cache hits.  Should be emperically determined.
 */
void Particle_Handler::SortParticles(Particle_Compare comp){
#ifndef PART_IN_CELL
  if(debug>0 && rank_MPI == 0) fprintf(stderr,"Sorting Particles.\n");
  std::sort(parts_.begin(),parts_.end(),comp);
#endif
}

void Particle_Handler::CountingSortParticles(Grid* grid){
#ifndef PART_IN_CELL
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
#endif
}


double Particle_Handler::computeCFLTimestep(Grid* grid){
  double maxV[3], dx[3];

  double mindt = 1e20;

  double v;
  grid->getSpacing(dx);
#ifdef PART_IN_CELL
  for (int i= 0; i<ncell_; i++) {
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
#else
  for (int i= 0; i<np_; i++) {
    for(int j = 0; j < 3; j++){
      v = parts_[i].v[j];
      if(v*v > 0){
        if(fabs(dx[j]/v) < mindt) mindt = fabs(dx[j]/v);
      } 
    }
  }
#endif
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
#ifdef PART_IN_CELL
  long size = 0;
  for(int i = 0; i < ncell_; i++) size += (long) parts_[i].size();
  parts_[ncell_].clear();
#else
  int nghost = 0;
  for(std::vector<Particle>::iterator iter = parts_.begin(); iter != parts_.end();){
    if(iter->isGhost){
      std::swap(*iter, parts_.back());
      parts_.pop_back();
      nghost +=1;
    } else {
      iter++;
    }
  }
  assert((long)parts_.size() == np_);
        if(debug>2)fprintf(stderr,"rank=%d: %d ghosts are cleared.\n",rank_MPI,nghost);
#endif
}

#ifdef PART_IN_CELL
void Particle_Handler::updateCellLists(Grid *grid,short justGhosts){
  long i = justGhosts ? ncell_ : 0;
  for(; i < ncell_ + 1; i++){
    for(std::vector<Particle>::iterator iter = parts_[i].begin(); iter != parts_[i].end();){
      long id = grid->getCellID_wGhost(iter->x[0],iter->x[1],iter->x[2]);
      id = id < 0 ? ncell_ : id;
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
  long ind = id < 0 ? ncell_ : id;
  parts_[ind].push_back(p);
}
#endif


void Particle_Handler::executeParticleBoundaryConditions(Grid* grid){
  for(int i = 0; i < 6; i++){
    // determine whether particles are ghost
        // place ghost particles
        // change the number of particles np_ in each domain
#ifdef PART_IN_CELL
    int inc = boundaries_[i]->computeParticleBCs(&parts_[ncell_]);
    updateCellLists(grid,1);
#else
    int inc = boundaries_[i]->computeParticleBCs(&parts_);
    incrementNParticles(inc);
#endif
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
#ifdef PART_IN_CELL
    for(int i = 0; i < ncell_; i++){
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
#else
    for(std::vector<Particle>::iterator iter = parts_.begin();iter!=parts_.end();++iter){
     if(iter->my_id < outputCount_){
        char fname[100];
        sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);  
        FILE *pout=fopen(fname,"w");
        fprintf(pout,"[1] time [2] x [3] y [4] z  [5] vx   [6] vy   [7] vz\n");
        fclose(pout);
      }
    }
#endif
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
#ifdef PART_IN_CELL
  for(int i =0; i<ncell_;i++){
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
#else
  for(std::vector<Particle>::iterator iter = parts_.begin();iter!=parts_.end();++iter){
    if(iter->my_id < outputCount_){
      sprintf(fname,"%strack_%d_%ld.dat",basename,iter->initRank,iter->my_id);  
                        if(debug>1)fprintf(stderr,"    track file name %s\n",fname);
      pout=fopen(fname,"a");
      assert(pout != NULL);
      fprintf(pout,"%e %.15e %.15e %.15e %.15e %.15e %.15e\n",t,
          iter->x[0],iter->x[1],iter->x[2],
          iter->v[0],iter->v[1],iter->v[2]);
      fclose(pout);
    }
  }
#endif
  if(debug)fprintf(stderr,"rank=%d: finish writing particle tracks!\n",rank_MPI);
}
#undef _USE_MATH_DEFINES
