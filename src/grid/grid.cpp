#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "grid.hpp"
#include "utils.hpp"

Grid::Grid(int *nxyz, int nGhosts, double *xyz0, double *Lxyz):
    nx_(nxyz[0]),
    ny_(nxyz[1]),
    nz_(nxyz[2]),
    nGhosts_(nGhosts),
    x0_(xyz0[0]),
    y0_(xyz0[1]),
    z0_(xyz0[2]),
    Lx_(Lxyz[0]),
    Ly_(Lxyz[1]),
    Lz_(Lxyz[2]),
    dx_(Lxyz[0]/nx_),
    dy_(Lxyz[1]/ny_),
    dz_(Lxyz[2]/nz_),
    idx_(1.0/dx_),
    idy_(1.0/dy_),
    idz_(1.0/dz_)
{
  Ex_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  Ey_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  Ez_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  Bx_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  By_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  Bz_ =new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  rho_=new_contiguous_3dArray<double>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
  pp_=new_contiguous_3dArray<PartProp>(nx_ + 2*nGhosts_, ny_ + 2*nGhosts_, nz_ + 2*nGhosts_);    
}

Grid::~Grid(){

}

void Grid::getSpacing(double* space){
  space[0] = dx_;  
  space[1] = dy_;  
  space[2] = dz_;  
}

void Grid::getInvSpacing(double* space){
  space[0] = dx_ == 0.0 ? 0 : 1.0/dx_;  
  space[1] = dy_ == 0.0 ? 0 : 1.0/dy_;  
  space[2] = dz_ == 0.0 ? 0 : 1.0/dz_;  
}

void Grid::getL0(double* L0){
  L0[0] = x0_;  
  L0[1] = y0_;  
  L0[2] = z0_;  
}

void Grid::getnxyz(int* nxyz){
  nxyz[0] = nx_;
  nxyz[1] = ny_;
  nxyz[2] = nz_;
}

void Grid::constE(double x, double y, double z){
  int i,j,k;  
  for(i = 0; i < nx_ + 2*nGhosts_; i++){
    for(j = 0; j < ny_ + 2*nGhosts_; j++){
      for(k = 0; k < nz_ + 2*nGhosts_; k++){
        Ex_[i][j][k] = x;
        Ey_[i][j][k] = y;
        Ez_[i][j][k] = z;
      }
    }
  }
}

void Grid::constB(double x, double y, double z){
  int i,j,k;  
  for(i = 0; i < nx_ + 2*nGhosts_; i++){
    for(j = 0; j < ny_ + 2*nGhosts_; j++){
      for(k = 0; k < nz_ + 2*nGhosts_; k++){
        Bx_[i][j][k] = x;
        By_[i][j][k] = y;
        Bz_[i][j][k] = z;
      }
    }
  }
}

void Grid::initializeFields(Input_Info_t *input_info){
  int restart = input_info->restart;
  if(restart==0 && strcmp(input_info->fields_init,"constant")==0){
    if(rank_MPI==0)printf("        Initializing fields to constants in input file...\n");
    double *E0 = input_info->E0;
    double *B0 = input_info->B0;
    constE(E0[0],E0[1],E0[2]); 
    constB(B0[0],B0[1],B0[2]); 
  }
}

void Grid::getCoupleID(long id, Couple* cp){
  int iz = id % nz_;
  int iy = ((id - iz) / nz_) % ny_;
  int ix = ((( id - iz)/ nz_ ) - iy) / ny_;

  getCouple(ix-1,iy-1,iz-1,cp);
}
void Grid::getCouple(int is, int js, int ks, Couple* cp){
  int i,j,k;
  int ic = is + nGhosts_;
  int jc = js + nGhosts_;
  int kc = ks + nGhosts_;

  for(i = 0; i < 3; i++){
    for(j = 0; j < 3; j++){
      for(k = 0; k < 3; k++){
        cp->Ex[i][j][k] = Ex_[ic+i][jc+j][kc+k];
        cp->Ey[i][j][k] = Ey_[ic+i][jc+j][kc+k];
        cp->Ez[i][j][k] = Ez_[ic+i][jc+j][kc+k];

        cp->Bx[i][j][k] = Bx_[ic+i][jc+j][kc+k];
        cp->By[i][j][k] = By_[ic+i][jc+j][kc+k];
        cp->Bz[i][j][k] = Bz_[ic+i][jc+j][kc+k];
      }
    }
  }
}

void Grid::setPartProp(int i, int j, int k, PartProp pp){
  pp_[i + nGhosts_][j + nGhosts_][k + nGhosts_] = pp;
}

void Grid::addPartProp(int i, int j, int k, PartProp pp){
  int ic = i + nGhosts_;
  int jc = j + nGhosts_;
  int kc = k + nGhosts_;

  pp_[ic][jc][kc].d  += pp.d;
  pp_[ic][jc][kc].jx += pp.jx;
  pp_[ic][jc][kc].jy += pp.jy;
  pp_[ic][jc][kc].jz += pp.jz;
}
