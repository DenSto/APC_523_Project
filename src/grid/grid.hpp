/*
 * Grid : Structure that represents the grid on which fields live.
 *         Contains grid spacing, lenghts and information on the ghost cells. 
 */
#ifndef GRID_HPP
#define GRID_HPP

#include "../globals.hpp"
#include "../IO/input.hpp"

typedef struct Couple {
  double Ex[3][3][3];
  double Ey[3][3][3];
  double Ez[3][3][3];

  double Bx[3][3][3];
  double By[3][3][3];
  double Bz[3][3][3];
} Couple;

typedef struct PartProp{
  double d;
  double jx;
  double jy;
  double jz;
}PartProp;

typedef struct PartCouple{
  PartProp pp[3][3][3];
} PartCouple;

inline int intFloor(double x){
  return (int)(x + 10000000) - 10000000;
}

class Grid {

public:
  Grid(int *nxyz, int nGhosts, double *xyz0, double *Lxyz);
  virtual ~Grid(); 

  inline int getCellID(double x, double y, double z){
    int ix = intFloor((x - x0_) * idx_);
    int iy = intFloor((y - y0_) * idy_);
    int iz = intFloor((z - z0_) * idz_);
    return (ny_*nz_)*ix + nz_*iy + iz;
  }

  inline int getCellID_wGhost(double x, double y, double z){
    int ix = intFloor((x - x0_) * idx_);
    if(ix < 0 || ix >= nx_) return -1;

    int iy = intFloor((y - y0_) * idy_);
    if(iy < 0 || iy >= ny_) return -2;

    int iz = intFloor((z - z0_) * idz_);
    if(iz < 0 || iz >= nz_) return -3;

    return (ny_*nz_)*ix + nz_*iy + iz;
  }

  void getCouple(int i, int j, int k, Couple* cp);
  void getCoupleID(long id,  Couple* cp);


  void setPartProp(int i, int j, int k, PartProp pp);
  void addPartProp(int i, int j, int k, PartProp pp);

  void getnxyz(int* nxyz);
  void getL0(double* L0);
  void getSpacing(double* space);
  void getInvSpacing(double* space);

  void constE(double x, double y, double z);
  void constB(double x, double y, double z);

  int  getGhost() {return nGhosts_;};

  double getx0() {return x0_;};
  double gety0() {return y0_;};
  double getz0() {return z0_;};

  double getidx() {return idx_;};
  double getidy() {return idy_;};
  double getidz() {return idz_;};
  double ***getRho() {return rho_;};

  void initializeFields(Input_Info_t* input);

#ifdef PARTICLE_GRID
  void pushParticle(Particle part);
  void updateParticleLocation();
#endif

protected:
  double ***Ex_;
  double ***Ey_;
  double ***Ez_;
  double ***Bx_;
  double ***By_;
  double ***Bz_;
  
  double ***rho_;

  PartProp ***pp_;

  int nx_;
  int ny_;
  int nz_;
  int nGhosts_;

  double x0_;
  double y0_;
  double z0_;
  double Lx_;
  double Ly_;
  double Lz_;

  double dx_;
  double dy_;
  double dz_;
  double idx_;
  double idy_;
  double idz_;

#ifdef PARTICLE_GRID
  std::vector<Particle>* parts_;
#endif

};

#endif
