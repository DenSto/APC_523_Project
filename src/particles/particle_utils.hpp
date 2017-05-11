#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include <stdlib.h>
#include "particle.hpp"
#include "../grid/grid.hpp"
#include "../globals.hpp"

#define LVEC 8

/*! Function to use are a comparator in std::sort for std::vec<Particle>
* Idea: Particle is sorted by outer index first 
* (i.e. particle at [2][12][43] should be closer to beginning of array than
* particle at [24][12][43])
*
* At the moment, implemented very slowly!!! Should be modified two ways:
*
*   1) Instead of comparing cell ID, compare i,j,k indice locations individually
*      to save time
*   2) Bring the code to calculating i,j,k into the comparison routine
*/
class Particle_Compare {
  public:
  Particle_Compare(Grid* grid)  
  {
    grid->getInvSpacing(idx_);
    grid->getL0(x0_);
  }

  bool operator()(Particle const a, Particle const b) const {
// Sort by individual i,j,k
    int na,nb;  

    na = (int)((a.x[0] - x0_[0])*idx_[0]);
    nb = (int)((b.x[0] - x0_[0])*idx_[0]);
    if(na < nb) return 0;
    if(na > nb) return 1;

    na = (int)((a.x[1] - x0_[1])*idx_[1]);
    nb = (int)((b.x[1] - x0_[1])*idx_[1]);
    if(na < nb) return 0;
    if(na > nb) return 1;

    na = (int)((a.x[2] - x0_[2])*idx_[2]);
    nb = (int)((b.x[2] - x0_[2])*idx_[2]);
    if(na <= nb) return 0;
    if(na > nb)  return 1;

    return 0;
  }

  private:
  
    double idx_[3];
    double x0_[3];
};


inline int cell(const double L0, const double x, const double idx, int *i, double*a){
   *a = (x - L0)*idx;
   *i = (int)(*a);
   if (((*a) -(*i)) < 0.5) return 0; 
   else return 1;
} 

inline void getwei_TSC(double L0[3], double idx[3], double x1, double x2, double x3,
                       double weight[3][3][3], int *is, int *js, int *ks,
                       int *ic, int *jc, int *kc)
{
  int i, j, k;
  double a, b, c, d;  /* grid coordinate for the position (x1,x2,x3) */
  double wei1[3], wei2[3], wei3[3];/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (idx[0] > 0.0) {
    cell(L0[0], x1, idx[0], &i, &a);    /* x1 index */
    *ic = i;
    *is = i - 1;        /* starting x1 index, wei[0] */
    d = a - i;
    wei1[0] = 0.5*SQR(1.0-d);     /* 0: left; 2: right */
    wei1[1] = 0.75-SQR(d-0.5);      /* one direction weight */
    wei1[2] = 0.5*SQR(d);
  }
  else { /* x1 dimension collapses */
    *ic = 0;
    *is = 0;
    wei1[1] = 0.0;  wei1[2] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (idx[1] > 0.0) {
    cell(L0[1], x2, idx[1], &j, &b);    /* x2 index */
    *jc = j;
    *js = j - 1;        /* starting x2 index */
    d = b - j;
    wei2[0] = 0.5*SQR(1.0-d);     /* 0: left; 2: right */
    wei2[1] = 0.75-SQR(d-0.5);      /* one direction weight */
    wei2[2] = 0.5*SQR(d);
  }
  else { /* x2 dimension collapses */
    *jc = 0;
    *js = 0;
    wei2[1] = 0.0;  wei2[2] = 0.0;
    wei2[0] = 1.0;
  }

/* x3 direction */
  if (idx[2] > 0.0) {
    cell(L0[2], x3, idx[2], &k, &c);    /* x3 index */
    *kc = k;
    *ks = k - 1;        /* starting x3 index */
    d = c - k;
    wei3[0] = 0.5*SQR(1.0-d);     /* 0: left; 2: right */
    wei3[1] = 0.75-SQR(d-0.5);      /* one direction weight */
    wei3[2] = 0.5*SQR(d);
  }
  else { /* x3 dimension collapses */
    *kc = 0;
    *ks = 0;
    wei3[1] = 0.0;  wei3[2] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


inline void interpolateFields(Couple *cp, double weight[3][3][3],Field_part *field){
  int i,j,k;

  field->e1=0.0;
  field->e2=0.0;
  field->e3=0.0;

  field->b1=0.0;
  field->b2=0.0;
  field->b3=0.0;

  double totalwei =0.0;
  
  for(i =0; i < 3; i++){
    for(j =0; j < 3; j++){
      for(k =0; k < 3; k++){
        totalwei+=weight[i][j][k];
        field->e1+=weight[i][j][k]*cp->Ex[i][j][k];
        field->e2+=weight[i][j][k]*cp->Ey[i][j][k];
        field->e3+=weight[i][j][k]*cp->Ez[i][j][k];

        field->b1+=weight[i][j][k]*cp->Bx[i][j][k];
        field->b2+=weight[i][j][k]*cp->By[i][j][k];
        field->b3+=weight[i][j][k]*cp->Bz[i][j][k];
      }
    }
  }
  double iw = 1.0/totalwei;
  field->e1 *= iw;
  field->e2 *= iw;
  field->e3 *= iw;

  field->b1 *= iw;
  field->b2 *= iw;
  field->b3 *= iw;
}



inline void depositRho(Grid* grid, std::vector<Particle>* nparts){

  static double* rhocells;

  int nghost = grid->getGhost();
  double idx[3];  grid->getInvSpacing(idx);
  int nxyz[3]; grid->getnxyz(nxyz);
  double L0[3]; grid->getL0(L0);

  double*** rho = grid->getRho();

  double ww, wwx, wwy,wwz;
  double x,y,z,invvol,q, wq0, wq, sxy, syy0, syy1, syy2, sxx0, sxx1, sxx2;
  double xint,yint,zint, xintsq,yintsq,zintsq;
  double sz0[LVEC], sz1[LVEC], sz2[LVEC];
  int ICELL[LVEC];

  double ww0[8*LVEC], www[LVEC*8];

  int j,k,l,ix,iy,iz,ic,n,nn,nv;
  int ncz, ncyz;
  long ip;


  invvol = idx[0]*idx[1]*idx[2];

  // need at least one ghost zone for scheme to make sense (so one ghost on each side)
  int ncells = nxyz[0]*nxyz[1]*nxyz[2];
  int ncellsTot = (nxyz[0]+2)*nxyz[1]*nxyz[2];

  if(rhocells == NULL) rhocells = (double*) malloc(sizeof(double)*ncellsTot*8);
  for(int i = 0; i < ncellsTot*8; i++) rhocells[i] = 0.0;

  ncz = nxyz[2] + 2;
  ncyz = ncz*nxyz[1];

  
#ifndef PART_IN_CELL
  ncells = 1;
#endif

  for(int i = 0; i < LVEC*8; i++) ww0[i] = 0.0;
    
  for(int iter = 0; iter< ncells; iter++ ){
    long np = nparts[iter].size();

    for(ip = 0; ip < np; ip += LVEC){
      for(n = 0; n < MIN(LVEC,np-ip); n++){
        nn= ip + n;
        //wq0 = p.q*invvol;
        x = (nparts[iter][nn].x[0]-L0[0])*idx[0];
        y = (nparts[iter][nn].x[1]-L0[0])*idx[1];
        z = (nparts[iter][nn].x[2]-L0[2])*idx[2];

        j = (int)(x + 100000.0) - 100000;
        k = (int)(y + 100000.0) - 100000;
        l = (int)(z + 100000.0) - 100000;
        
        ICELL[n]= (1+l) + k*ncz + j*ncyz;

        xint = x-j;
        yint = y-k;
        zint = z-l;
        xintsq=SQR(xint);
        yintsq=SQR(yint);
        zintsq=SQR(zint);

        syy0  =0.5*SQR(0.5 - yint);
        syy1  =(0.75-yintsq);
        syy2  =0.5*SQR(0.5+yint);
        sxx0  =0.5*SQR(0.5-xint);
        sxx1  =(0.75-xintsq);
        sxx2  =0.5*SQR(0.5+xint);
        sz0[n]=0.5*SQR(0.5-zint);
        sz1[n]=(0.75-zintsq);
        sz2[n]=0.5*SQR(0.5+zint);
        www[8*n + 0] = syy0*sxx0;
        www[8*n + 1] = syy1*sxx0;
        www[8*n + 2] = syy2*sxx0;
        www[8*n + 3] = syy0*sxx1;
        www[8*n + 4] = syy2*sxx1;
        www[8*n + 5] = syy0*sxx2;
        www[8*n + 6] = syy1*sxx2;
        www[8*n + 7] = syy2*sxx2;
        sxy=syy1*sxx1;
        ww0[8*n + 0]=sxy*sz0[n];
        ww0[8*n + 1]=sxy*sz1[n];
        ww0[8*n + 2]=sxy*sz2[n];
        for(nv = 0; nv < 8; nv++){
          wq0+= www[8*n + nv];
          wq0+= ww0[8*n + nv];
        }
	double fac = nparts[iter][nn].q / wq0;
        for(nv = 0; nv < 8; nv++){
          www[8*n + nv] *= fac;
          ww0[8*n + nv] *= fac;
        }
      }
      j += nghost; k += nghost; l += nghost;

      for(n = 0; n < MIN(LVEC,np-ip); n++){
#pragma simd
        for(nv = 0; nv < 8; nv++){
          ww=www[8*n + nv];
          rhocells[8*(ICELL[n]-1) + nv] += ww*sz0[n];
          rhocells[8*(ICELL[n]  ) + nv] += ww*sz1[n];
          rhocells[8*(ICELL[n]+1) + nv] += ww*sz2[n];
        }
        for(nv = -1; nv <= 1; nv++){
            rho[j][k][l+nv] += ww0[8*n + nv + 1];
        }
      }
    }
  }
  // REDUCTION OVER ENTIRE GRID
  for(ix = nghost; ix < nghost + nxyz[0] ; ix++){
    for(iy = nghost; iy < nghost + nxyz[1]; iy++){
#pragma simd
      for(iz = nghost - 1; iz < nghost + nxyz[2] + 1; iz++){
        ic=iz+(iy-nghost)*ncz+(ix-nghost)*ncyz;
        rho[ix-1][iy-1][iz] += rhocells[8*ic + 0];
        rho[ix-1][iy][iz]   += rhocells[8*ic + 1];
        rho[ix-1][iy+1][iz] += rhocells[8*ic + 2];
        rho[ix][iy-1][iz]   += rhocells[8*ic + 3];
        rho[ix][iy+1][iz]   += rhocells[8*ic + 4];
        rho[ix+1][iy-1][iz] += rhocells[8*ic + 5];
        rho[ix+1][iy][iz]   += rhocells[8*ic + 6];
        rho[ix+1][iy+1][iz] += rhocells[8*ic + 7];
      }
    }
  } 
}
#endif
