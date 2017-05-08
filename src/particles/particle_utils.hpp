#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include "particle.hpp"
#include "../grid/grid.hpp"

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
    idx_ = grid->getidx();
    idy_ = grid->getidy();
    idz_ = grid->getidz();

    x0_ = grid->getx0();
    y0_ = grid->gety0();
    z0_ = grid->getz0();
  }

  bool operator()(Particle const a, Particle const b) const {
// Sort by individual i,j,k
    int na,nb;  

    na = (int)((a.x[0] - x0_)*idx_);
    nb = (int)((b.x[0] - x0_)*idx_);
    if(na < nb) return 0;
    if(na > nb) return 1;

    na = (int)((a.x[1] - y0_)*idy_);
    nb = (int)((b.x[1] - y0_)*idy_);
    if(na < nb) return 0;
    if(na > nb) return 1;

    na = (int)((a.x[2] - z0_)*idz_);
    nb = (int)((b.x[2] - z0_)*idz_);
    if(na <= nb) return 0;
    if(na > nb)  return 1;

    return 0;

// Sort by cell ID (probably overkill)
/*
    int id1 = grid_->getCellID(a.x[0],a.x[1],a.x[2]);
    int id2 = grid_->getCellID(b.x[0],b.x[1],b.x[2]);

    // if all indices are equal, return 0;
    if(id1 <= id2)
      return 0;
    else
      return 1;
*/
  }

  private:
  
    double idx_,idy_,idz_;
    double x0_,y0_,z0_;
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


inline void interpolateFields(Couple cp, double weight[3][3][3],Field_part *field){
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
        field->e1+=weight[i][j][k]*cp.Ex[i][j][k];
        field->e2+=weight[i][j][k]*cp.Ey[i][j][k];
        field->e3+=weight[i][j][k]*cp.Ez[i][j][k];

        field->b1+=weight[i][j][k]*cp.Bx[i][j][k];
        field->b2+=weight[i][j][k]*cp.By[i][j][k];
        field->b3+=weight[i][j][k]*cp.Bz[i][j][k];
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
#endif
