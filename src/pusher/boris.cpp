#include "boris.hpp"
#include "../globals.hpp"

Boris::Boris(){

}

Boris::~Boris(){
}

void Boris::Step(double* x, double* v, double q_m, Field_part field, double dt, double* ret){
  double v_mx,v_my,v_mz,v_px,v_py,v_pz;
  double v_ix,v_iy,v_iz;
  double bx,by,bz,ex,ey,ez,bsqr;
  double tx,ty,tz,sx,sy,sz,tsqr;
  double x_new, y_new,z_new, vx_new, vy_new,vz_new;
  double m,q,q_p;



  q_p = (0.5 * dt * q_m);
  q_p *= UNIT_ACC; // multiply unit of acceleration

  ex = field.e1;
  ey = field.e2;
  ez = field.e3;

  bx = field.b1;
  by = field.b2;
  bz = field.b3;

  tx = q_p*bx;
  ty = q_p*by;
  tz = q_p*bz;

  bsqr = bx*bx + by*by + bz*bz;
  tsqr = 2.0/(1.0 + q_p*q_p*bsqr);

  sx = tsqr*tx;
  sy = tsqr*ty;
  sz = tsqr*tz;

  v_mx = v[0] + q_p*ex; 
  v_my = v[1] + q_p*ey; 
  v_mz = v[2] + q_p*ez; 

  v_ix = v_mx + v_my*tz - v_mz*ty; 
  v_iy = v_my + v_mz*tx - v_mx*tz; 
  v_iz = v_mz + v_mx*ty - v_my*tx; 

  v_px = v_mx + v_iy*sz - v_iz*sy; 
  v_py = v_my + v_iz*sx - v_ix*sz; 
  v_pz = v_mz + v_ix*sy - v_iy*sx; 

  vx_new = v_px + q_p*ex; 
  vy_new = v_py + q_p*ey; 
  vz_new = v_pz + q_p*ez; 

  x_new = x[0] + dt*vx_new;
  y_new = x[1] + dt*vy_new;
  z_new = x[2] + dt*vz_new;

  /*
  part->x[0] = x_new;
  part->x[1] = y_new;
  part->x[2] = z_new;

  part->v[0] = vx_new;
  part->v[1] = vy_new;
  part->v[2] = vz_new;
*/
  ret[0] = x_new;
  ret[1] = y_new;
  ret[2] = z_new;

  ret[3] = vx_new;
  ret[4] = vy_new;
  ret[5] = vz_new;
}

