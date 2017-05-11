#ifndef GLOBALS_HPP
#define GLOBALS_HPP


#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (y)) ? (x) : (y))

typedef struct vector_t{
	double x1, x2, x3;
}vector;

// Units, see units.pdf for detail
// unit of time 1->33.4ps
#define UNIT_TIME 1

// unit of rad frequency 1/2*pi*33.4ps, in unit of GHz
#define UNIT_FRAD 1

// unit for acceleration in Newton's equation
#define UNIT_ACC 1

// unit of charge and current density in StatC/cm^3,
// equals to charge density of 1 electron in 1cc cell
#define UNIT_RHOJ 1

// unit for electric field 1->299.79 KV/cm
#define UNIT_EFIELD 1

// 1eV electron thermal velocity, in unit of c
#define UNIT_VTH 1

// electron plasma frequency f_{pe}=\omega_{pe}/(2\pi)
// in unit of KHz, assuming density = 1 cc
#define UNIT_FPE 1

// electron gyro frequency f_{ce}=\Omega_{ce}/(2\pi)
// in unit of GHz, assuming B=1->1KG
#define UNIT_FCE 1

	
enum fieldID {E_X, E_Y,E_Z, B_X, B_Y, B_Z};

#ifdef MAIN_CPP
int rank_MPI, size_MPI;
int useVecDeposition;
int debug; // printf debug flag
double time_phys; // current physical time in simulation
double dt_phys; // current physical time step

#else // MAIN_CPP

extern int rank_MPI, size_MPI;
extern int useVecDeposition;
extern int debug;
extern double time_phys, dt_phys;

#endif // MAIN_CPP
#endif // GLOBALS_HPP
