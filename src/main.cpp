//! Driver program for EMOOPIC
/*! ****************************************************** \n
 * To run the program use the following syntex             \n
 *   Serial version:                                       \n
 *       ./EMOOPIC <inputfile>                             \n
 *   MPI version                                           \n
 *       mpirun -np <nproc> ./EMOOPIC <inputfile>          \n
 *                                                         \n
 * The inputs are specified in <inputfile>                 \n
 **********************************************************/
#define MAIN_CPP
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include <iostream>

#include <sys/time.h>

#include "./globals.hpp"
#include "./domain/domain.hpp"
#include "./IO/input.hpp"
#include "./grid/grid.hpp"


#include "./particles/particle.hpp"
#include "./particles/particle_handler.hpp"
#include "./particles/particle_utils.hpp"

#include "./boundaries/particle_bc_factory.hpp"

#include "./pusher/pusher.hpp"
#include "./pusher/boris.hpp"

#if USE_MPI
    #include "mpi.h"  
#else
    #include<time.h>
#endif


static void stopwatch(short lap, FILE* fp){
  static struct timeval tv_prev, tv_curr;
  static short init = 1;
  static double tot = 0.0;

  assert(fp != NULL);

  if(init){
    init = 0;
    gettimeofday(&tv_prev,NULL);
#ifdef PART_IN_CELL
    fprintf(fp,"#[1]Output  [2]Interp [3]Push   [4]Cell [5]Bound [6]Deposit [7]Total\n");
#else
    fprintf(fp,"#[1]Output  [2]Sort   [3]Interp [4]Push [5]Bound [6]Deposit [7]Total\n");
#endif
    return;
  }

  gettimeofday(&tv_curr,NULL);
  double time = (double)(tv_curr.tv_sec - tv_prev.tv_sec) +
      		       1.0e-6*(double)(tv_curr.tv_usec - tv_prev.tv_usec);
  tv_prev = tv_curr;
  tot += time;
  fprintf(fp,"%f  ", time);
  if(lap){ 
    fprintf(fp,"%f \n",tot);
    fflush(fp);
    tot = 0.0;
  }
  return;
}


int main(int argc, char *argv[]){

    int size,rank=0;
    int sort = 0;

    /* Initialize *****************************************/
#if USE_MPI
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    double begin = MPI_Wtime();
    if(rank==0){
        printf("This simulation uses MPI domain decomposition.\n");
    }
#else
    size=1;
    rank=0;
    clock_t begin=clock();
    printf("This simulation is serial.\n");
#endif
    size_MPI=size;
    rank_MPI=rank;

    
    /* Read and check command line input ******************/
#ifdef PART_IN_CELL
    if(argc<3 && rank==0){
      fprintf(stderr,"The correct usage is:\n");
      fprintf(stderr,"%s <inputfile> <0- standard depo 1-vec depo>\n",argv[0]);
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,1);
#else
#endif 
#else
    if(argc<4 && rank==0){
      fprintf(stderr,"The correct usage is:\n");
      fprintf(stderr,"%s <inputfile> <0- standard depo 1-vec depo> <1-quicksort 2-counting>\n",argv[0]);
#endif
      exit(1);
    }


    useVecDeposition=atoi(argv[2]);
    
#ifndef PART_IN_CELL
    sort=atoi(argv[3]);
#endif

    double CFL = 1.0;

#ifdef PART_IN_CELL
    if(argc==4) CFL = atof(argv[3]);
#else
    if(argc==5) CFL = atof(argv[4]);
#endif

    char fname[50];
#ifdef VECTORIZE
    const char* vec = "_vec";
#else
    const char* vec = "";
#endif

#ifdef PART_IN_CELL
    sprintf(fname,"PIC_dep%d%s.dat",useVecDeposition,vec);
#else
    sprintf(fname,"PL_sort%d_dep%d%s.dat",sort,useVecDeposition,vec);
#endif
    static FILE *timer = fopen(fname,"w");
    fprintf(timer,"# %s \n",fname);

    /* Read and broadcast input file **********************/
    Input *input =  new Input();
    // Master read input file 
    if(rank==0){
      printf("Master reading input file...\n");
      int err = input->readinfo(argv[1]);
      // Check input self-consistency and load physical units to input
      err += input->ProcessInfo(); // should not be commented out
      if(err!=0) {
        std::cerr << "Input Error. Terminating..." << std::endl;
#if USE_MPI
         MPI_Abort(MPI_COMM_WORLD,1);
#else
         exit(1);
#endif
      }
    }

    Input_Info_t *input_info = input->getinfo();
#if USE_MPI
    // Master broadcast input info
    if(rank==0)printf("Master broadcasting input infomation...\n");
    input->passinfo();
#endif
    int restart = input_info->restart;
    debug = input_info->debug; // global debug flag
    if(debug>1) checkinput(input_info);

    /***************************************************************************/
    /* Initial setup                                                           */
    /***************************************************************************/
    if(rank==0)printf("Initial set up...\n");
    // Domain decomposition
    Domain *domain = new Domain(input_info);
    if(debug>1) checkdomain(domain);


    // Initialize particles and pusher
    Particle_Handler *part_handler = new Particle_Handler(); 
    part_handler->setPusher(new Boris());

    // Set up particle boundary conditions
    BC_Particle** bc = Part_BC_Factory::getInstance().constructConditions(domain,input_info);
    part_handler->setParticleBoundaries(bc);
    if(debug) fprintf(stderr,"rank=%d:Finish assigning particle boundary condition\n",rank);

    // Initialize grid
    Grid *grid;
    grid = new Grid(domain->getnxyz(),1,domain->getxyz0(),domain->getLxyz()); 

    // Load particles, allow restart
    if(rank==0)printf("    Loading particles...\n");
    part_handler->Load(input_info,domain,grid);
	  Particle_Compare* compare = new Particle_Compare(grid);
     if(debug) fprintf(stderr,"rank=%d: Finish loading particles\n",rank);   

    // Deposit charge and current from particles to grid
    if(restart==0 && strcmp(input_info->fields_init,"poisson")==0){
        part_handler->depositRhoJ(grid);
#if USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // Initialize fields by solving poisson or read files
    if(rank==0)printf("    Initializing fields...\n");
    grid->initializeFields(input_info);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing fields\n",rank);

    // Interpolate fields from grid to particle
    // Prepare initial push of particles
    part_handler->InterpolateEB(grid);
    if(debug) fprintf(stderr,"rank=%d: Finish initializing interpolation\n",rank);

    // particle output
    //part_handler->outputParticles(dir.c_str(),0,input_info); 
    if(debug) fprintf(stderr,"rank=%d: Finish writing initial particle output\n",rank);


    // prepare time step
    int nt = input_info->nt; //number of steps to run
    time_phys = input_info->t0; //initial time
    dt_phys = 0.1;
    dt_phys = CFL*part_handler->computeCFLTimestep(grid);
    if(debug) fprintf(stderr,"rank=%d: Finish preparing time step\n",rank);


    struct timeval tv_curr, tv_prev;
    double step_time;	

    gettimeofday(&tv_prev,NULL);

    std::string inputname = argv[1];
    std::string dir = inputname.substr(0,inputname.find_last_of('/')+1);

    stopwatch(0,timer);

    /***************************************************************************/
    /* Advance time step                                                       */
    /***************************************************************************/
    if(rank==0)printf("Advancing time steps with dt = %f ps\n",dt_phys*UNIT_TIME);
    for(int ti=0;ti<nt;ti++){

		gettimeofday(&tv_curr,NULL);
		step_time = (double)(tv_curr.tv_sec - tv_prev.tv_sec) +
				1.0e-6*(double)(tv_curr.tv_usec - tv_prev.tv_usec);

		tv_prev = tv_curr;
       if(rank==0 && ti%1==0)fprintf(stderr,"ti=%d\t\t t=%f ps steptime=%e s\n",ti,time_phys*UNIT_TIME,step_time);   


#ifndef PART_IN_CELL
    stopwatch(0,timer);
       if(input_info->nstep_sort > 0 && ti % input_info->nstep_sort == 0){
         switch(sort){
          case 0 : break;
          case 1 : 
            printf("Quicksorting particles\n");
            part_handler->SortParticles(*compare);
            break;
          case 2:
            printf("Countsorting particles\n");
            part_handler->CountingSortParticles(grid);
            break;

         }
       }
#endif
		
    stopwatch(0,timer);
       // Interpolate fields from grid to particle
       part_handler->InterpolateEB(grid);

    stopwatch(0,timer);
       /* push particles ***********************/
       part_handler->Push(dt_phys);

       // Pass particle across MPI boundaries, or implement physical boundary conditions
       // All particles are in physical cells, no particle lives in ghost cell
#ifdef PART_IN_CELL
    stopwatch(0,timer);
       part_handler->updateCellLists(grid,0);
#endif
    stopwatch(0,timer);
       part_handler->executeParticleBoundaryConditions(grid);

       // remove any particles left in the ghost cells
       part_handler->clearGhosts();

       // only deposit particles in physical cells
    stopwatch(0,timer);
       part_handler->depositRhoJ(grid);


    stopwatch(1,timer);
       time_phys += dt_phys;

       part_handler->outputParticles(dir.c_str(),ti+1,input_info); 

     } // timestep loop

     if(rank==0) fprintf(stderr,"ti=%d\t\t t=%f ps\n",nt,time_phys*UNIT_TIME);   
     if(rank==0) printf("***Timestep loop complete***\n");

    /***************************************************************************/
    /* output, finalize                                                        */
    /***************************************************************************/
    if(rank==0)printf("Writing final output files...\n");
    //writeoutput(grid,part_handler); //MPI
    if(debug) fprintf(stderr,"rank=%d: Finish final output\n",rank);

    // free memory
    delete domain;
    delete [] bc; // particle boundary condition
    delete part_handler;
	delete compare;
    delete grid;
    delete input;
    if(debug) fprintf(stderr,"rank=%d: Finish free\n",rank);

#if USE_MPI
    double time = MPI_Wtime()-begin;
    double maxtime;
    MPI_Reduce(&time,&maxtime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#else
    clock_t end = clock();
    double maxtime = (double)(end-begin)/CLOCKS_PER_SEC;
#endif

    if(rank==0){
        printf("Program completed successfully!\n");
        printf("Elapsed: %f seconds\n",maxtime);
    }
#if USE_MPI
	MPI_Finalize();
#endif
    fclose(timer);

}
