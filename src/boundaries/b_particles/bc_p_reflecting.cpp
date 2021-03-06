#include "../../globals.hpp"
#include "../particles_boundary.hpp"
#include "../particle_bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_Reflecting : public BC_Particle {
	public:
		BC_P_Reflecting(Domain* domain, int dim_Index, short isRight, std::string type);
		~BC_P_Reflecting();
		int computeParticleBCs(std::vector<Particle> *pl);
		int completeBC(std::vector<Particle> *pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;
};

BC_P_Reflecting::BC_P_Reflecting(Domain* domain, int dim_Index, short isRight, std::string type) 
	:	dim_index_(dim_Index),
		isRight_(isRight),
		type_(type)
{
	xMin_ = domain->getxyz0()[dim_index_];
        xMax_ = xMin_+domain->getLxyz()[dim_index_];
	if(debug>1)fprintf(stderr,"rank=%d:dim=%d,isRight=%d,reflect_BC,xMin=%f,xMax=%f\n",
                                   rank_MPI,dim_index_,isRight_,xMin_,xMax_); 	

}

BC_P_Reflecting::~BC_P_Reflecting(){
}

int BC_P_Reflecting::completeBC(std::vector<Particle> *pl){
	// Nothing to do
	return 0;
}


inline int BC_P_Reflecting::particle_BC(Particle* p){
	if(p->x[dim_index_] > xMax_ && isRight_){
		p->x[dim_index_] = 2.0*xMax_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	} 
	if(p->x[dim_index_] < xMin_ && !isRight_){
		p->x[dim_index_] = 2.0*xMin_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	}
	return 0;
}

// cycle through all particles
int BC_P_Reflecting::computeParticleBCs(std::vector<Particle> *pl) {
  long size = pl->size();
#pragma simd
  for(long i = 0; i < size; i++){
    (*pl)[i].isGhost = (*pl)[i].isGhost || 
      particle_BC(&(*pl)[i]);
  }
  return completeBC(pl);
}

// Registers bounary condition into BC_Factory dictionary
static RegisterParticleBoundary instance("reflecting", makeBCParticle<BC_P_Reflecting>);
