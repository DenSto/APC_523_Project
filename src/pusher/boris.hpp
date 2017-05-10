#ifndef BORIS_HPP
#define BORIS_HPP
#include "pusher.hpp"

class Boris : public Pusher {
    public:
        Boris();
        ~Boris();
        void Step(double *x, double *v, double q_m, Field_part field, double dt, double* ret);
		
};

#endif
