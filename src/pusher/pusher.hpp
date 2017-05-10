#ifndef PUSHER_HPP
#define PUSHER_HPP

#include "../particles/particle.hpp"

class Pusher {
    public:
        virtual ~Pusher() {};
        virtual void Step(double *x, double *v, double q_m, Field_part field, double dt, double* ret) = 0;
};



#endif
