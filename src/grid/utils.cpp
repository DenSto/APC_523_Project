#include <stdlib.h>
#include <stdio.h>
#include "assert.h"

double **new_contiguous_2dArray(int nx, int ny){
    double *newarr;
    double **ret;
    int i;

    newarr = (double*) malloc(nx*ny*sizeof(double));
    ret= (double**) malloc(nx*sizeof(double*));
    assert(newarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
        ret[i] = &(newarr[i*ny]);
    }
    return ret;
}

