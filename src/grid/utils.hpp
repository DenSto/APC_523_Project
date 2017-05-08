/*
 * Utils : Helpful utilities.
 *
 * 			Specifically: constructor of 2D array that's guaranteed to be contiguous
 * 			              in memory (unnecessary?).
 */
#ifndef UTILS_H
#define UTILS_H


double **new_contiguous_2dArray(int nx, int ny);
void free_contiguous_2dArray(double **array);

//double ***new_contiguous_3dArray(int nx, int ny, int nz);
//void free_contiguous_2dArray(double ***array);

/*
double ***new_contiguous_3dArray(int nx, int ny, int nz){
    double ***ret;
    double **nextarr;
    double *newarr;
    int i,j;

    newarr = (double*) malloc(nx*ny*nz*sizeof(double));
    nextarr = (double**) malloc(nx*ny*sizeof(double*));
    ret = (double***) malloc(nx*sizeof(double**));
    assert(newarr != NULL && nextarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
        nextarr[ny*i + j] = &newarr[i*(ny*nz)+j*nz];
      }
    }
    for(i = 0; i < nx; i++){
    ret[i] = &nextarr[i*ny];
    }
    return ret;
}

void free_contiguous_3dArray(double ***array){
    free(array[0][0]);
    free(array[0]);
    free(array);
}
*/


using namespace std; 

template <typename T> 
T*** new_contiguous_3dArray(int nx, int ny, int nz){
    T ***ret;
    T **nextarr;
    T *newarr;
    int i,j;

    newarr = (T*) malloc(nx*ny*nz*sizeof(T));
    nextarr = (T**) malloc(nx*ny*sizeof(T*));
    ret = (T***) malloc(nx*sizeof(T**));
    assert(newarr != NULL && nextarr != NULL && ret != NULL);

    for(i = 0; i < nx; i++){
      for(j = 0; j < ny; j++){
        nextarr[ny*i + j] = &newarr[i*(ny*nz)+j*nz];
      }
    }
    for(i = 0; i < nx; i++){
    ret[i] = &nextarr[i*ny];
    }
    return ret;
}

template <typename T> 
void free_contiguous_3dArray(T ***array){
    free(array[0][0]);
    free(array[0]);
    free(array);
}
#endif
