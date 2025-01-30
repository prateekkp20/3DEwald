#include "libinclude.h"
#include "const.h"

double selfe(int n_atomtype, int *natoms_type, float *chargs, double betaa){
    double self_energy=0;
    for (int i = 0; i < n_atomtype; i++){
        self_energy+=natoms_type[i]*pow(chargs[i],2);
    }
    return self_energy*(-betaa)/sqrt(M_PI);
}