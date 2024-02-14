#include "libinclude.h"
#include "fundec.h"

double real_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box){
    double real_energy=0;
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < i; j++){
            if(i!=j){
                double modR=dist(PosIons,i,j,box);
                // cout<<modR<<"\n";
                real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                // cout<<real_energy<<"\n";
            }
        }
    }
    
    return real_energy;
}