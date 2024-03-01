#include "libinclude.h"
#include "fundec.h"
#include "omp.h"
#define NUM_THREADS 10
#define PAD 8

/* false sharing */
// double real_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box){
//     int nthreads;
//     double sum[NUM_THREADS][PAD],real_energy=0;
//     omp_set_num_threads(NUM_THREADS);

//     #pragma omp parallel
//     {
//         int id, nthrds;
//         id = omp_get_thread_num();
//         nthrds = omp_get_num_threads();
//         if(id==0) nthreads=nthrds;
//         int i;
//         for (i = id,sum[id][0]=0; i < natoms; i+=nthrds){
//             for (int j = 0; j < i; j++){
//                 if(i!=j){
//                     double modR=dist(PosIons,i,j,box);
//                     sum[id][0]+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
//                 }
//             }
//         }
//     }
//     for(int j=0;j<nthreads;j++){
//         real_energy+=sum[j][0];
//     }
//     return real_energy;
// }

/* synchromization construct critical */
double real_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box){
    int nthreads;
    double real_energy=0;
    omp_set_num_threads(NUM_THREADS);

    #pragma omp parallel
    {
        int i, id, nthrds;
        double sum;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        if(id==0) nthreads=nthrds;
        for (i = id,sum=0; i < natoms; i+=nthrds){
            for (int j = 0; j < i; j++){
                if(i!=j){
                    double modR=dist(PosIons,i,j,box);
                    sum+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                }
            }
        }
        #pragma omp critical
            real_energy+=sum;
    }
    return real_energy;
}

/*Original loop, no parallelization*/
// double real_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box){
//     double real_energy=0;
//     for (int i = 0; i < natoms; i++){
//         for (int j = 0; j < i; j++){
//             if(i!=j){
//                 double modR=dist(PosIons,i,j,box);
//                 real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
//             }
//         }
//     }
    
//     return real_energy;
// }