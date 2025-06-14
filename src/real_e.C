#include "libinclude.h"
#include "fundec.h"
#include "const.h"

#define PAD 8
//Uncomment only one of them
// #define NAIVE 1
// #define FALSE_SHARING_AND_PADDING 2
#define REDUCTION 3
// #define SYNCHRONIZATION_CONSTRUCT 4


#if defined FALSE_SHARING_AND_PADDING
//* false sharing and padding
    double real_energy(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
        int nthreads;
        double sum[NUM_THREADS][PAD],real_energy=0;
        #pragma omp parallel
        {
            int id, nthrds;
            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            if(id==0) nthreads=nthrds;
            int i;
            for (i = id,sum[id][0]=0; i < natoms; i+=nthrds){
                for (int j = 0; j < i; j++){
                        double modR=dist(PosIons,i,j,box);
                        if(modR>cutoff)continue;

                        /* Erfc Approximations*/
                        // double val = betaa*modR;
                        // double exp_x2 = exp(-val*val);
                        // double t, t1 =  t  = 1/(1+0.3275911*val);
                        // double erfcx = exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
                        // sum[id][0]+=(ion_charges[i]*ion_charges[j]*erfcx)/modR;

                        sum[id][0]+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                }
            }
        }
        for(int j=0;j<nthreads;j++){
            real_energy+=sum[j][0];
        }
        return real_energy;
    }

#elif defined SYNCHRONIZATION_CONSTRUCT
//* synchromization construct critical
    double real_energy(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
        int nthreads;
        double real_energy=0;
        #pragma omp parallel
        {
            int i, id, nthrds;
            double sum;

            id = omp_get_thread_num();
            nthrds = omp_get_num_threads();
            if(id==0) nthreads=nthrds;

            for (i = id,sum=0; i < natoms; i+=nthrds){
                for (int j = 0; j < i; j++){
                        double modR=dist(PosIons,i,j,box); // calculating the minimum distance between the given i and j atoms
                        if(modR>cutoff)continue;

                        /* Erfc Approximations*/
                        // double val = betaa*modR;
                        // double exp_x2 = exp(-val*val);
                        // double t, t1 =  t  = 1/(1+0.3275911*val);
                        // double erfcx = exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
                        // sum+=(ion_charges[i]*ion_charges[j]*erfcx)/modR;
                        
                        sum+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                }
            }
            #pragma omp critical // this section stops the all other threads and first performs this action
                real_energy+=sum;
        }
        return real_energy;
    }

#elif defined NAIVE
//*Original loop, no parallelization
    double real_energy(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
        double real_energy=0;
        for (int i = 0; i < natoms; i++){
            for (int j = 0; j < i; j++){
                    double modR=dist(PosIons,i,j,box);
                    if(modR>cutoff)continue;

                    /* Erfc Approximations*/
                    // double val = betaa*modR;
                    // double exp_x2 = exp(-val*val);
                    // double t, t1 =  t  = 1/(1+0.3275911*val);
                    // double erfcx = exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
                    // real_energy+=(ion_charges[i]*ion_charges[j]*erfcx)/modR;

                    real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
            }
        }
        
        return real_energy;
    }

#elif defined REDUCTION
//* For reduction construct
    double real_energy(double *PosIons, double *ion_charges, int natoms, double betaa, double **box, double cutoff){
        double real_energy=0;
        #pragma omp parallel for simd schedule(runtime) reduction(+: real_energy)
            for (int i = 0; i < natoms; i++){
                #pragma omp SIMD
                for (int j = 0; j < i; j++){
                        double modR=dist(PosIons,i,j,box);
                        if(modR>cutoff)continue;

                        // Erfc Approximations
                        double val = betaa*modR;
                        double exp_x2 = exp(-val*val);
                        double t, t1 =  t  = 1/(1+0.3275911*val);
                        double erfcx = exp_x2*(0.254829592*t - 0.284496736*(t*=t1) + 1.421413741*(t*=t1) - 1.453152027*(t*=t1) + 1.061405429*(t*=t1));
                        real_energy+=(ion_charges[i]*ion_charges[j]*erfcx)/modR;

                        // real_energy+=(ion_charges[i]*ion_charges[j]*erfc(betaa*modR))/modR;
                }
            }
        
        return real_energy;
    }

#else
    #error "Please define the method of parallelization"

#endif