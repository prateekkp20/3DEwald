#include "libinclude.h"
#include "fundec.h"
#include "const.h"
#include "omp.h"
#define NUM_THREADS 10

//*Original loop, no parallelization
// double reci_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
//     double reci_energy=0;
//     for (int kx = -K; kx < K+1; kx++){
//         for (int ky = -K; ky < K+1; ky++){
//             for (int kz = -K; kz < K+1; kz++){
//                 if((kx==0) && (ky==0) && (kz==0))continue;
//                 complex<double> sg=0;
//                 complex<double> t(0,1);
//                 double G[3]={2*M_PI*kx/box[0][0], 2*M_PI*ky/box[1][1], 2*M_PI*kz/box[2][2]};
//                 for (int  i = 0; i < natoms; i++){
//                     double G_dot_r=G[0]*PosIons[i][0]+G[1]*PosIons[i][1]+G[2]*PosIons[i][2];
//                     complex<double> charge(ion_charges[i],0.0);
//                     sg+=charge*(cos(G_dot_r)+t*sin(G_dot_r));
//                 }
//                 double norm_sg = norm(sg);
//                 double mod_g = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
//                 reci_energy+=(1/pow(mod_g,2))*exp(-pow((mod_g)/(betaa*2),2))*norm_sg;
//             }
//         }
//     }
//     reci_energy*=(2*M_PI)/(box[0][0]*box[1][1]*box[2][2]);
//     return reci_energy;
// }

//* For reduction construct */
/*There is no workaround for the complex reduction, we have to do the ugly work*/
// double reci_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
//     double reci_energy=0;
//     for (int kx = -K; kx < K+1; kx++){
//         for (int ky = -K; ky < K+1; ky++){
//             for (int kz = -K; kz < K+1; kz++){
//                 if((kx==0) && (ky==0) && (kz==0))continue;
//                 complex<double> sg=0;
//                 complex<double> t(0,1);
//                 double G[3]={2*M_PI*kx/box[0][0], 2*M_PI*ky/box[1][1], 2*M_PI*kz/box[2][2]};
//                 #pragma omp parallel for reduction(+: sg)
//                     for (int  i = 0; i < natoms; i++){
//                         double G_dot_r=G[0]*PosIons[i][0]+G[1]*PosIons[i][1]+G[2]*PosIons[i][2];
//                         complex<double> charge(ion_charges[i],0.0);
//                         sg+=charge*(cos(G_dot_r)+t*sin(G_dot_r));
//                     }
//                 double norm_sg = norm(sg);
//                 double mod_g = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
//                 reci_energy+=(1/pow(mod_g,2))*exp(-pow((mod_g)/(betaa*2),2))*norm_sg;
//             }
//         }
//     }
//     reci_energy*=(2*M_PI)/(box[0][0]*box[1][1]*box[2][2]);
//     return reci_energy;
// }

//* synchronization construct critical
double reci_energy(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K){
    int nthreads;
    double reci_energy=0;
    omp_set_num_threads(NUM_THREADS);

    for (int kx = -K; kx < K+1; kx++){
        for (int ky = -K; ky < K+1; ky++){
            for (int kz = -K; kz < K+1; kz++){
                if((kx==0) && (ky==0) && (kz==0))continue;
                complex<double> sg=0;
                complex<double> t(0,1);
                double G[3]={2*M_PI*kx/box[0][0], 2*M_PI*ky/box[1][1], 2*M_PI*kz/box[2][2]};
                
                // Parallel Threading begins for the structure factor loop with goes over 0 to "natoms"
                #pragma omp parallel
                {
                    int i, id, nthrds;
                    complex<double> sum;
                    id = omp_get_thread_num();
                    nthrds = omp_get_num_threads();
                    if(id==0) nthreads=nthrds;

                    for (i = id,sum=0; i < natoms; i+=nthrds){
                        double G_dot_r=G[0]*PosIons[i][0]+G[1]*PosIons[i][1]+G[2]*PosIons[i][2];
                        complex<double> charge(ion_charges[i],0.0);
                        sum+=charge*(cos(G_dot_r)+t*sin(G_dot_r));
                    }
                    #pragma omp critical // this section stops the all other threads and first performs this action
                        sg+=sum;
                }
                double norm_sg = norm(sg);
                double mod_g = sqrt(G[0]*G[0]+G[1]*G[1]+G[2]*G[2]);
                reci_energy+=(1/pow(mod_g,2))*exp(-pow((mod_g)/(betaa*2),2))*norm_sg;
            }
        }
    }
    reci_energy*=(2*M_PI)/(box[0][0]*box[1][1]*box[2][2]);
    return reci_energy;
}