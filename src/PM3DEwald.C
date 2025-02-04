#include "libinclude.h"
#include "const.h"
#include "fundec.h"
#include "header.h"

#define REAL 0
#define IMAG 1
// Disable this declaration if openmp parallelization is not required, would not be helpful for smaller systems
#define ENABLE_OMP 11

double PM3DEwald(double **PosIons, double *ion_charges, int natoms, double betaa, double **box, int Grid, int M, int n){
    // n: order of b-spline interpolation
   // initializing the new variables
    // double G[3][3];
    double **u,**x_direc, **y_direc, **z_direc;
    u= new double * [natoms];
    x_direc= new double * [natoms];
    y_direc= new double * [natoms];
    z_direc= new double * [natoms];

    for (int  i = 0; i < natoms; i++){
        u[i] = new double  [3];
        x_direc[i] = new double  [Grid];
        y_direc[i] = new double  [Grid];
        z_direc[i] = new double  [Grid];
    }

    double L1 = box[0][0];
    double L2 = box[1][1];
    double L3 = box[2][2];
    int n_max=1;

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid*Grid);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *Grid*Grid*Grid);
    fftw_plan p;
    p = fftw_plan_dft_3d(Grid,Grid,Grid, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for
    #endif
    // Calculating the fractional coordinates
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=Grid*dotProduct(PosIons[i],G[j]);
        }
    }

    #if defined ENABLE_OMP
        #pragma omp parallel for
    #endif
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < Grid; k1++){
            x_direc[i][k1]=0;
            for (int  n1 = -n_max; n1 < n_max+1; n1++){
                x_direc[i][k1]+=M_n(u[i][0]-k1-n1*Grid,n);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < Grid; k2++){
            y_direc[i][k2]=0;
            for (int  n2 = -n_max; n2 < n_max+1; n2++){
                y_direc[i][k2]+=M_n(u[i][1]-k2-n2*Grid,n);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < Grid; k3++){
            z_direc[i][k3]=0;
            for (int  n3 = -n_max; n3 < n_max+1; n3++){
                z_direc[i][k3]+=M_n(u[i][2]-k3-n3*Grid,n);
            }
        }
    }

    // initializing the "in" vector with zero values
    for (int tx = 0; tx < Grid; tx++){
        for (int ty = 0; ty < Grid; ty++){
            for (int tz = 0; tz < Grid; tz++){
                in[tx * (Grid * Grid) + ty * Grid + tz][0] = 0.0;
            }
        }
    }

    #if defined ENABLE_OMP
        #pragma omp parallel for 
    #endif
    // Final Q Matrix
    for (int j = 0; j < natoms; j++){
        if (ion_charges[j] == 0)continue;
        for (int tx = 0; tx < Grid; tx++){
            if (x_direc[j][tx] == 0)continue;

            for (int ty = 0; ty < Grid; ty++){
                if (y_direc[j][ty] == 0)continue;

                for (int tz = 0; tz < Grid; tz++){
                    if (z_direc[j][tz] == 0)continue;
                    #if defined ENABLE_OMP
                        #pragma omp atomic update
                    #endif
                    in[tx * (Grid * Grid) + ty * Grid + tz][0] += ion_charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
                }
            }
        }
    }

    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup();

    double energy=0;
    // double constant=(M_PI*M_PI)/(betaa*betaa);
    // collapse doesn't makes a difference here much; dynamic and runtime give the same time 
    int ii,jj,kk;
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: energy) collapse(3)
    #endif
    for (int i = -M; i < M+1; i++){
        for (int j = -M; j< M+1; j++){
            for (int k = -M; k < M+1; k++){
                if(i==0&&j==0&&k==0)continue;
                double m[3];
                for (int t = 0; t < 3; t++){
                    m[t]=i*G[0][t]+j*G[1][t]+k*G[2][t];    
                }
                double m2=dotProduct(m,m);
                int ic,jc,kc;
                if(i<0) {ii=Grid+i;ic=(2*M+1)+i;}
                else {ii=i;ic=i;}
                if(j<0) {jj=Grid+j;jc=(2*M+1)+j;}
                else  {jj=j;jc=j;}
                if(k<0) {kk=Grid+k;kc=(2*M+1)+k;}
                else  {kk=k;kc=k;}

                int temp=ii * (Grid * Grid) + jj * Grid + kk;
                int tempexpfactor = ic * ((2*M+1) * (2*M+1)) + jc * (2*M+1) + kc;

                double norm_FQ=out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];

                // energy += norm_FQ*exp(-m2*constant)*norm(B(i,n,Grid)*B(j,n,Grid)*B(k,n,Grid))/m2;
                // energy += norm_FQ*ExpFactor[tempexpfactor]*norm(B(i,n,Grid)*B(j,n,Grid)*B(k,n,Grid));
                // energy += norm_FQ*ExpFactor[tempexpfactor]*norm(CoeffX[ic]*CoeffY[jc]*CoeffZ[kc]);
                energy += norm_FQ*ExpFactor[tempexpfactor];
            }
        }
    }
    energy/=(2*M_PI*volume);
    return energy;
}