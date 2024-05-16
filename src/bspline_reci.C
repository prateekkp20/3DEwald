#include "libinclude.h"
#include "const.h"

#define REAL 0
#define IMAG 1
const long double pi = M_PI;

long double M_n(long double u, int n){
    if(n<2)return 0;
    else if(n==2){
        if(u<0 || u>n) return 0;
        else{
            return 1-abs(u-1);
        }
    }
    else{
        if(u<0 || u>n) return 0;
        else{
            return (u*M_n(u,n-1)/(n-1))+((n-u)*M_n(u-1,n-1)/(n-1));
        }
    }
}

complex<long double>B(int m, int n, int K){
    const complex<long double> t(0.0, 1.0);
    complex<long double> bi_mi=exp((2*pi*pi*(n-1))/K*t);
    complex<long double> denox;
    for (int f = 0; f < n-1; f++){
    denox+=M_n(f+1,n)*exp((2*pi*m*f)/K*t);
    }
    bi_mi/=denox;
    return bi_mi;
}

long double dir(long double u, int t, int K, int n, int n_max){
    long double direc=0;
    for (int  n1 = -n_max; n1 < n_max+1; n1++){
        direc+=M_n(u-t-n1*K,n);
    }
    return direc;
}

long double dotProductu(double *v1,long double *v2) {
    double result = 0.0;
    for (int i = 0; i < 3; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

long double dotProduct(long double *v1,long double *v2) {
    double result = 0.0;
    for (int i = 0; i < 3; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

void crossProduct(float *v_A, float *v_B,long double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

double bspline(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int M, int n){
    // n: order of b-spline interpolation
    // initializing the new variables
    // fftw_init_threads();
    cout<<fixed<<setprecision(10);
    long double G[3][3], m[3];
    long double **u,**x_direc, **y_direc, **z_direc;
    // long double G[3][3], u[natoms][3], x_direc[natoms][K], y_direc[natoms][K], z_direc[natoms][K], m[3];
    u= new long double * [natoms];
    x_direc= new long double * [natoms];
    y_direc= new long double * [natoms];
    z_direc= new long double * [natoms];

    for (int  i = 0; i < natoms; i++){
        u[i] = new long double  [3];
        x_direc[i] = new long double  [K];
        y_direc[i] = new long double  [K];
        z_direc[i] = new long double  [K];
    }

    float L1 = box[0][0];
    float L2 = box[1][1];
    float L3 = box[2][2];
    int n_max=2;

    fftw_complex *in;   // input variable using standard fftw syntax
    fftw_complex *out;	// output variable

    float L[3]={L1,L2,L3};
    long double volume = L1*L2*L3;
    // omp_set_num_threads(thread::hardware_concurrency());
    // omp_set_num_threads(8);
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    fftw_plan p;
    // #pragma omp critical (make_plan)
    p = fftw_plan_dft_3d(K,K,K, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_plan_with_nthreads(thread::hardware_concurrency());

    // Calculating the reciprocal vectors
    crossProduct(box[1],box[2],G[0]);
    crossProduct(box[2],box[0],G[1]);
    crossProduct(box[0],box[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;

    // Calculating the fractional coordinates
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=K*dotProductu(PosIons[i],G[j]);
        }
    }

    // Calculating the Q Matrix
    // for (int i = 0; i < natoms; i++){
    //     // for X direction
    //     for (int  k1 = 0; k1 < K; k1++){
    //         x_direc[i][k1]=0;
    //         for (int  n1 = -n_max; n1 < n_max+1; n1++){
    //             x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K,n);
    //         }
    //     }

    //     // for Y direction
    //     for (int  k2 = 0; k2 < K; k2++){
    //         y_direc[i][k2]=0;
    //         for (int  n2 = -n_max; n2 < n_max+1; n2++){
    //             y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K,n);
    //         }
    //     }

    //     // for Z direction
    //     for (int  k3 = 0; k3 < K; k3++){
    //         z_direc[i][k3]=0;
    //         for (int  n3 = -n_max; n3 < n_max+1; n3++){
    //             z_direc[i][k3]+=M_n(u[i][2]-k3-n3*K,n);
    //         }
    //     }
    // }

    // initializing the "in" vector with zero values
    for (int tx = 0; tx < K; tx++){
        for (int ty = 0; ty < K; ty++){
            for (int tz = 0; tz < K; tz++){
                in[tx * (K * K) + ty * K + tz][0] = 0.0;
            }
        }
    }
    for (int j = 0; j < natoms; j++){
    if (natoms == 0)
        continue;
        for (int tx = 0; tx < K; tx++){
            long double x_dir = dir(u[j][0],tx, K, n, n_max);
            if (x_dir == 0)continue;
            // if (x_direc[j][tx] == 0)continue;

            for (int ty = 0; ty < K; ty++){
                long double y_dir = dir(u[j][1],ty, K, n, n_max);
                if (y_dir == 0)continue;

                for (int tz = 0; tz < K; tz++){
                    long double z_dir = dir(u[j][2],tz, K, n, n_max);
                    if (z_dir == 0)continue;

                    // in[tx * (K * K) + ty * K + tz][0] += ion_charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
                    in[tx * (K * K) + ty * K + tz][0] += ion_charges[j] * x_dir * y_dir * z_dir;
                }
            }
        }
    }

    fftw_execute(p);
    #pragma omp critical (FFTW)
    fftw_destroy_plan(p);
    fftw_cleanup();
    // fftw_cleanup_threads();
    long double energy=0;
    double constant=(M_PI*M_PI)/(betaa*betaa);
    int i,j,k,ii,jj,kk;
    // collapse doesn't makes a difference here much; dynamic and runtime give the same time
    // #pragma omp parallel for reduction(+: energy)
    // #pragma omp parallel for schedule(runtime) reduction(+: energy)
    // #pragma omp parallel for schedule(runtime) reduction(+: energy) collapse(3)
    for (i = -M; i < M+1; i++){
        for (j = -M; j< M+1; j++){
            for (k = -M; k < M+1; k++){
                if(i<0) ii=K+i;
                else ii=i;
                if(j<0) jj=K+j;
                else  jj=j;
                if(k<0) kk=K+k;
                else  kk=k;
                if(i==0&&j==0&&k==0)continue;
                m[0]=i*G[0][0]+j*G[1][0]+k*G[2][0];
                m[1]=i*G[0][1]+j*G[1][1]+k*G[2][1];
                m[2]=i*G[0][2]+j*G[1][2]+k*G[2][2];
                long double m2=dotProduct(m,m);
                int temp=ii * (K * K) + jj * K + kk;
                long double norm_FQ=out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];
                energy += norm_FQ*exp(-m2*constant)*norm(B(i,n,K)*B(j,n,K)*B(k,n,K))/m2;
            }
        }
    }
    energy/=(2*M_PI*volume);
    return energy;
}