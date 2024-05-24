#include "libinclude.h"
#include "const.h"

#define REAL 0
#define IMAG 1
// Disable this declaration if openmp parallelization is not required, would not be helpful for smaller systems
#define ENABLE_OMP 11
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
    long double G[3][3], m[3];
    long double **u,**x_direc, **y_direc, **z_direc;
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
    
    // Volume Calculations
    long double A[3];
    long double C[3]={box[2][0],box[2][1],box[2][2]};
    crossProduct(box[0],box[1],A);
    long double volume = dotProduct(A,C);
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    fftw_plan p;
    p = fftw_plan_dft_3d(K,K,K, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Calculating the reciprocal vectors
    crossProduct(box[1],box[2],G[0]);
    crossProduct(box[2],box[0],G[1]);
    crossProduct(box[0],box[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;
    #if defined ENABLE_OMP
        omp_set_num_threads(thread::hardware_concurrency());
        #pragma omp parallel for simd
    #endif
    // Calculating the fractional coordinates
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=K*dotProductu(PosIons[i],G[j]);
        }
    }
    #if defined ENABLE_OMP
        #pragma omp parallel for simd
    #endif
    // Calculating the cofficients in the x,y and z directions for the Q Matrix
    for (int i = 0; i < natoms; i++){
        // for X direction
        for (int  k1 = 0; k1 < K; k1++){
            x_direc[i][k1]=0;
            for (int  n1 = -n_max; n1 < n_max+1; n1++){
                x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K,n);
            }
        }
        // for Y direction
        for (int  k2 = 0; k2 < K; k2++){
            y_direc[i][k2]=0;
            for (int  n2 = -n_max; n2 < n_max+1; n2++){
                y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K,n);
            }
        }
        // for Z direction
        for (int  k3 = 0; k3 < K; k3++){
            z_direc[i][k3]=0;
            for (int  n3 = -n_max; n3 < n_max+1; n3++){
                z_direc[i][k3]+=M_n(u[i][2]-k3-n3*K,n);
            }
        }
    }
    #if defined ENABLE_OMP
        #pragma omp for collapse(3) 
    #endif
    // initializing the "in" vector with zero values
    for (int tx = 0; tx < K; tx++){
        for (int ty = 0; ty < K; ty++){
            for (int tz = 0; tz < K; tz++){
                in[tx * (K * K) + ty * K + tz][0] = 0.0;
            }
        }
    }
    #if defined ENABLE_OMP
        #pragma omp for simd
    #endif
    // Final Q Matrix
    for (int j = 0; j < natoms; j++){
    if (natoms == 0)
        continue;
        for (int tx = 0; tx < K; tx++){
            if (x_direc[j][tx] == 0)continue;

            for (int ty = 0; ty < K; ty++){
                if (y_direc[j][ty] == 0)continue;

                for (int tz = 0; tz < K; tz++){
                    if (z_direc[j][tz] == 0)continue;

                    in[tx * (K * K) + ty * K + tz][0] += ion_charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
                }
            }
        }
    }

    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup();
    long double energy=0;
    double constant=(M_PI*M_PI)/(betaa*betaa);
    // collapse doesn't makes a difference here much; dynamic and runtime give the same time 
    int ii,jj,kk;
    #if defined ENABLE_OMP
        #pragma omp parallel for schedule(runtime) reduction(+: energy) collapse(3)
    #endif
    for (int i = -M; i < M+1; i++){
        for (int j = -M; j< M+1; j++){
            for (int k = -M; k < M+1; k++){
                if(i==0&&j==0&&k==0)continue;
                for (int t = 0; t < 3; t++)                {
                    m[t]=i*G[0][t]+j*G[1][t]+k*G[2][t];    
                }                
                long double m2=dotProduct(m,m);
                if(i<0) ii=K+i;
                else ii=i;
                if(j<0) jj=K+j;
                else  jj=j;
                if(k<0) kk=K+k;
                else  kk=k;
                int temp=ii * (K * K) + jj * K + kk;
                long double norm_FQ=out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];
                energy += norm_FQ*exp(-m2*constant)*norm(B(i,n,K)*B(j,n,K)*B(k,n,K))/m2;
            }
        }
    }
    energy/=(2*M_PI*volume);
    return energy;
}