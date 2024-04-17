#include "libinclude.h"
#include "const.h"

double M_n(double u, int n){
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

complex<double>B(int m, int n, int K){
    // const double M_PI = acos(-1.0);
    const complex<double> t(0.0, 1.0);
    complex<double> bi_mi=exp((2*M_PI*m*(n-1))/K*t);
    complex<double> denox;
    for (int f = 0; f < n-1; f++){
    denox+=M_n(f+1,n)*exp((2*M_PI*m*f)/K*t);
    }
    bi_mi/=denox;
    return bi_mi;
}

double dotProduct(double *v1, double *v2) {
    double result = 0.0;
    for (int i = 0; i < 3; ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

void crossProduct(double *v_A, double *v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

double bspline(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int M, int order){

    double **G,**u,**x_direc, **y_direc, **z_direc;
    G= new double * [3];
    u= new double * [natoms];
    x_direc= new double * [natoms];
    y_direc= new double * [natoms];
    z_direc= new double * [natoms];

    for (int  i = 0; i < natoms; i++){
        G[i] = new double  [3];
        u[i] = new double  [3];
        x_direc[i] = new double  [3];
        y_direc[i] = new double  [3];
        z_direc[i] = new double  [3];
    }

    float L1 = box[0][0];
    float L2 = box[1][1];
    float L3 = box[2][2];
    int n_max=1;
    int n=5; //order of b-spline interpolation

    fftw_complex *in;   //input variable using standard fftw syntax
    fftw_complex *out;	// output variable

    float L[3]={L1,L2,L3};
    double G[3]={2*M_PI/box[0][0], 2*M_PI/box[1][1], 2*M_PI/box[2][2]};

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K*K*K);
    fftw_plan p;
    p = fftw_plan_dft_3d(K,K,K, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Calculating the fractional coordinates
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=K*dotProduct(PosIons[i],G[j]);
        }
    }

    // Calculating the Q Matrix
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
                y_direc[i][k2]+=M_n(u[i][0]-k2-n2*K,n);
            }
        }

        // for Z direction
        for (int  k3 = 0; k3 < K; k3++){
            x_direc[i][k3]=0;
            for (int  n3 = -n_max; n3 < n_max+1; n3++){
                z_direc[i][k3]+=M_n(u[i][0]-k3-n3*K,n);
            }
        }
    }
    
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


}