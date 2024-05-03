#include "libinclude.h"
#include "const.h"

#define REAL 0
#define IMAG 1

void scalarProductMat(double **mat, double k){
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            mat[i][j] *= k;    
}

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

void crossProduct(float *v_A, float *v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

double bspline(double **PosIons, float *ion_charges, int natoms, double betaa, float **box, int K, int M, int n){
    // initializing the new variables
    double **G,**u,**x_direc, **y_direc, **z_direc, *m;
    m= new double [3];
    G= new double * [3];
    u= new double * [natoms];
    x_direc= new double * [natoms];
    y_direc= new double * [natoms];
    z_direc= new double * [natoms];

    for (int  i = 0; i < 3; i++){
        G[i] = new double  [3];
    }
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
    int n_max=2;
    // int n=order; //order of b-spline interpolation

    fftw_complex *in;   //input variable using standard fftw syntax
    fftw_complex *out;	// output variable

    float L[3]={L1,L2,L3};
    double volume = L1*L2*L3;

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
            z_direc[i][k3]=0;
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

    crossProduct(box[1],box[2],G[0]);
    crossProduct(box[2],box[0],G[1]);
    crossProduct(box[0],box[1],G[2]);
    scalarProductMat(G,1/volume);

    double energy=0;
    double constant=(M_PI*M_PI)/(betaa*betaa);
    int ii,jj,kk;
    for (int i = -M; i < M+1; i++){
        if(i<0) ii=K+i;
        else ii=i;
        for (int j = -M; j< M+1; j++){
            if(j<0) jj=K+j;
            else  jj=j;
            for (int k = -M; k < M+1; k++){
                if(k<0) kk=K+k;
                else  kk=k;
                if(i==0&&j==0&&k==0)continue;
                m[0]=i*G[0][0]+j*G[1][0]+k*G[2][0];
                m[1]=i*G[0][1]+j*G[1][1]+k*G[2][1];
                m[2]=i*G[0][2]+j*G[1][2]+k*G[2][2];
                double m2=dotProduct(m,m);
                int temp=ii * (K * K) + jj * K + kk;
                // cout<<temp<<"\n";
                double norm_FQ=out[temp][REAL]*out[temp][REAL]+out[temp][IMAG]*out[temp][IMAG];
                cout<<norm_FQ<<"\n";
                energy += norm_FQ*exp(-m2*constant)*norm(B(i,n,K)*B(j,n,K)*B(k,n,K))/m2;
            }
        }
    }
    energy/=(2*M_PI*volume);

    // for (int  i = 0; i < 3; i++){
    //     delete [] G[i];
    // }
    // for (int  i = 0; i < natoms; i++){
    //     delete [] u[i];
    //     delete [] x_direc[i];
    //     delete [] y_direc[i];
    //     delete [] z_direc[i];
    // }
    // delete [] G;
	// delete [] u;
	// delete [] x_direc;
	// delete [] y_direc;
	// delete [] z_direc;
    // cout<<energy;
    return energy;

}