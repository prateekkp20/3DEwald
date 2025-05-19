#include "libinclude.h"
#include "const.h"
#include "fundec.h"

const complex<double> t(0.0, 1.0);

template<typename T>
double dotProduct(T v1,T v2) {
    double result = 0.0;
    for (int i = 0; i < 3; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}
template double dotProduct<float*>(float* v1, float* v2);
template double dotProduct<double*>(double* v1, double* v2);

void crossProduct(double *v_A, double *v_B, double *out){
   out[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   out[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   out[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
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
    complex<double> bi_mi=exp((2*M_PI*M_PI*(n-1)/K)*t);
    complex<double> denox;
    for (int f = 0; f < n-1; f++){
    denox+=M_n(f+1,n)*exp((2*M_PI*m*f/K)*t);
    }
    bi_mi/=denox;
    return bi_mi;
}

double error(double a, double b){
    return abs(a-b)/a;
}