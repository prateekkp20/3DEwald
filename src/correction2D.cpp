#include "libinclude.h"
#include "fundec.h"
#include "const.h"

double dipoleCorrection(double **PosIons, float *ion_charges, int natoms, float **box){
    // Volume Calculations
    long double A[3];
    long double C[3]={box[2][0],box[2][1],box[2][2]};
    crossProduct(box[0],box[1],A);
    long double volume = dotProduct(A,C);

    long double M[3]={0,0,0};
    for (int i = 0; i < natoms; i++){
        for(int j = 0; j < 3; j++){
            M[j]+=ion_charges[i]*PosIons[i][j];
        }
    }
    double J = M[2]*M[2]*2*M_PI/volume;
    // double J = dotProduct(M,M)*2*M_PI/(volume*3);
    return J;
}
