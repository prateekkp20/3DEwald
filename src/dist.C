//***********************************************
// function to calculate distances according to minimum image convention 
//***********************************************

#include "libinclude.h" 

double dist(double **PosIons, int atom1, int atom2, float **box){

	double Dx, Dy, Dz;
	double Dx1, Dy1, Dz1;
    double distatoms;
	
	Dx = PosIons[atom1][0] - PosIons[atom2][0];
	Dy = PosIons[atom1][1] - PosIons[atom2][1];
	Dz = PosIons[atom1][2] - PosIons[atom2][2];
	
    Dx1 = Dx - box[0][0]*ceil(Dx/box[0][0]-0.5) - box[1][0]*ceil(Dy/box[1][1]-0.5) - box[2][0]*ceil(Dz/box[2][2]-0.5);
	Dy1 = Dy - box[0][1]*ceil(Dx/box[0][0]-0.5) - box[1][1]*ceil(Dy/box[1][1]-0.5) - box[2][1]*ceil(Dz/box[2][2]-0.5);
	Dz1 = Dz - box[0][2]*ceil(Dx/box[0][0]-0.5) - box[1][2]*ceil(Dy/box[1][1]-0.5) - box[2][2]*ceil(Dz/box[2][2]-0.5);
	
	distatoms = sqrt(pow(Dx1,2) + pow(Dy1,2) + pow(Dz1,2));
	return distatoms;
	
}
