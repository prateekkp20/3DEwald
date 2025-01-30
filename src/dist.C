//**********************************************************************
// function to calculate distances according to minimum image convention 
//**********************************************************************

#include "libinclude.h" 

double dist(double **PosIons, int atom1, int atom2, double **box){

	double Dx, Dy, Dz;
	double Dx1, Dy1, Dz1;
    double distatoms;
	
	Dx = PosIons[atom1][0] - PosIons[atom2][0];
	Dy = PosIons[atom1][1] - PosIons[atom2][1];
	Dz = PosIons[atom1][2] - PosIons[atom2][2];

	double a  = ceil(Dx/box[0][0]-0.5), b = ceil(Dy/box[1][1]-0.5), c=ceil(Dz/box[2][2]-0.5);
	
    Dx1 = Dx - box[0][0]*a - box[1][0]*b - box[2][0]*c;
    Dy1 = Dy - box[0][1]*a - box[1][1]*b - box[2][1]*c;
    Dz1 = Dz - box[0][2]*a - box[1][2]*b - box[2][2]*c;
	
	distatoms = sqrt(Dx1*Dx1 + Dy1*Dy1 + Dz1*Dz1);
	return distatoms;
	
}