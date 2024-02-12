//***********************************************
// function to print coordinates in CONTCAR format!!
//***********************************************

#include "libinclude.h"
#include "const.h"
#include "mathdec.h"
#include "fundec.h"

// print coordinates in readable format

void print_coor(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep, char printmode, string filename){

	FILE *pFile;

	pFile = fopen(filename.c_str(), &printmode);

	if (printtrj == 0){

		//	cout<<"****Printing coordinates file****"<<endl;   //change by krishna
		fprintf(pFile, "System \n          1\n");

		for (int i = 0; i < 3; i++){
			fprintf(pFile, "%11.6f %11.6f %11.6f \n", boxcell[i][0], boxcell[i][1], boxcell[i][2]);
		}

		for (int i = 0; i < n_atomtype; i++){
			fprintf(pFile, "  %s", atomtype[i].c_str());
		}
		fprintf(pFile, "\n");
		for (int i = 0; i < n_atomtype; i++){
			fprintf(pFile, "  %d", natoms_type[i]);
		}
		fprintf(pFile, "\n");
	}

	fprintf(pFile, "Direct  configuration=        %d\n", MDstep);

	//	fprintf(pFile, "Cartesian  configuration=        %d\n", MDstep);
	// convert cartesian to direct coordinates

	double dir_coor[natoms][3];

	for (int i = 0; i < natoms; i++){
		dir_coor[i][0] = PosIons[i][0] / boxcell[0][0] - PosIons[i][1] * boxcell[0][1] / boxcell[0][0] / boxcell[1][1] + PosIons[i][2] * (boxcell[0][1] * boxcell[1][2] - boxcell[0][2] * boxcell[2][2]) / boxcell[0][0] / boxcell[1][1] / boxcell[2][2];
		dir_coor[i][1] = PosIons[i][1] / boxcell[1][1] - PosIons[i][2] * boxcell[1][2] / boxcell[1][1] / boxcell[2][2];
		dir_coor[i][2] = PosIons[i][2] / boxcell[2][2];

		fprintf(pFile, "%11.6f %11.6f %11.6f \n", dir_coor[i][0], dir_coor[i][1], dir_coor[i][2]);
	}

	/*

		for(int i=0;i<natoms;i++)
		{
			fprintf(pFile, "%11.6f %11.6f %11.6f \n",PosIons[i][0],PosIons[i][1],PosIons[i][2]);
		}

	  */

	fclose(pFile);
}

void print_carcoor(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, string *atomtype, int printtrj, int MDstep, char printmode, string filename){

	cout << "****Printing CONTCAR****" << endl;

	FILE *pFile;

	pFile = fopen(filename.c_str(), &printmode);

	if (printtrj == 0){
		fprintf(pFile, "System \n          1\n");

		for (int i = 0; i < 3; i++){
			fprintf(pFile, "%11.6f %11.6f %11.6f \n", boxcell[i][0], boxcell[i][1], boxcell[i][2]);
		}

		for (int i = 0; i < n_atomtype; i++){
			fprintf(pFile, "  %s", atomtype[i].c_str());
		}
		fprintf(pFile, "\n");
		for (int i = 0; i < n_atomtype; i++){
			fprintf(pFile, "  %d", natoms_type[i]);
		}
		fprintf(pFile, "\n");
	}

	//	fprintf(pFile, "Direct  configuration=        %d\n", MDstep);

	fprintf(pFile, "Cartesian  configuration=        %d\n", MDstep);

	for (int i = 0; i < natoms; i++){
		fprintf(pFile, "%11.6f %11.6f %11.6f \n", PosIons[i][0], PosIons[i][1], PosIons[i][2]);
	}

	fclose(pFile);
}

// Print coordinates in outputfile

void printCoor(double **PosIons, int natoms, string *type){
	// cout<<"*****Coordinates of atoms in Angstrom****"<<endl<<endl;  //changes by krishna

	for (int i = 0; i < natoms; i++){
		printf("%d  %s  %11.6f %11.6f %11.6f \n", i + 1, type[i].c_str(), PosIons[i][0], PosIons[i][1], PosIons[i][2]);
	}

	cout << endl;
}

void printVel(double **Vel, int natoms, string *type){
	// cout<<"*****Velocities of atoms in Angs/fs******"<<endl<<endl;  //changes by krishna

	for (int i = 0; i < natoms; i++){
		printf("%d  %s  %11.6f %11.6f %11.6f \n", i + 1, type[i].c_str(), Vel[i][0], Vel[i][1], Vel[i][2]);
	}

	cout << endl;
}

void printFor(double **ForceIons, int natoms, string *type){
	// cout<<"*****Forces on atoms in eV/Angs****"<<endl<<endl;    //changes by krishna

	for (int i = 0; i < natoms; i++){
		printf("%d  %s  %11.6f %11.6f %11.6f \n", i + 1, type[i].c_str(), ForceIons[i][0] / KJEV, ForceIons[i][1] / KJEV, ForceIons[i][2] / KJEV);
	}

	cout << endl;
}

void printVel(double **vel, int natoms){
	ofstream VelOut;
	// cout<<"****Printing velocities*****"<<endl;			//changes by krishna

	VelOut.open("velocities.out");

	for (int i = 0; i < natoms; i++){
		for (int j = 0; j < 3; j++){
			VelOut << vel[i][j] << "  \t";
		}
		VelOut << endl;
	}

	VelOut.close();
}

void printprobVel(double **vel, int natoms, float *mass, float Temp){
	ofstream VelOut;

	VelOut.open("veldist.out");

	double ldist = 0.0001;
	double RMax = 0.05;
	double RMin = -0.05;
	int tmp;

	int ndist = (RMax - RMin) / ldist + 2;

	int Dist[ndist];

	for (int i = 0; i < ndist; i++){
		Dist[i] = 0;
	}

	for (int i = 0; i < natoms; i++){
		for (int j = 0; j < 3; j++){

			tmp = ceil((vel[i][j] - RMin) / ldist);

			if (tmp < 0){
				Dist[0] = Dist[0] + 1;
			}
			else if (tmp >= (ndist - 1)){
				Dist[ndist - 1] = Dist[ndist - 1] + 1;
			}
			else
				Dist[tmp] = Dist[tmp] + 1;
		}
	}

	for (int i = 0; i < ndist; i++){
		VelOut << RMin - ldist / 2 + i * ldist << "  \t" << Dist[i] << endl;
	}

	VelOut.close();
}