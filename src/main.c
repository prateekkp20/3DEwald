#include "libinclude.h" 
#include "potdec.h"
#include "mathdec.h"
#include "fundec.h"
#include "const.h"
#include "fftw3.h"
#include "complex"
#include "header.h"

double *ExpFactor, *PreExpFactor;
double G[3][3];
double volume;
complex<double> *CoeffX,*CoeffY,*CoeffZ;

int main(int argc, char **argv){

	string filename;
	string posfile;
	string garbage, garbage1;
	char garbage2[200];
	int tmp;
	int i, j, MC, MD;
	int  nMD, md_print;
	float Temp;
	int  MDtime, MDeq;
	float dt;
	int genpairlist;
	int randnum;
	int hyb, NMCMD, rdvel;
	float tmp2;
	double unitzer=(AVG*COL*CHARELEC*CHARELEC/rtoa)/1000/4.184;

	//get the name of the file for input variables, default is input.in 
	if(argc <2)
		filename = "input.in";
	else
		filename = argv[1];

	///////////////////////////////////////////////
	//////initialize variables////////////////////

	ifstream InputIn(filename.c_str(),ios::in);
	if(!InputIn){
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	getline(InputIn, garbage);

	cout<<"Project Name: "<<garbage<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>posfile;

	InputIn>>garbage>>garbage;
	InputIn>>genpairlist;

	InputIn>>garbage>>garbage;
	InputIn>>Temp;

	// cout<<"Temperature of the system:\t\t"<<Temp<<" K"<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>hyb;

	InputIn>>garbage>>garbage;
	InputIn>>NMCMD;

	InputIn>>garbage>>garbage;
	InputIn>>rdvel;

	bool comdens;

	InputIn>>garbage>>garbage;
	InputIn>>comdens;

	// cout<<"Compute Density: "<<comdens<<"\n";

	double lim1, lim2, dx;

	InputIn>>garbage>>garbage;
	InputIn>>lim1>>lim2>>dx;

	// cout<<"Range of x-values considered for computing densities: "<<lim1<<"\t"<<lim2<<endl;

	//Parameters for simulated annealing		
	bool simann;

	InputIn>>garbage>>garbage;
	InputIn>>simann;

	float TempB, TempE;

	InputIn>>garbage>>garbage;
	InputIn>>TempB>>TempE;

	int rampsize, rampstep;

	InputIn>>garbage>>garbage;
	InputIn>>rampsize>>rampstep;

	int ncycle;

	InputIn>>garbage>>garbage;
	InputIn>>ncycle;

	if(simann==1){	
		cout<<"Perform Simulated Annealing"<<endl;
		cout<<"Start Temperature: "<<TempB<<" K"<<endl;
		cout<<"End Temperature: "<<TempE<<" K"<<endl;
		cout<<"Delta T:"<<rampsize<<" K"<<endl;
		cout<<"Increase temperature after: "<<rampstep<<" MD steps"<<endl;
		cout<<"ncycles: "<<ncycle<<endl;
	}

	bool minimize;

	InputIn>>garbage>>garbage;
	InputIn>>minimize;                      //whether to minimize or not
	if(minimize)cout<<"****Perform minimization****"<<endl;

	InputIn.close();

	/////////////////////////////// ewald input file /////////

	ifstream EWALDIn("ewald.in",ios::in);
	if(!EWALDIn){
		cerr << "File EWALDIn could not be opened" <<endl;
		exit(1);
	}

	int Kvec[3], Grid[3], Order[3];

	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>Kvec[i];
	}

	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>Grid[i];
	}

	for (int i = 0; i < 3; i++){
		EWALDIn>>garbage>>garbage;
		EWALDIn>>Order[i];
	}                             

	EWALDIn.close();

	// cout<<"Initial positions read from the file: "<<posfile<<endl;

	//////////////////////get positions from the POSCAR file//////////////////////////////////////////////////////////

	ifstream PosIn(posfile.c_str(),ios::in);
	if(!PosIn){
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	getline(PosIn, garbage);
	istringstream StrStream(garbage);

	int n_atomtype=0;
	// Counting through the types of the atoms present in the cell and storing in the n_atomtype variable
	while(StrStream){
		getline(StrStream, garbage1, ' ');
		if(garbage1.compare("") != 0)
			n_atomtype = n_atomtype + 1;
	}
	
	string *atomtype;
	atomtype=new string [n_atomtype];

	istringstream strstream(garbage);

	tmp = 0;

	while(strstream){
		getline(strstream, garbage1, ' ');
		if(garbage1.compare("") != 0){
			atomtype[tmp]=garbage1;
			tmp = tmp + 1;
		}
	}
	
	int *natoms_type, natoms;  //number of atoms for each type and total atoms in the unit cell

	natoms_type=new int [n_atomtype];
	natoms=0;

	//Creating 3x3 boxcell array
	double **boxcell;
	boxcell=new double * [3];
	for(i = 0; i<3;i++){
		boxcell[i]=new double [3];
	}

	getline(PosIn, garbage); //read overall scaling factor

	//get box cell vectors
	for(i = 0; i<3;i++){
		for(j=0;j<3;j++){
			PosIn>>boxcell[i][j];
		}
	}
	getline(PosIn, garbage);//read the atoms again
	getline(PosIn, garbage1);

	//get number of atoms for each type
	for(i=0;i<n_atomtype;i++){
		PosIn>>natoms_type[i];
		natoms=natoms+natoms_type[i];
	}

	getline(PosIn, garbage);
	getline(PosIn, garbage);

	double *PosIons, *ForceIons, *vel;
	int *fixatoms;
	double *mass;

	mass=new double [natoms];
	PosIons=new double [natoms*3];
	ForceIons=new double [natoms*3];
	vel=new double [natoms*3];
	fixatoms=new int [natoms*3];

	//getting positions of each atom
	for(i=0;i<natoms;i++){
		PosIn>>PosIons[3*i]>>PosIons[3*i+1]>>PosIons[3*i+2];
	}

	PosIn.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// charge input file ////////////////////////////////////////////////////////////

	ifstream CHARGEIn("charge.in",ios::in);
	if(!CHARGEIn){
		cerr << "File CHARGEIn could not be opened" <<endl;
		exit(1);
	}

	float *chg;
	chg=new float [n_atomtype];

	for(int i=0; i<n_atomtype; i++){
		CHARGEIn>>chg[i];
	}

	CHARGEIn.close();

	/*Checking Charge Neutrality for the system*/
	double total_charge = 0;
	for (int i = 0; i < n_atomtype; i++){
		total_charge += chg[i]* natoms_type[i];
	}
	if(total_charge){
		cout<<"Error: System is not charge neutral"<<endl;
		return 0;
	}

	//creating the charge array for each atom present in the unit cell ////////////////////////////////////
	double *ion_charges;
	ion_charges=new double [natoms];
	int c=0;
	for (int i = 0; i < n_atomtype; i++){
		for (int j = 0; j < natoms_type[i]; j++){
			ion_charges[c]=chg[i];
			c++;
		}		
	}

	// Volume Calculations
    double A[3];
    double C[3]={boxcell[2][0],boxcell[2][1],boxcell[2][2]};
    crossProduct(boxcell[0],boxcell[1],A);
    volume = dotProduct(A,C);

	// Calculating the reciprocal vectors
    crossProduct(boxcell[1],boxcell[2],G[0]);
    crossProduct(boxcell[2],boxcell[0],G[1]);
    crossProduct(boxcell[0],boxcell[1],G[2]);
    for (int x = 0; x < 3; x++)
        for (int q = 0; q < 3; q++)
            G[x][q] /= volume;

	// print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");
	// print_coor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "COOR");
	// print_lammps_input_file(PosIons, chg, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "out.data");

	/*Check for orthogonality of sides*/
	/*Ewald method only works for unit cell with orthogonal sides*/
	if(dotProduct(boxcell[1],boxcell[0]) || dotProduct(boxcell[2],boxcell[0]) || dotProduct(boxcell[1],boxcell[2])){
		cout<<"Error: Unit Cell with Non Orthogonal Sides"<<endl;
		return 0;
	}

	double Lmin=min(boxcell[0][0],min(boxcell[1][1],boxcell[2][2]));
	double a=5.42/Lmin;
	double cutoff = Lmin/2;

	/*Useful configuration independent Computations for the reciprocal space summation*/
	/*Cofficients Bi[mi] of the bspline interpolation*/ //Refer to Essmann et al.
	CoeffX = new complex<double> [2*Kvec[0]+1];
    CoeffY = new complex<double> [2*Kvec[1]+1];
    CoeffZ = new complex<double> [2*Kvec[2]+1];
	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[0]+1)+i;}
        else {ic=i;}
        CoeffX[ic] = B(i,Order[0],Grid[0]);
    }
	for (int i = -Kvec[1]; i < Kvec[1]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[1]+1)+i;}
        else {ic=i;}
        CoeffY[ic] = B(i,Order[1],Grid[1]);
    }
	for (int i = -Kvec[2]; i < Kvec[2]+1; i++){
        int ic;
        if(i<0) {ic=(2*Kvec[2]+1)+i;}
        else {ic=i;}
        CoeffZ[ic] = B(i,Order[2],Grid[2]);
    }

	/* B(m1,m2,m3)*Exp(-|G|)/|G| term in the reciprocal loop*/
	double constant=(M_PI*M_PI)/(a*a);
	ExpFactor = new double [(2*Kvec[0]+1)*(2*Kvec[1]+1)*(2*Kvec[2]+1)];

	for (int i = -Kvec[0]; i < Kvec[0]+1; i++){
        for (int j = -Kvec[1]; j< Kvec[1]+1; j++){
            for (int k = -Kvec[2]; k < Kvec[2]+1; k++){
				int ii,jj,kk;
				if(i<0) ii=(2*Kvec[0]+1)+i;
                else ii=i;
                if(j<0) jj=(2*Kvec[1]+1)+j;
                else  jj=j;
                if(k<0) kk=(2*Kvec[2]+1)+k;
                else  kk=k;
				int temp=ii * ((2*Kvec[2]+1) * (2*Kvec[1]+1)) + jj * (2*Kvec[2]+1) + kk;
				double m[3];
                for (int t = 0; t < 3; t++){
                    m[t]=i*G[0][t]+j*G[1][t]+k*G[2][t];    
                }
                double m2=dotProduct(m,m);
				// ExpFactor[temp] = exp(-m2*constant)/m2;
				ExpFactor[temp] = norm(CoeffX[ii]*CoeffY[jj]*CoeffZ[kk])*exp(-m2*constant)/m2;
			}
		}
	}

	int Kvec2[3] = {6,6,6};
	PreExpFactor = new double [(2*Kvec2[0]+1)*(2*Kvec2[1]+1)*(2*Kvec2[2]+1)];
	for (int i = -Kvec2[0]; i < Kvec2[0]+1; i++){
		for (int j = -Kvec2[1]; j< Kvec2[1]+1; j++){
			for (int k = -Kvec2[2]; k < Kvec2[2]+1; k++){
				int ii,jj,kk;
				if(i<0) ii=(2*Kvec2[0]+1)+i;
				else ii=i;
				if(j<0) jj=(2*Kvec2[1]+1)+j;
				else  jj=j;
				if(k<0) kk=(2*Kvec2[2]+1)+k;
				else  kk=k;
				int temp=ii * ((2*Kvec2[2]+1) * (2*Kvec2[1]+1)) + jj * (2*Kvec2[2]+1) + kk;
				double m[3];
				for (int t = 0; t < 3; t++){
					m[t]=i*G[0][t]+j*G[1][t]+k*G[2][t];    
				}
				double m2=dotProduct(m,m);
				PreExpFactor[temp] = exp(-m2*constant)/m2;
			}
		}
	}
	
	double selfenergy=selfe(n_atomtype, natoms_type, chg, a)*unitzer;
	cout<<fixed<<setprecision(5)<<"Self Energy: "<<selfenergy<<" Kcal/mol"<<"\n\n";

	chrono::time_point<std::chrono::system_clock> start1, end1;
	start1 = chrono::system_clock::now();
	double recienergy=reci_energy(PosIons, ion_charges, natoms, a, boxcell, Kvec2)*unitzer;
	cout<<fixed<<setprecision(5)<<"Reciprocal Energy: "<<recienergy<<" Kcal/mol"<<"\n";
	end1 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds1 = end1- start1;
    time_t end_time1 = std::chrono::system_clock::to_time_t(end1);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds1.count() << " sec\n\n";

	chrono::time_point<std::chrono::system_clock> start2, end2;
	start2 = chrono::system_clock::now();
	double realenergy=real_energy(PosIons, ion_charges, natoms, a, boxcell,cutoff)*unitzer;
	cout<<fixed<<setprecision(5)<<"Real Energy: "<<realenergy<<" Kcal/mol"<<"\n";
	end2 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds2 = end2 - start2;
    time_t end_time2 = std::chrono::system_clock::to_time_t(end2);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds2.count() << " sec\n\n";

	chrono::time_point<std::chrono::system_clock> start3, end3;
	start3 = chrono::system_clock::now();
	double recienergy_bs=PM3DEwald(PosIons, ion_charges, natoms, a, boxcell, Grid, Kvec, Order)*unitzer;
	cout<<fixed<<setprecision(5)<<"Reciprocal Energy FFTW: "<<recienergy_bs<<" Kcal/mol"<<"\n";
	end3 = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds3 = end3 - start3;
    time_t end_time3 = std::chrono::system_clock::to_time_t(end3);
	cout<<fixed<<setprecision(8)<< "Elapsed time: " << elapsed_seconds3.count() << " sec\n\n";

	// chrono::time_point<std::chrono::system_clock> start4, end4;
	// start4 = chrono::system_clock::now();
	// double correctionTerm=dipoleCorrection(PosIons, ion_charges, natoms, boxcell)*unitzer;
	// // cout<<fixed<<setprecision(5)<<correctionTerm<<",";
	// cout<<fixed<<setprecision(5)<<"J(M,S): "<<correctionTerm<<" Kcal/mol"<<"\n";
	// end4 = chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds4 = end4 - start4;
    // time_t end_time4 = std::chrono::system_clock::to_time_t(end4);
	// // cout<< "Elapsed time: " << elapsed_seconds4.count() << " sec\n\n";

	// cout<<fixed<<setprecision(10)<<"Relative Error with FFTW: "<< error(recienergy,recienergy_bs)<<"\n";

/* using std::chrono::duration_cast; */
/* using HR = std::chrono::high_resolution_clock; */
/* using HRTimer = HR::time_point; */
/* using std::chrono::microseconds; */
/* using std::chrono::seconds; */

/*   HRTimer start = HR::now(); */
/*  recienergy_bs=bspline(PosIons, ion_charges, natoms, a, boxcell,60,6,5)*unitzer; */
/* 	cout<<fixed<<setprecision(5)<<"Reciprocal Energy FFTW: "<<recienergy_bs<<" Kcal/mol"<<"\n"; */
/*   HRTimer end = HR::now(); */
/*   auto duration = duration_cast<microseconds>(end - start).count(); */
/* 	cout<< "Elapsed time: " << duration << " usec"; */

	// delete dynamic variables 
	for(i=0;i<3;i++){
		delete [] boxcell[i]; 
	}

	delete [] PosIons;
	delete [] boxcell;
	delete [] atomtype;
	delete [] natoms_type;
	delete [] ForceIons;
	delete [] vel;
	delete [] fixatoms;
    delete [] chg;

	return 0;
}
//end of main loop
