//----------------------------------------------------------------
// V. Agarwal 
// Starting date: 21.06.19
// Last Modified:
// code for computing rates in a liquid
//-----------------------------------------------------------------

#include "libinclude.h" 
#include "potdec.h"
#include "mathdec.h"
#include "fundec.h"
#include "const.h"
#include "fftw3.h"
#include "complex"

int main(int argc, char **argv){

	// cout<<"\n******************************************************************"<<endl;
	// cout<<"******************************************************************"<<endl;
	// cout<<"      Code for adsorption/desorption near a sold surface          "<<endl;
	// cout<<"                    Written by Vishal Agarwal                     "<<endl;
	// cout<<"******************************************************************"<<endl;
	// cout<<"******************************************************************"<<endl<<endl<<endl;

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

	//////get the name of the file for input variables, default is input.in 
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

	cout<<"Temperature of the system:\t\t"<<Temp<<" K"<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>hyb;

	InputIn>>garbage>>garbage;
	InputIn>>NMCMD;

	InputIn>>garbage>>garbage;
	InputIn>>rdvel;

	bool comdens;

	InputIn>>garbage>>garbage;
	InputIn>>comdens;

	cout<<"Compute Density: "<<comdens<<"\n";

	double lim1, lim2, dx;

	InputIn>>garbage>>garbage;
	InputIn>>lim1>>lim2>>dx;

	cout<<"Range of x-values considered for computing densities: "<<lim1<<"\t"<<lim2<<endl;

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

	int gx, gy, gz, nd;

	EWALDIn>>garbage>>garbage;
	EWALDIn>>gx;

	EWALDIn>>garbage>>garbage;
	EWALDIn>>gy;                                    

	EWALDIn>>garbage>>garbage;
	EWALDIn>>gz;                             

	EWALDIn>>garbage>>garbage;
	EWALDIn>>nd;                             

	EWALDIn.close();

	cout<<"Initial positions read from the file: "<<posfile<<endl;

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
	
	int *natoms_type, natoms;  //number of atoms for each type

	natoms_type=new int [n_atomtype];
	natoms=0;

	//Creating 3x3 boxcell array
	float **boxcell;
	boxcell=new float * [3];
	for(i = 0; i<3;i++){
		boxcell[i]=new float [3];
	}

	getline(PosIn, garbage); //read overall scaling factor

	//get box cell vectors
	for(i = 0; i<3;i++){
		for(j=0;j<3;j++){
			PosIn>>boxcell[i][j];
		}
	}

	//get number of atoms for each type
	for(i=0;i<n_atomtype;i++){
		PosIn>>natoms_type[i];
		natoms=natoms+natoms_type[i];
	}

	getline(PosIn, garbage);
	getline(PosIn, garbage);

	double **PosIons, **ForceIons, **vel;
	int **fixatoms;
	float *mass;

	mass=new float [natoms];
	PosIons=new double * [natoms];
	ForceIons=new double *[natoms];
	vel=new double *[natoms];
	fixatoms=new int *[natoms];

	for(i=0;i<natoms;i++){
		PosIons[i]=new double [3];
		ForceIons[i]=new double [3];
		vel[i]=new double [3];		
		fixatoms[i]=new int [3];
	}

	//getting positions of each atom
	for(i=0;i<natoms;i++){
		PosIn>>PosIons[i][0]>>PosIons[i][1]>>PosIons[i][2];
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

	////////////////////////////////////////////////////////////

	// print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");
	float a=5.42/boxcell[0][0];
	cout<<selfe(n_atomtype, natoms_type, chg, a)<<"\n";
	
	// delete dynamic variables 

	for(i=0;i<3;i++){
		delete [] boxcell[i]; 
	}

	for(i=0;i<natoms;i++)
	{
		delete [] PosIons[i];
		delete [] ForceIons[i];
		delete [] vel[i];
		delete [] fixatoms[i];
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