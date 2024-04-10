// g++ Bspline_pm_copy_2.cpp -o Bspline_pm_copy_2 -lfftw3 -lm && ./Bspline_pm_copy_2
#include<cmath>
#include<vector>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<numeric>
#include<string>
#include<algorithm>
#include<iomanip>
#include<complex>
#include<fftw3.h>
#include <chrono>
#include <ctime>

#define ld long double
#define REAL 0
#define IMAG 1

ld epsilon=8.854187817e-12;
ld chargeelec=1.60217662e-19;               //charge of electron
ld Na=6.02214076e+23;                       //Avogadro Number
ld Ke=8.9875517923e+9;                      //Coulomb constant (1/4piepsilono)   unit N-m2C(-2)
ld rtoa=1e-10;                              //angstrom conversion
ld unitzer=(Na*Ke*chargeelec*chargeelec/rtoa)/1000/4.184; // Conversion factor =(1/4*pi*epsilon0)*(q1q2/r)  Dividing by 4.184e+03 making unit from kJ/mol to "kcal/mol"
const ld pi = acos(-1.0);

using namespace std;

ld M_n(ld u, long int n){
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
    // return 0;
}

struct Atoms {
    vector<vector<ld>> coordinates;
    vector<int> charges;
    ld Volume;
    vector<vector<ld>> lattice_vectors;
    vector<vector<ld>> reciprocal_lattice_vectors;
};

ld dotProduct(const vector<ld> v1,const vector<ld> v2) {
    ld result = 0.0;
    for (int i = 0; i < v1.size(); ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

vector<ld> crossProduct(const vector<ld> v_A,const vector<ld> v_B) {
    vector<ld> c_P(3, 0);
   c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
   c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
   c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
   return c_P;
}

Atoms read_file(string filename){

    Atoms atomdata;
    vector<vector<ld>> coordinates;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file." << endl;
        abort();
    }

    string line, line2;
    getline(file, line); // defining atoms names
    getline(file, line2);
    ld scale_factor = stod(line2); // scaling factor

    atomdata.lattice_vectors.resize(3, vector<ld>(3));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j){
            file >> atomdata.lattice_vectors[i][j];
            }

    int cations, anions;
    file >> cations >> anions;
    int total_atoms = cations + anions;
    atomdata.coordinates.resize(total_atoms, vector<ld>(3));
    string coordinate_type;
    file >> coordinate_type; // Read the coordinate type (either 'Direct' or 'Cartesian')

    if (coordinate_type == "Direct") {
        for (int i = 0; i < total_atoms; ++i){
            for (int j = 0; j < 3; ++j){
                file >> atomdata.coordinates[i][j];
                // cout<<atomdata.coordinates[i][j];
            }
        }
        // Convert from direct to Cartesian coordinates
        for (int i = 0; i < total_atoms; ++i)
            for (int j = 0; j < 3; ++j)
                atomdata.coordinates[i][j] *= scale_factor * atomdata.lattice_vectors[j][j];
    }

    else if (coordinate_type == "Cartesian") {
        for (int i = 0; i < total_atoms; ++i){
            for (int j = 0; j < 3; ++j){
                file >> atomdata.coordinates[i][j];
                }
            }
    }
    else {
        cerr << "Invalid coordinate type in file." << endl;
        abort();
    }

    vector<int> pos_charges(cations, 1);
    vector<int> neg_charges(anions, -1);
    atomdata.charges.reserve(total_atoms);
    atomdata.charges.insert(atomdata.charges.end(), pos_charges.begin(), pos_charges.end());
    atomdata.charges.insert(atomdata.charges.end(), neg_charges.begin(), neg_charges.end());

    atomdata.Volume=dotProduct(crossProduct(atomdata.lattice_vectors[1],atomdata.lattice_vectors[2]),atomdata.lattice_vectors[0]);
    atomdata.reciprocal_lattice_vectors.resize(3, vector<ld>(3));
    atomdata.reciprocal_lattice_vectors[0]=crossProduct(atomdata.lattice_vectors[1],atomdata.lattice_vectors[2]);
    atomdata.reciprocal_lattice_vectors[1]=crossProduct(atomdata.lattice_vectors[2],atomdata.lattice_vectors[0]);
    atomdata.reciprocal_lattice_vectors[2]=crossProduct(atomdata.lattice_vectors[0],atomdata.lattice_vectors[1]);
    vector<double> myarray;
    double myconstant=atomdata.Volume;
    for (int i = 0; i < 3; i++){
        transform(atomdata.reciprocal_lattice_vectors[i].begin(), atomdata.reciprocal_lattice_vectors[i].end(), atomdata.reciprocal_lattice_vectors[i].begin(), [&myconstant](auto& c){return c/myconstant;});
    }

    file.close();
    return atomdata;
}

void printAtomData(const Atoms& atomData) {
    cout << "Coordinates:\n";
    for (const auto& coord : atomData.coordinates) {
        for (const auto& val : coord) {
            cout << val << " ";
        }
        cout << "\n";
    }

    cout << "Charges:\n";
    for (const auto& charge : atomData.charges) {
        cout << charge << " ";
    }
    cout << "\n";
}

complex<ld>B(int m,int n,int K){
    const ld pi = acos(-1.0);
    const complex<ld> t(0.0, 1.0);
    complex<ld> bi_mi=exp((2*pi*m*(n-1))/K*t);
    complex<ld> denox;
    for (int f = 0; f < n-1; f++){
    denox+=M_n(f+1,n)*exp((2*pi*m*f)/K*t);
    }
    bi_mi/=denox;
    return bi_mi;
}

// complex<ld> B_FQ(int m1,int m2, int m3,vector<vector<ld>> u,vector<int> charges,int n_max, int n,vector<int> K){
//     complex<ld> B_FQ;

// 	fftw_complex *in;   //input variable using standard fftw syntax
// 	fftw_complex *out;	// output variable

//     ld two_pi_m1=2*pi*m1;
//     ld two_pi_m2=2*pi*m2;
//     ld two_pi_m3=2*pi*m3;
//     const complex<ld> t(0.0, 1.0);

// 	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K[0]*K[1]*K[2]);
// 	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K[0]*K[1]*K[2]);
// 	fftw_plan p;
// 	p = fftw_plan_dft_3d(K[0],K[1],K[2], in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//     vector<vector<ld>> x_direc;
//     x_direc.resize(charges.size(), vector<ld>(K[0]));
//     vector<vector<ld>> y_direc;
//     y_direc.resize(charges.size(), vector<ld>(K[1]));
//     vector<vector<ld>> z_direc;
//     z_direc.resize(charges.size(), vector<ld>(K[2]));

//     for (int i = 0; i < charges.size(); i++){
//         // for X direction
//         for (int k1 = 0;  k1 < K[0]; k1++){
//             x_direc[i][k1]=0;
//             for (int n1 = -n_max; n1 < n_max+1; n1++){
//                             x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K[0],n);
//             }
//         }

//         // for Y direction
//         for (int k2 = 0;  k2 < K[1]; k2++){
//             y_direc[i][k2]=0;
//             for (int n2 = -n_max; n2 < n_max+1; n2++){
//                             y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K[1],n);
//             }
//         }

//         // for Z direction
//         for (int k3 = 0;  k3 < K[2]; k3++){
//             z_direc[i][k3]=0;
//             for (int n3 = -n_max; n3 < n_max+1; n3++){
//                             z_direc[i][k3]+=M_n(u[i][2]-k3-n3*K[2],n);
//             }
//         }
//         // initializing the "in" vector with zero values
//     	for (int tx = 0; tx < K[0]; tx++){
// 		    for (int ty = 0; ty < K[1]; ty++){
// 			    for (int tz = 0; tz < K[2]; tz++){
// 				    in[tx * (K[2] * K[1]) + ty * K[2] + tz][0] = 0.0;
// 			    }
// 		    }
// 	    }

// 	    for (int j = 0; j < charges.size(); j++){
// 		    for (int tx = 0; tx < K[0]; tx++){
// 			    if (x_direc[j][tx] == 0)continue;

// 			    for (int ty = 0; ty < K[1]; ty++){
// 				    if (y_direc[j][ty] == 0)continue;

// 				    for (int tz = 0; tz < K[2]; tz++){
// 					    if (z_direc[j][tz] == 0)continue;

// 					    in[tx * (K[2] * K[1]) + ty * K[2] + tz][0] += charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
// 				    }
// 			    }
// 		    }
// 	    }
//     }

//     fftw_execute(p);
// 	fftw_destroy_plan(p);
// 	fftw_cleanup();
//     return B_FQ;
// }

ld pm_reciprocal_energy(Atoms atomdata){
    const ld pi = acos(-1.0);
    ld L1 = atomdata.lattice_vectors[0][0];
    ld L2 = atomdata.lattice_vectors[1][1];
    ld L3 = atomdata.lattice_vectors[2][2];
    ld beta=5.42/L1;
    vector<ld> L = { L1, L2, L3 };
    vector<int> K={60,60,60}; //Number of grid points in each direction
    int n=7; //order of b-spline interpolation
    int n_max=1;
    int total_atoms=atomdata.charges.size();

    // Calculating the fractional coordinates
    ld u[total_atoms][3];
    vector<int> M={1,1,1};
    // vector<int> M={6,6,6};
    for (int i = 0; i < total_atoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=K[j]*dotProduct(atomdata.coordinates[i],atomdata.reciprocal_lattice_vectors[j]);
        }
    }

    fftw_complex *in;   //input variable using standard fftw syntax
    fftw_complex *out;	// output variable

    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K[0]*K[1]*K[2]);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *K[0]*K[1]*K[2]);
    fftw_plan p;
    p = fftw_plan_dft_3d(K[0],K[1],K[2], in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    vector<vector<ld>> x_direc;
    x_direc.resize(atomdata.charges.size(), vector<ld>(K[0]));
    vector<vector<ld>> y_direc;
    y_direc.resize(atomdata.charges.size(), vector<ld>(K[1]));
    vector<vector<ld>> z_direc;
    z_direc.resize(atomdata.charges.size(), vector<ld>(K[2]));

    for (int i = 0; i < atomdata.charges.size(); i++){
        // for X direction
        for (int k1 = 0;  k1 < K[0]; k1++){
            x_direc[i][k1]=0;
            for (int n1 = -n_max; n1 < n_max+1; n1++){
                            x_direc[i][k1]+=M_n(u[i][0]-k1-n1*K[0],n);
            }
        }

        // for Y direction
        for (int k2 = 0;  k2 < K[1]; k2++){
            y_direc[i][k2]=0;
            for (int n2 = -n_max; n2 < n_max+1; n2++){
                            y_direc[i][k2]+=M_n(u[i][1]-k2-n2*K[1],n);
            }
        }

        // for Z direction
        for (int k3 = 0;  k3 < K[2]; k3++){
            z_direc[i][k3]=0;
            for (int n3 = -n_max; n3 < n_max+1; n3++){
                            z_direc[i][k3]+=M_n(u[i][2]-k3-n3*K[2],n);
            }
        }
    }
    // initializing the "in" vector with zero values
    for (int tx = 0; tx < K[0]; tx++){
        for (int ty = 0; ty < K[1]; ty++){
            for (int tz = 0; tz < K[2]; tz++){
                in[tx * (K[2] * K[1]) + ty * K[2] + tz][0] = 0.0;
            }
        }
    }

    for (int j = 0; j < total_atoms; j++){
    if (total_atoms == 0)
        continue;
        for (int tx = 0; tx < K[0]; tx++){
            if (x_direc[j][tx] == 0)continue;

            for (int ty = 0; ty < K[1]; ty++){
                if (y_direc[j][ty] == 0)continue;

                for (int tz = 0; tz < K[2]; tz++){
                    if (z_direc[j][tz] == 0)continue;

                    in[tx * (K[2] * K[1]) + ty * K[2] + tz][0] += atomdata.charges[j] * x_direc[j][tx] * y_direc[j][ty] * z_direc[j][tz];
                }
            }
        }
    }

    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup();

    ld energy = 0;
    ld constant=(pi*pi)/(beta*beta);
    int ii,jj,kk;
    for (int i = -M[0]; i < M[0]+1; i++){
        if(i<0) ii=M[0]+i;
        else  ii=i;
        for (int j = -M[1]; j < M[1]+1; j++){
            if(j<0) jj=M[1]+j;
            else  jj=j;
            for (int k = -M[2]; k < M[2]+1; k++){
                if(k<0) kk=M[2]+k;
                else  kk=k;
                if(i==0&&j==0&&k==0)continue;
                vector<ld> m = {i*atomdata.reciprocal_lattice_vectors[0][0]+j*atomdata.reciprocal_lattice_vectors[1][0]+k*atomdata.reciprocal_lattice_vectors[2][0],i*atomdata.reciprocal_lattice_vectors[0][1]+j*atomdata.reciprocal_lattice_vectors[1][1]+k*atomdata.reciprocal_lattice_vectors[2][1],i*atomdata.reciprocal_lattice_vectors[0][2]+j*atomdata.reciprocal_lattice_vectors[1][2]+k*atomdata.reciprocal_lattice_vectors[2][2]};
                ld m2=dotProduct(m,m);
                ld norm_FQ=out[ii * (K[2] * K[1]) + jj * K[2] + kk][REAL]*out[ii * (K[2] * K[1]) + jj * K[2] + kk][REAL]+out[ii * (K[2] * K[1]) + jj * K[2] + kk][IMAG]*out[ii * (K[2] * K[1]) + jj * K[2] + kk][IMAG];
                energy += norm_FQ*exp(-m2*constant)*norm(B(i,n,K[0])*B(j,n,K[1])*B(k,n,K[2]))/m2; //have to multiply the factor of F(Q)
                cout<<i<<" "<<j<<" "<<k<<" "<<norm_FQ<<"\n";
                // cout<<i<<" "<<j<<" "<<k<<" ("<<out[ii * (K[2] * K[1]) + jj * K[2] + kk][REAL]<<","<<out[ii * (K[2] * K[1]) + jj * K[2] + kk][IMAG] <<")\n";
                // cout<<i<<" "<<j<<" "<<k<<" ("<<out[ii * (K[2] * K[1]) + jj * K[2] + kk][REAL]<<","<<out[ii * (K[2] * K[1]) + jj * K[2] + kk][IMAG] <<")\n";
                // energy += exp((-pi2*m2)/(beta2))*norm(B(i,n,K[0])*B(j,n,K[1])*B(k,n,K[2]))/m2;
                // energy += exp((-pi2*m2)/(beta2))*norm(B_FQ(i,j,k,u,atomdata.charges,n_max,n,K))/m2;
            }
        }
    }
    energy*=unitzer;
    energy/=(2*pi*atomdata.Volume);
    return energy;
}

ld error(ld pm, ld ewald){
    return abs(pm-ewald)/ewald;
}

int main(){
    // clock_t start, end;
    // string in="/home/prateek/Documents/Prateek/3DEWALD/run/fifty/POSCAR.7";
    string in="/home/prateek/Documents/Prateek/3d_ewald/lampss_files/3D EWALD/random_generator/20atoms";
    Atoms atomData = read_file(in);
    // printAtomData(atomData);
    // start=clock();
    chrono::time_point<std::chrono::system_clock> start, end;
    start = chrono::system_clock::now();

    ld energy = pm_reciprocal_energy(atomData);

    end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout<<fixed<<setprecision(8)<<"Elapsed time: " <<elapsed_seconds.count()<<" sec\n";
    // end=clock();
    cout <<fixed<<setprecision(5)<<"Energy: "<<energy<<" kcal/mol"<<"\n" ;
    cout << fixed << setprecision(15) <<"Error: "<<error(energy,5.330360132800516e+02)<<"\n";
    // ld time_taken= double(end-start)/double(CLOCKS_PER_SEC);
    // cout << fixed << setprecision(10)<<"Time Taken: "<<time_taken<<" Sec"<<"\n" ;
    return 0;
}

// int main(){
//     ofstream filer("Bspline4.txt");
//     if (!filer.is_open()) {
//         cerr << "Error opening file!" << endl;
//         return 1; // Return an error code
//     }
//     vector<ld> EWALD={1.404731514289580e+03,8.354139318895441e+02,6.662950659666857e+02,4.273425704945629e+02,2.773380683582925e+02,3.541768938144129e+02};
//     for (int interator = 0; interator < 6; interator+=1){
//         clock_t start, end;
//         cout<<"go"<<interator<<"\n";
//         string in="3d_ewald\\lampss_files\\3D EWALD\\random_generator\\cell";
//         Atoms atomData = read_file(in+to_string(10+5*interator));
//         // filer <<"cell"<<10+5*interator<<"\n";
//         start=clock();
//         ld energy = pm_reciprocal_energy(atomData);
//         end=clock();
//         filer << fixed << setprecision(8) <<"Energy: "<<energy<<" kcal/mol"<<"\n";
//         // filer <<setprecision(15)<<error(energy,EWALD[interator])<<"\n";*
//         filer <<"Error: "<<setprecision(15)<<error(energy,EWALD[interator])<<"\n";
//         ld time_taken = double(end-start)/double(CLOCKS_PER_SEC);
//         // filer << fixed <<setprecision(2)<<time_taken<<endl;
//         filer << fixed <<"Time Taken: "<<setprecision(2)<<time_taken<<" Sec\n"<<endl;
//         cout<<"done"<<interator<<"\n";
//     }
//     filer.close();
//     return 0;
// }