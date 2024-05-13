// g++ Bspline_pm.cpp -o Bspline_pm -lfftw3 -lm -fopenmp && ./Bspline_pm
#include<cmath>
#include<vector>
#include<stdio.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include<algorithm>
#include <cmath>
#include <iomanip>
#include <complex>
#include "omp.h"
#include <thread>

#define ld long double
#define REAL 0
#define IMAG 1

ld epsilon=8.854187817e-12;
ld chargeelec=1.60217662e-19;               //charge of electron
ld Na=6.02214076e+23;                       //Avogadro Number
ld Ke=8.9875517923e+9;                      //Coulomb constant (1/4piepsilono)   unit N-m2C(-2)
ld rtoa=1e-10;                              //angstrom conversion
// Conversion factor =(1/4*pi*epsilon0)*(q1q2/r)  Dividing by 4.184e+03 making unit from kJ/mol to "kcal/mol"
ld unitzer=(Na*Ke*chargeelec*chargeelec/rtoa)/1000/4.184;

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

    if (coordinate_type == "Direct"){
        for (int i = 0; i < total_atoms; ++i){
            for (int j = 0; j < 3; ++j){
                file >> atomdata.coordinates[i][j];
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

complex<ld> B_FQ(int m1,int m2, int m3,vector<vector<ld>> u,vector<int> charges,int n_max, int n,vector<int> K){
    complex<ld> B_FQ;
    const ld pi = acos(-1.0);
    ld two_pi_m1=2*pi*m1;
    ld two_pi_m2=2*pi*m2;
    ld two_pi_m3=2*pi*m3;
    const complex<ld> t(0.0, 1.0);
    
    for (int i = 0; i < charges.size(); i++){
        // for X direction
        complex<ld> x_dir;
        for (int k1 = 0;  k1 < K[0]; k1++){
            complex<ld> xyz;
            for (int n1 = -n_max; n1 < n_max+1; n1++){
                            xyz+=M_n(u[i][0]-k1-n1*K[0],n);
            }
            xyz*=exp((two_pi_m1*k1)/K[0]*t);
            x_dir+=xyz;
        }

        // for Y direction
        complex<ld> y_dir;
        for (int k2 = 0;  k2 < K[1]; k2++){
            complex<ld> xyz;
            for (int n2 = -n_max; n2 < n_max+1; n2++){
                            xyz+=M_n(u[i][1]-k2-n2*K[1],n);
            }
            xyz*=exp((two_pi_m2*k2)/K[1]*t);
            y_dir+=xyz;
        }

        // for Z direction
        complex<ld> z_dir;
        for (int k3 = 0;  k3 < K[2]; k3++){
            complex<ld> xyz;
            for (int n3 = -n_max; n3 < n_max+1; n3++){
                            xyz+=M_n(u[i][2]-k3-n3*K[2],n);
            }
            xyz*=exp((two_pi_m3*k3)/K[2]*t);
            z_dir+=xyz;
        }
        
        const complex<ld> charge(charges[i], 0.0);
        B_FQ+=x_dir*y_dir*z_dir*charge;
    }
    return B_FQ;
}

ld pm_reciprocal_energy(Atoms atomdata){
    const ld pi = acos(-1.0);
    ld L1 = atomdata.lattice_vectors[0][0];
    ld L2 = atomdata.lattice_vectors[1][1];
    ld L3 = atomdata.lattice_vectors[2][2];
    ld beta=5.42/L1;
    vector<ld> L = { L1, L2, L3 };
    vector<int> K={20,20,20}; //Number of grid points in each direction
    int n=7; //order of b-spline interpolation
    int n_max=1;
    int total_atoms=atomdata.charges.size();
    vector<vector<ld>> u;
    u.resize(total_atoms, vector<ld>(3));
    // vector<int> M={1,1,1};
    vector<int> M={6,6,6};
    for (int i = 0; i < total_atoms; i++){
        for (int j = 0; j < 3; j++){
            u[i][j]=K[j]*dotProduct(atomdata.coordinates[i],atomdata.reciprocal_lattice_vectors[j]);
        }
    }
    ld energy = 0;
    ld pi2=pi*pi;
    ld beta2=beta*beta;
    omp_set_num_threads(thread::hardware_concurrency());
    // #pragma omp parallel for schedule(runtime) reduction(+: energy) collapse(3)
    // #pragma omp parallel for schedule(runtime) reduction(+: energy) 
    for (int i = -M[0]; i < M[0]+1; i++){
        for (int j = -M[1]; j < M[1]+1; j++){
            for (int k = -M[2]; k < M[2]+1; k++){
                if(i==0&&j==0&&k==0)continue;
                vector<ld> m = {i*atomdata.reciprocal_lattice_vectors[0][0]+j*atomdata.reciprocal_lattice_vectors[1][0]+k*atomdata.reciprocal_lattice_vectors[2][0],i*atomdata.reciprocal_lattice_vectors[0][1]+j*atomdata.reciprocal_lattice_vectors[1][1]+k*atomdata.reciprocal_lattice_vectors[2][1],i*atomdata.reciprocal_lattice_vectors[0][2]+j*atomdata.reciprocal_lattice_vectors[1][2]+k*atomdata.reciprocal_lattice_vectors[2][2]};
                ld m2=dotProduct(m,m);
                energy += exp((-pi2*m2)/(beta2))*norm(B(i,n,K[0])*B(j,n,K[1])*B(k,n,K[2]))*norm(B_FQ(i,j,k,u,atomdata.charges,n_max,n,K))/m2;
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
    clock_t start, end;
    string in="/home/prateek/Documents/Prateek/3DEwald/run/bench/bench1.POSCAR";
    // string in="/home/prateek/Documents/Prateek/3d_ewald/lampss_files/3D EWALD/random_generator/20atoms";
    Atoms atomData = read_file(in);
    // printAtomData(atomData);
    start=clock();
    ld energy = pm_reciprocal_energy(atomData);
    end=clock();
    ld time_taken= double(end-start)/double(CLOCKS_PER_SEC);
    cout << fixed << setprecision(2)<<"Elapsed Taken: "<<time_taken<<" Sec"<<"\n" ;
    cout << fixed << setprecision(5) <<"Energy: "<<energy<<" kcal/mol"<<"\n" ;
    // cout << fixed << setprecision(15) <<"Error: "<<error(energy,5.330360132800516e+02)<<"\n";
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
//             string in="/home/prateek/Documents/Prateek/3d_ewald/lampss_files/3D EWALD/random_generator/cell";
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