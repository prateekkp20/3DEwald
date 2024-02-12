double lj_pot(double Dist, double *Par);
double morse_pot(double Dist, double *Par);

void lj_force(int atom1, int atom2, double Dist, double **PosIons, double **ForceIons, double *Par, float **boxcell);

void morse_force(int atom1, int atom2, double Dist, double **PosIons, double **ForceIons, double *Par, float **boxcell);


void fixatoms_forces(int natoms, double **ForceIons, int **fixatoms);

double selfinteraction(int natoms,int fullcharge,double betaa);

