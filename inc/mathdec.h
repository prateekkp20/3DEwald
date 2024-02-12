
double getnorm(double *val, int len);
double getmax(double *val, int len);
double getabsmax(double *val, int len);

double getMatmax(double **mat, int rows, int columns);
double getMatabsmax(double **mat, int rows, int columns);
double getMatRMS(double **mat, int rows, int columns);

double factorial(int num);
////////////////////////////////////////////////////////////////////
double vec_dot(double *vec1, double *vec2, int len);

int vec_dot_mat(double *vec, int veclen, double **mat, int rows, int columns, double *vecres);


void  vec_X_vecT(double *vec1, double *vec2, int len1, int len2, double **matres );

int mat_X_vec(double *vec, int veclen, double **mat, int rows, int columns, double *vecres);

void mat_X_mat(double **mat1, int rows1, int columns1, double **mat2, int rows2, int columns2, double **matres);

void matT_X_mat(double **mat1, int rows1, int columns1, double **mat2, int rows2, int columns2, double **matres);

//////////////////////////////////////////////////////////////////////////////////////
void mat_transpose(double **mat, int rows, int columns, double **mat_trans);

void vec_equate(double *vec1, double *vec2, int len);

void mat_equate(double **mat1, double **mat2, int rows, int columns);

void invmat(double **mat, int rows, double **invmat);

void eye(int rows, int columns, double **mat);

void zeros(int rows, int columns, double **mat);


void zeros_int(int rows, int columns, int **mat);
double sum2mat(double **mat, int rows, int columns);


double DotVecPos(double **PosIons, int atom1, int atom2);
/////////////////////////////////////////////////////////
void print_vec(double *vec, int len, string str);

void print_mat(double **mat, int rows, int columns, string str);

void print_mat_int(int **mat, int rows, int columns, string str);

void print_pos_forces(double *Pos, double *Jac, double len);


void print_prog(double maxf, double normx, double delPot, float maxF, float delPos, float delV );


///////////////////////////////////////////////////
double dist(double **PosIons, int atom1, int atom2, float **boxcell);
void ini_vel(double **vel, int natoms, float Temp, float *mass, int numrand, int **fixatoms);
void ini_mass(float *mass, int natoms, string *type);


void readvel(double **vel, int natoms);

void printVel(double **vel, int natoms);

void getCOM(double **PosIons, int natoms, int natom1, int natom2, float **boxcell, float *mass, double *RCOM);

void printprobVel(double **vel, int natoms, float *mass, float Temp);

double enfft(double **PosIons, float **boxcell, int natoms, int *fullcharge, double betaa, int nd, int gx, int gy, int gz, double hx, double hy, int *gxnew, int *gznew, double wx, double wy, double wz);


