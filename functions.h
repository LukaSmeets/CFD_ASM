#define sqr(a)      ((a)*(a))
#define max2(a,b)   (a > b ? a : b)
#define max3(a,b,c) (a > b ? max2(a,c) : max2(b,c))
#define min2(a,b)   (a < b ? a : b)
#define min3(a,b,c) (a < b ? min2(a,c) : min2(b,c))

extern int *int_1D_array(int np);
extern double *double_1D_array(int np);
extern double **double_2D_matrix(int nm, int np);
extern void free_double_2D_matrix (double **m, int nm);

extern void grid(void);
extern void gridbound(void);
extern void init(void);
extern void fixedbound(void);
int * fixedBoundary(void);

extern void bound(void);
extern void globcont(void);
extern void derivatives(void);

extern void ucoeff(double **aE,double **aW,double **aN,double **aS,double **aP,double **b,int Istart,int Iend,int Jstart,int Jend, int *fixed);
extern void vcoeff(double **aE,double **aW,double **aN,double **aS,double **aP,double **b,int Istart,int Iend,int Jstart,int Jend, int *fixed);
extern void pccoeff(double **aE,double **aW,double **aN,double **aS,double **aP,double **b,int Istart,int Iend,int Jstart,int Jend);
extern void Tcoeff(double **aE,double **aW,double **aN,double **aS,double **aP,double **b,int Istart,int Iend,int Jstart,int Jend);
extern void Alphacoeff(double **aE,double **aW,double **aN,double **aS,double **aP,double **b,int Istart,int Iend,int Jstart,int Jend, int *fixed);

extern void conv(void);
extern void solve   (double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend,int ITER,double ratio);
extern void solveGS (double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend,int ITER,double ratio);
extern void solveSOR_ADI(double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend);
extern void solveSOR(double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend);
extern void solveSIP(double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend);
extern void solveSIP_new(double **fi,double **b,double **aE,double **aW,double **aN,double **aS,double **aP,int Istart,int Iend,int Jstart,int Jend);
extern void velcorr(void);
extern void storeresults(void);

extern void density(void);
extern void viscosity(void);
extern void conductivity(void);

extern void output(void);
extern void print(void);
extern void print2(void);
extern void memalloc(void);

extern void printConv(void);

void outputmanager(int number);
void readconstants(void);
void defaultconstants(void);
void calcmean(void);
