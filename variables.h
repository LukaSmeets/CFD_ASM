double *x;               /* x coordinate on pressure points [m] */
double *x_u;             /* x coordinate on u-velocity points [m] */
double *y;               /* y coordinate on pressure points [m] */
double *y_v;             /* y coordinate on v-velocity points [m] */
double Dt;
int    NDt;

double **u;
double **v;
double **pc;
double **p;
double **T;
double **rho;
double **mu;
double **Gamma;
double **Cp;
double **Alpha;

double **u_old;
double **v_old;
double **pc_old;
double **T_old;
double **rho_old;
double **Alpha_old;

double **usum;
double **vsum;
double **psum;
double **rhosum;
double **musum;
double **alphasum;

double **umean;
double **vmean;
double **pmean;
double **rhomean;
double **mumean;
double **alphamean;

double **dudx;
double **dudy;
double **dvdx;
double **dvdy;

double **aE;
double **aW;
double **aN;
double **aS;
double **aP;
double **b;

int    *Istart;
int    *Iend;
int    *Jstart;
int    *Jend;

int    Jinmin;
int    Jinmax;
int    Joutmin;
int    Joutmax;
int    Isolmin;
int    Isolmax;
int    Jsolmin;
int    Jsolmax;
int    Ipref;
int    Jpref;

double SAVG;
double SMAX;
double DELTA;

double **SP;
double **Su;

double **F_u;
double **F_v;

double **d_u;
double **d_v;

double m_in;
double m_out;

double *relax;
int *solve_fi;
int *print_fi;

double omega;

char outputdir[200];
double IinLeft;
double IinRight;

int NPI;//        100        /* number of grid cells in x-direction [-] */
int NPJ;//        50        /* number of grid cells in y-direction [-] */

int Jfix11;
int Jfix12;
int Ifix11;
int Ifix12;
int Ifix21;
int Ifix22;

int IJfix[6];
int *fixed;

double XMAX;//       1.5      /* width of the domain [m] */ 
double YMAX;//       0.5      /* height of the domain [m] */
double RHOL;//       1000.    /* rho liquid [kg/m^3] */
double RHOG;//       1.2      /* rho gas [kg/m^3] */
double USLIP;//      0.2      
double ALPHAIN;//    0.0417
int DUMPFREQ;//   100
double TOTAL_TIME;
int TOTAL_DUMPS;
double current_time;
int NNOZ;
double DNOZ;
double V_IN;
int CONFIG;
