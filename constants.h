//#define NPI        100        /* number of grid cells in x-direction [-] */
//#define NPJ        50        /* number of grid cells in y-direction [-] */
//#define XMAX       1.5      /* width of the domain [m] */ 
//#define YMAX       0.5      /* height of the domain [m] */
#define PI         3.1415927 /* value of pi [-] */
#define MAX_ITER   1      /* maximum number of outer iterations [-] */
#define U_ITER     5        /* number of Newton iterations for u equation [-] */
#define V_ITER     5        /* number of Newton iterations for u equation [-] */
#define PC_ITER    30       /* number of Newton iterations for pc equation [-] */
#define T_ITER     5        /* number of Newton iterations for T equation [-] */
#define ALPHA_ITER 5        /* number of Newton iterations for T equation [-] */
#define SMAXneeded 1E-6      /* maximum accepted error in mass balance [kg/s] */
#define SAVGneeded 1E-6      /* maximum accepted average error in mass balance [kg/s] */
#define LARGE      1E30      /* arbitrary very large value [-] */
#define SMALL      1E-30     /* arbitrary very small value [-] */
#define P_ATM      101000.   /* athmospheric pressure [Pa] */
#define U_IN       0.08      /* in flow velocity [m/s] */
//#define V_IN       1.        /* in flow velocity [m/s] */
#define Cdrag      50.E3
//#define RHOL       1000.
//#define RHOG       1.2
//#define USLIP      0.2
//#define ALPHAIN    0.0417
#define CZERO      0.1
#define BIT        1
#define CMUBIT     0.6
#define DBUB       0.004
//#define DUMPFREQ   100

#define NU         1
#define NV         2
#define NPC        3
#define NP         4
#define NT         5
#define NRHO       6
#define NMU        7
#define NGAMMA     8
#define NALPHA     9
#define NFIMAX     9

#define TRUE       1
#define FALSE      0
#define NPRINT     10
