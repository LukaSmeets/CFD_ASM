/***** Solves: Unsteady, compressible convection-diffusion problems.

****** Description:
****** This program solves unsteady convection-diffusion problems	
****** using the transient simple algorithm described in ch. 8.7.1 in "Computational 
****** Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
****** variables follow exactly the notations in this reference, and all 
****** equations cited are from this reference unless mentioned otherwise.

****** References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. 
******			    Malalasekera, Longman Group Ltd, 1995
******/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables.h"
#include "constants.h"
#include "functions.h"

/* ################################################################# */
int main(int argc, char *argv[])
/* ################################################################# */
{
	int    iter;
	int	  *fixed;
	double du, dv;
	FILE *fp;
	char logfile[200];// = "logfile.dat";

	init();
	fixed = fixedBoundary();
	sprintf (logfile, "%slogfile.dat", outputdir);

	fixedbound(); /* set boundary values that remain untouched during */
	              /* the iteration process */

	if ((fp = fopen(logfile,"w"))==NULL) {
		printf("could not open file: %s \n", logfile);
		exit(1);
	}
	while ((fp = fopen(logfile,"w"))==NULL);

	fprintf (fp, "#Time\t\tu[%d][%d]\tv[%d][%d]\tAlpha[%d][%d]\n", NPI/2, NPJ/2, NPI/2, NPJ/2, NPI/2, NPJ/2);
	fprintf (fp, "%11.4e\t%12.4e\t%12.4e\t%12.4e\n", current_time, u[NPI/2][NPJ/2], v[NPI/2][NPJ/2], Alpha[NPI/2][NPJ/2]);
	fclose(fp);

	outputmanager(NDt); /* dump initial conditions */
	
	for (current_time = Dt; current_time <= TOTAL_TIME; current_time += Dt) {
		iter = 1;
		NDt++;
		while (iter <= MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded) { /* outer iteration loop */

			if (solve_fi[NU] == TRUE || solve_fi[NV] == TRUE ) derivatives();

			if (solve_fi[NU] == TRUE) {
				ucoeff(aE, aW, aN, aS, aP, b, Istart[NU], Iend[NU], Jstart[NU], Jend[NU], fixed);
				solveGS(u, b, aE, aW, aN, aS, aP, Istart[NU], Iend[NU], Jstart[NU], Jend[NU], U_ITER, 0.25);
			} /* if */

			if (solve_fi[NV] == TRUE) {
				vcoeff(aE, aW, aN, aS, aP, b, Istart[NV], Iend[NV], Jstart[NV], Jend[NV], fixed);
				solveGS(v, b, aE, aW, aN, aS, aP, Istart[NV], Iend[NV], Jstart[NV], Jend[NV], V_ITER, 0.25);
			} /* if */

			if (solve_fi[NPC] == TRUE) {
				bound();
				pccoeff(aE, aW, aN, aS, aP, b, Istart[NPC], Iend[NPC], Jstart[NPC], Jend[NPC]);
				solve(pc, b, aE, aW, aN, aS, aP, Istart[NPC], Iend[NPC], Jstart[NPC], Jend[NPC], PC_ITER, 0.1);
			} /* if */

			velcorr(); /* Correct pressure and velocity */

			if (solve_fi[NT] == TRUE) {
				Tcoeff(aE, aW, aN, aS, aP, b, Istart[NT], Iend[NT], Jstart[NT], Jend[NT]);
				solveGS(T, b, aE, aW, aN, aS, aP, Istart[NT], Iend[NT], Jstart[NT], Jend[NT], T_ITER, 0.25);
			} /* if */

			if (solve_fi[NALPHA] == TRUE) {
				Alphacoeff(aE, aW, aN, aS, aP, b, Istart[NALPHA], Iend[NALPHA], Jstart[NALPHA], Jend[NALPHA], fixed);
				solveGS(Alpha, b, aE, aW, aN, aS, aP, Istart[NALPHA], Iend[NALPHA], Jstart[NALPHA], Jend[NALPHA], ALPHA_ITER, 0.1);
			} /* if */

			if (solve_fi[NRHO]   == TRUE) density();

			if (solve_fi[NMU]    == TRUE) viscosity();

			if (solve_fi[NGAMMA] == TRUE) conductivity();
		
			if (current_time == Dt) {
				printf ("Time\t\tu[%d][%d]\tv[%d][%d]\tAlpha[%d][%d]\n", NPI/2, NPJ/2, NPI/2, NPJ/2, NPI/2, NPJ/2);
			} /* if */
			if (NDt % NPRINT == 0) {
				du = d_u[NPI/2][NPJ/2]*(pc[NPI/2-1][NPJ/2] - pc[NPI/2][NPJ/2]);
				dv = d_v[NPI/2][NPJ/2]*(pc[NPI/2][NPJ/2-1] - pc[NPI/2][NPJ/2]);
				printf ("%11.4e\t%12.4e\t%12.4e\t%12.4e\n", current_time, u[NPI/2][NPJ/2], v[NPI/2][NPJ/2], Alpha[NPI/2][NPJ/2]);
			} /* if */

//			if ((fp = fopen(logfile,"w"))==NULL) {
//				printf("could not open file: %s \n", logfile);
//				exit(1);
//			}
			while ((fp = fopen(logfile,"a"))==NULL);
			fprintf (fp, "%11.4e\t%12.4e\t%12.4e\t%12.4e\n", current_time, u[NPI/2][NPJ/2], v[NPI/2][NPJ/2], Alpha[NPI/2][NPJ/2]);
			fclose(fp);
			
			bound();
			storeresults(); /* Store data at current time level in arrays for "old" data*/
			calcmean(); /* Calculate time averaged data fields */

			iter++;
		} /* for outer iteration loop */

		/* reset SMAX and SAVG */
		SMAX = LARGE;
		SAVG = LARGE;

		/* create dumpfile if necessary */
		if(NDt%DUMPFREQ == 0) outputmanager(NDt/DUMPFREQ);
//		if(((int)TOTAL_TIME/Dt)%NDt == DUMPFREQ) outputmanager(NDt/DUMPFREQ);
	} /* for Dt */
//	output();

	return 0;

} /* main */


/* ################################################################# */
void grid(void)
/* ################################################################# */
{
/***** Purpose: Defining the geometrical variables ******/
/*****          See fig. 6.2-6.4 in ref. 1 ******/
	int    I, J, i, j;
	double Dx, Dy;

	/* Length of volume element */

	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;
	DELTA = sqrt(Dx*Dy);

	/* Length variable for the scalar points in the x direction */

	x[0] = 0.;
	x[1] = 0.5*Dx;

	for (I = 2; I <= NPI; I++)
		x[I] = x[I-1] + Dx;

	x[NPI+1] = x[NPI] + 0.5*Dx;

	/* Length variable for the scalar points fi[i][j] in the y direction */

	y[0] = 0.;
	y[1] = 0.5*Dy;

	for (J = 2; J <= NPJ; J++)
		y[J] = y[J-1] + Dy;

	y[NPJ+1] = y[NPJ] + 0.5*Dy;

	/* Length variable for the velocity components u[i][j] in the x direction */

	x_u[0] = 0.;
	x_u[1] = 0.;

	for (i = 2; i <= NPI + 1; i++)
		x_u[i] = x_u[i-1] + Dx;

	/* Length variable for the velocity components v[i][j] in the y direction */

	y_v[0] = 0.;
	y_v[1] = 0.;
	for (j = 2; j <= NPJ + 1; j++)
		y_v[j] = y_v[j-1] + Dy;

} /* grid */


/* ################################################################# */
void gridbound(void)
/* ################################################################# */
{
/***** Purpose: Definding the geometrical boundary coordinates ******/

	Jinmin = 1;
	Jinmax = NPJ;
	Joutmin = 1;
	Joutmax = NPJ;

} /* gridbound */

/* ################################################################# */
int * fixedBoundary(void)
/* ################################################################# */
{
/***** Purpose: Defining the boundaries ******/
	//int IJfix[6];

	if (CONFIG == 1) {
		Jfix11 = 0;
		Jfix12 = 100;
		Ifix11 = 0;
		Ifix12 = 11;
		Ifix21 = 54;
		Ifix22 = 65;
	}
	
	if (CONFIG == 2) {
	    Jfix11 = 50;
		Jfix12 = 100;
		Ifix11 = 0;
		Ifix12 = 20;
		Ifix21 = 44;
		Ifix22 = 65;
	}

	if (CONFIG == 3) {
		Jfix11 = 75;
		Jfix12 = 150;
		Ifix11 = 0;
		Ifix12 = 15;
		Ifix21 = 49;
		Ifix22 = 65;
	}

	IJfix[0] = Jfix11;
	IJfix[1] = Jfix12;
	IJfix[2] = Ifix11;
	IJfix[3] = Ifix12;
	IJfix[4] = Ifix21;
	IJfix[5] = Ifix22;

	return IJfix;
} /* fixedBoundary */


/* ################################################################# */
void init(void)
/* ################################################################# */
{
/***** Purpose: To initilise all parameters. ******/
	int    I, J, i, j, NFI;
	int    inputfile_used = 1;

	if (inputfile_used == 1) {
		readconstants();
	} else {
		defaultconstants();
		}

	memalloc();
	grid();
	gridbound(); /* define the coordinates at the boundary nodes */

	/* Initialising all variables  */

	omega = 1.0; /* Over-relaxation factor for SOR solver */

	/* Initialize convergence parameters at large values */

	SMAX = LARGE;
	SAVG = LARGE;

	m_in  = 1.;
	m_out = 1.;

	/* calculate dumpfrequency [1/timesteps] */

	DUMPFREQ = (int) TOTAL_TIME/(Dt*TOTAL_DUMPS);
	if (DUMPFREQ == 0) DUMPFREQ = 1; /* minimum dumpfrequency is equal to 1 */

	for (I = 0; I <= NPI + 1; I++) {
		i = I;
		for (J = 0; J <= NPJ + 1; J++) {
			j = J;
			u      [i][J] = 0.;       /* Velocity in x-direction */
			v      [I][j] = 0.;       /* Velocity in y-direction */
			T      [I][J] = 273.;     /* Temperature */
			Alpha  [I][J] = 0.;       /* Gas fraction */
			rho    [I][J] = (1. - Alpha[I][J])*RHOL + Alpha[I][J]*RHOG;      /* Density */
//			p      [I][J] = rho[I][J]*9.81*(NPI-I)/NPI*XMAX;       /* Relative pressure */
			p      [I][J] = rho[I][J]*9.81*(NPJ-J)/NPJ*YMAX;       /* Relative pressure */
			pc     [I][J] = 0.;       /* Pressure correction (equivalent to p´ in ref. 1). */
			mu     [I][J] = 1.E-3;    /* Viscosity */
			Cp     [I][J] = 1013.;    /* J/(K*kg) Heat capacity - assumed constant for this problem */
			Gamma  [I][J] = 0.0315/Cp[I][J]; /* Thermal conductivity */
			d_u    [i][J] = 0.;       /* Variable d[i][J] to calculate pc defined in 6.23 */
			d_v    [I][j] = 0.;       /* Variable d[I][j] to calculate pc defined in 6.23 */
			b      [I][J] = 0.;	  /* The general constant */
			SP     [I][J] = 0.;       /* Source term */
			Su     [I][J] = 0.;	  /* Source term */
			u_old  [i][J] = u[i][J];  /* Velocity in x-direction old timestep */
			v_old  [I][j] = v[I][j];  /* Velocity in y-direction old timestep */
			pc_old [I][J] = pc[I][J]; /* Pressure correction old timestep */
			T_old  [I][J] = T[I][J];  /* Temperature old timestep */
			rho_old[I][J] = rho[I][J];/* Density old timestep */
			Alpha_old[I][J] = 0.;       /* Gas fraction old timestep */
		} /* for J */
	} /* for I */

	
	if (solve_fi[NU] == TRUE)
		for (J = Joutmin; J <= Joutmax; J++) /* Important to avoid crash!! */
//			u[NPI][J] = U_IN/1000;
			u[NPI][J] = 0.;
			/* Othervise m_out calculated in subroutine globcont */
			/* would be zero at first iteration=>m_in/m_out =INF */

	/* Initialising the logical parameter for which variable fi to solve */
	/* and to print results for  */

	for (NFI = 1; NFI <= NFIMAX; NFI++) {
		solve_fi[NFI] = FALSE; 
		print_fi[NFI] = TRUE;
	} /* for */
	
	solve_fi[NU]     = TRUE;
	solve_fi[NV]     = TRUE;
	solve_fi[NPC]    = TRUE;
	solve_fi[NALPHA] = TRUE;
	solve_fi[NRHO]   = TRUE;
	solve_fi[NMU ]   = TRUE;

	/* Setting the relaxation parameters */

	relax[NU]   = 0.8;             /* See eq. 6.36 */
	relax[NV]   = relax[NU];       /* See eq. 6.37 */
	relax[NPC]  = 1.1 - relax[NU]; /* See eq. 6.33 */
	relax[NT]   = 1.0;  /* Relaxation factor for temperature */
	relax[NRHO] = 0.1;  /* Relaxation factor for density */
	relax[NALPHA] = 0.25;  /* Relaxation factor for gas volume fraction */


//	relax[NU]   = 0.8;             /* See eq. 6.36 */
//	relax[NV]   = relax[NU];       /* See eq. 6.37 */
//	relax[NPC]  = 1.1 - relax[NU]; /* See eq. 6.33 */
//	relax[NT]   = 1.0;  /* Relaxation factor for temperature */
//	relax[NRHO] = 0.1;  /* Relaxation factor for density */
//	relax[NALPHA] = 0.25;  /* Relaxation factor for gas volume fraction */

	/* Definition of first internal node for the different fi variables. */
	/* See fig. 9.1. */

	Istart[NU] = 2;
	Iend  [NU] = NPI;
	Jstart[NU] = 1;
	Jend  [NU] = NPJ;

	Istart[NV] = 1;
	Iend  [NV] = NPI;
	Jstart[NV] = 2;
	Jend  [NV] = NPJ;

	Istart[NPC] = 1;
	Iend  [NPC] = NPI;
	Jstart[NPC] = 1;
	Jend  [NPC] = NPJ;

	Istart[NT] = 1;
	Iend  [NT] = NPI;
	Jstart[NT] = 1;
	Jend  [NT] = NPJ;

	Istart[NALPHA] = 1;
	Iend  [NALPHA] = NPI;
	Jstart[NALPHA] = 1;
	Jend  [NALPHA] = NPJ;

} /* init */

/* ################################################################# */
void fixedbound(void)
/* ################################################################# */
{
/***** Purpose: Specify the boundary values that remain untouched during *****/
/*****          the iteration process ******/
	int I, J, n;
	double distance;


//	for (J = 0; J <= NPJ + 1; J++) {
//		u[1][J] = U_IN*1.5*(1.-sqr(2.*(y[J]-YMAX/2.)/YMAX)); /* inlet - EDIT: SET INLET X VELOCITY (05-04-2019) */
//	} /* for J */

//	for (I = 0; I <= NPI + 1; I++) {
//		v[I][1] = V_IN; //*1.5*(1.-sqr(2.*(y[I]-YMAX/2.)/YMAX)); /* inlet - EDIT: SET INLET Y VELOCITY (05-04-2019) */
//	} /* for I */

	for (J = 0; J <= NPJ + 1; J++) {
		u[0][J]   = 0.; /* inlet */
		u[NPI][J] = 0.; /* outlet */
	} /* for J */

	for (I = 0; I <= NPI + 1; I++) {
		v[I][0]     = 0.; /* bottom wall */
		v[I][NPJ+1] = 0.; /* top wall */
	} /* for I */

//	for (J = JinLeft*NPJ; J <= JinRight*NPJ; J++) {
//		Alpha[0][J] = ALPHAIN; /* inlet */
//	} /* for J */
	
/*** Distributes the nozzles with diameter on input length ***/
	distance = (IinRight - IinLeft-NNOZ*DNOZ)/(NNOZ-1);
	for ( n = 0; n <= NNOZ-1; n++ ) {
		for ( I = ( IinLeft + n*DNOZ + n*distance )*NPI + 1; I <=  (IinLeft + (n + 1)*DNOZ + n*distance )*NPI; I++) {
			Alpha[I][0] = ALPHAIN;
            v[I][1] = V_IN;
		}
	}

	//for (I = IinLeft*NPI + 1; I <= IinRight*NPI; I++) {
		//Alpha[I][0] = ALPHAIN; /* inlet */
	//}
	//} /* for I */

	for (I = 0; I <= NPI + 1; I++) {
		/* Temperature at the walls in Kelvin */
		T[I][0]     = 273.; /* bottom wall */
		T[I][NPJ+1] = 273.; /* top wall */
	} /* for I */

} /* fixedbound */

/* ################################################################# */
void bound(void)
/* ################################################################# */
{
/***** Purpose: Specify boundary conditions for a calculation ******/
	int    J, j, I, i;

//	if (solve_fi[NU] == TRUE) globcont();

	/* Velocity and temperature gradient at outlet = zero: */

	for (J = Joutmin; J <= Joutmax; J++) {
		j = J;
		/* Correction factor m_in/m_out is used to satisfy global continuity */
//		if (solve_fi[NU] == TRUE) u[NPI+1][J] = u[NPI][J]*m_in/m_out;
//		if (solve_fi[NV] == TRUE) v[NPI+1][j] = v[NPI][j];
	} /* for J */

	for (I = 0; I <= NPI+1; I++) {
		i = I;
		if (solve_fi[NT] == TRUE) T[I][NPJ+1] = T[I][NPJ];
		if (solve_fi[NALPHA] == TRUE) Alpha[I][NPJ+1] = Alpha[I][NPJ];
		if (solve_fi[NU] == TRUE) u[i][NPJ+1] = u[i][NPJ];
	} /* for I */

} /* bound */

/* ################################################################# */
void globcont(void)
/* ################################################################# */
{
/***** Purpose: Calculate mass in and out of the calculation domain to ******/
/*****          correct for the continuity at outlet. ******/
	int    J, j;
	double AREAw;

	conv();

	m_in = 0.;
	m_out = 0.;

	for (J = Jinmin; J <= Jinmax; J++) {
		j = J;
		AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
		m_in += F_u[1][J]*AREAw;
	} /* for J */

	for (J = Joutmin; J <= Joutmax; J++) {
		j = J;
		AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
		m_out += F_u[NPI][J]*AREAw;
	} /* for J */

} /* globcont */

/* ################################################################# */
void derivatives(void)
/* ################################################################# */
{
/***** Purpose: To calculate derivatives ******/
	int    I, J, i, j;
	
	for (I = Istart[NU] - 1; I <= Iend[NU]; I++) {
		i = I;
		for (J = Jstart[NU] - 1; J <= Jend[NU]; J++) {
			j = J;
			dudx[I][J] = (u[i+1][J] - u[i][J]) / (x_u[i+1] - x_u[i]);
		} /* for J */
	} /* for I */
	
	for (I = Istart[NU] - 1; I <= Iend[NU]; I++) {
		i = I;
		for (J = Jstart[NU]; J <= Jend[NU]; J++) {
			j = J;
			dudy[i][j] = (u[i][J]   - u[i][J-1]) / (y[J] - y[J-1]);
		} /* for J */
	} /* for I */
	
	for (I = Istart[NV]; I <= Iend[NV]; I++) {
		i = I;
		for (J = Jstart[NV] - 1; J <= Jend[NV]; J++) {
			j = J;
			dvdx[i][j] = (v[I][j]   - v[I-1][j]) / (x[I] - x[I-1]);
		} /* for J */
	} /* for I */

	for (I = Istart[NV] - 1; I <= Iend[NV]; I++) {
		i = I;
		for (J = Jstart[NV] - 1; J <= Jend[NV]; J++) {
			j = J;
			dvdy[I][J] = (v[I][j+1] - v[I][j]) / (y_v[j+1] - y_v[j]);
		} /* for J */
	} /* for I */
}

/* ################################################################# */
void solve(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP, int Istart, int Iend, int Jstart, int Jend, int ITER, double ratio)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. ******/
	int    I, J, space, iter;
	double *Ari, *Cmri, Cri;
	double residual, res0;
	
/* TDMA along a horizontal row from west to east. equation to solve: */

	/* - aW*fiW + aP*fiP - aE*fiE = aS*fiS + aN*fiN + b */

	/* equivalences with variables in eq. 7.1-7.6: */
	/* BETA = aW[I][J] Def. in eq. 7.2 */
	/* D    = aP[I][J] Def. in eq. 7.2 */
	/* ALFA = aE[I][J] Def. in eq. 7.2 */
	/* A    = Ari[I]	 Def. in eq. 7.6b */
	/* C    = Cri	 The right side assumed temporarily known (see eq. 7.8) */
	/* C´   = Cmri[I]  Def. in eq. 7.6c */
	/* b    = b[I][J]	 Def. in eq. 7.7 */

	space = max2((Iend - Istart + 3),(Jend - Jstart + 3));
	Ari   = double_1D_array(space);
	Cmri  = double_1D_array(space);

	residual = LARGE;
	res0 = SMALL;
	iter = 0;

	/* determine initial summed residual */

	for (J = Jstart; J <= Jend; J++)
		for (I = Istart; I <= Iend; I++)
			res0 += fabs((aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J]
	      		     - aP[I][J]*fi[I  ][J  ]));

	/* main inner iteration loop (aka Newton loop) */

	do {
		/* Solving the (e-w) lines from the south */

		for (J = Jstart; J <= Jend; J++) {
			/* At the inlet boundary: */
			Ari [Istart-1] = 0;
			Cmri[Istart-1] = fi[Istart-1][J];

			for (I = Istart; I <= Iend; I++) { /* Forward substitution */
				Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
				Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
				Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6c */
			}

			for (I = Iend; I >= Istart; I--)  /* Back substitution */
				fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
		} /* for J from south */

		/* Solving the (e-w) lines from the north */

		for (J = Jend-1; J >= Jstart; J--) {
			/* At the inlet boundary: */
			Ari [Istart-1] = 0;
			Cmri[Istart-1] = fi[Istart-1][J];

			for (I = Istart; I <= Iend; I++) { /* Forward substitution */
				Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
				Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
				Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]);  /* eq. 7.6c */
			}

			for (I = Iend; I >= Istart; I--)  /* Back substitution */
				fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
		} /* for J from north */

	/* TDMA along a vertical column from south to north. equation to solve: */

		/* - aS*fiW + aP*fiP - aN*fiE = aW*fiS + aE*fiN + b (eq. 7.8) */

		/* equivalences with variables in eq. 7.1-7.6: */
		/* BETA = aS[I][J] Def. in eq. 7.2 */
		/* D    = aP[I][J] Def. in eq. 7.2 */
		/* ALFA = aN[I][J] Def. in eq. 7.2 */
		/* A    = Ari[I]	 Def. in eq. 7.6b */
		/* C    = Cri      The right side assumed temporarily known (see eq. 7.8) */
		/* C´   = Cmri[I]  Def. in eq. 7.6c */
		/* b    = b[I][J]	 Def. in eq. 7.7 */

		/* Solving (n-s) lines from the west */

		for (I = Istart; I <= Iend; I++) {
			/* At the bottom boundary: */
			Ari[Jstart-1] = 0;
			Cmri[Jstart-1] = fi[I][Jstart-1];

			for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
				Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
				Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
				Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
			}

			for (J = Jend; J >= Jstart; J--) /* Back substitution */
				fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
		} /* for I from west */

		/* Solving (n-s) lines from the east */

		for (I = Iend - 1; I >= Istart; I--) {
			/* At the bottom boundary: */
			Ari[Jstart-1] = 0;
			Cmri[Jstart-1] = fi[I][Jstart-1];

			for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
				Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
				Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
				Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
			}

			for (J = Jend; J >= Jstart; J--) /* Back substitution */
				fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
		} /* for I from east */

		/* reset residual */
		residual =0.;
		
		for (J = Jstart; J <= Jend; J++)
			for (I = Istart; I <= Iend; I++)
				residual += fabs((aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
				           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J]
				           - aP[I][J]*fi[I  ][J  ]));

		iter++;
	} while(iter < ITER && residual/res0 > ratio);

	/* free memory */
	
	free(Ari);
	free(Cmri);

//	printf("iter = %2d\n",iter);
}

/* ################################################################# */
void solveGS(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP, int Istart, int Iend, int Jstart, int Jend, int ITER, double ratio)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel ******/
	int    I, J, iter;
	double residual, res0;

	residual = LARGE;
	res0 = SMALL;
	iter = 0;

	/* determine initial summed residual */

	for (J = Jstart; J <= Jend; J++)
		for (I = Istart; I <= Iend; I++)
			res0 += fabs((aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J]
	      		     - aP[I][J]*fi[I  ][J  ]));

	/* main inner iteration loop (aka Newton loop) */

	do {
		for (I = Istart; I <= Iend; I++)
			for (J = Jstart; J <= Jend; J++)
				fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
				           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			      	     /aP[I][J];

		/* reset residual */
		residual =0.;
		
		for (J = Jstart; J <= Jend; J++)
			for (I = Istart; I <= Iend; I++)
				residual += fabs((aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
				           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J]
				           - aP[I][J]*fi[I  ][J  ]));

		iter++;
	} while(iter < ITER && residual/res0 > ratio);

} /* solveGS */

/* ################################################################# */
void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP, int Istart, int Iend, int Jstart, int Jend)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel ******/
	int    I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           /aP[I][J] * omega - (omega - 1)*fi[I][J];
} /* solveSOR */

/* ################################################################# */
void solveSIP_new(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP, int Istart, int Iend, int Jstart, int Jend)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Stone's Strong Implicit Solver (SIP) ******/
	int	I, J, it, MAXIT;
      double **LW,**LS,**LPR;
      double **UE,**UN,**RES;
//      double LW[NPI+1][NPJ+1],LS[NPI+1][NPJ+1],LPR[NPI+1][NPJ+1];
//      double UE[NPI+1][NPJ+1],UN[NPI+1][NPJ+1],RES[NPI+1][NPJ+1];
      double P1, P2, alpha, RES0, RESit, RSM, RESMAX;

	it = 0;
	RES0 = 0.;
	RESMAX = 1E-6;
	MAXIT = 2;
	alpha = 0.92;
	RSM = 1.;

	LW  = double_2D_matrix((Iend + 2),(Jend + 2));
	LS  = double_2D_matrix((Iend + 2),(Jend + 2));
	LPR = double_2D_matrix((Iend + 2),(Jend + 2));
	UE  = double_2D_matrix((Iend + 2),(Jend + 2));
	UN  = double_2D_matrix((Iend + 2),(Jend + 2));
	RES = double_2D_matrix((Iend + 2),(Jend + 2));

	for (I = Istart; I <= Iend; I++) {
		for (J = Jstart; J <= Jend; J++) {
			LW[I][J]  = 0.;
			LS[I][J]  = 0.;
			UE[I][J]  = 0.;
			UN[I][J]  = 0.;
			RES[I][J] = 0.;
		}
	}

/* CALCULATE ELEMENTS OF [L] AND [U] MATRICES */

	for (I = Istart; I <= Iend; I++) {
		for (J = Jstart; J <= Jend; J++) {
			LW[I][J] = -aW[I][J]/(1 + alpha*UN[I-1][J  ]);
			LS[I][J] = -aS[I][J]/(1 + alpha*UE[I  ][J-1]);
			P1 = alpha*LW[I][J]*UN[I-1][J  ];
			P2 = alpha*LS[I][J]*UE[I  ][J-1];
			LPR[I][J] = 1./(aP[I][J] + P1 + P2 - LW[I][J]*UE[I-1][J] - LS[I][J]*UN[I][J-1]);
			UN[I][J] = (-aN[I][J] - P1)*LPR[I][J];
			UE[I][J] = (-aE[I][J] - P2)*LPR[I][J];
		}
	}

/* CALCULATE RESIDUAL AND AUXILLIARY VECTORS; INNER ITERATION LOOP */

	while(RSM > RESMAX && it < MAXIT) {
		RESit = 0.;
		for (I = Istart; I <= Iend; I++) {
			for (J = Jstart; J <= Jend; J++) {
				RES[I][J] = b[I][J] - aP[I][J]*fi[I][J] + aN[I][J]*fi[I][J+1]
				            + aS[I][J]*fi[I][J-1] + aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J];

				RESit += fabs(RES[I][J]);
				RES[I][J] = (RES[I][J] - LS[I][J]*RES[I][J-1] - LW[I][J]*RES[I-1][J])*LPR[I][J];
			}
		}

		if(it == 0) RES0 = RESit;

/* CALCULATE INCREMENT AND CORRECT VARIABLE */

		for (I = Iend - 1; I >= Istart; I--) { /* Back substitution */
			for (J = Jend; J >= Jstart; J--) { /* Back substitution */
				RES[I][J] += -UN[I][J]*RES[I][J+1] - UE[I][J]*RES[I+1][J];
				fi [I][J] += RES[I][J];
			}
		}

/* CONVERGENCE CHECK */

		RSM = RESit/(RES0 + SMALL);
		it++;
	} /* while it */
/*	printf("Sweep %d, RSM = %e\n", it, RSM);*/
} /* solveSIP_new */

/* ################################################################# */
void conv(void)
/* ################################################################# */
{
/***** Purpose: To calculate the convective mass flux component pr. unit ******/
/*****          area defined in eq. 5.7 ******/
	int    I, J, i, j;

	for (I = 1; I <= NPI + 1; I++) {
		i = I;
		for (J = 1; J <= NPJ + 1; J++) {
			j = J;
//			F_u[i][J] = (rho[I-1][J  ]*(x[I] - x_u[i]) + rho[I][J]*(x_u[i] - x[I-1]))*u[i][J]/(x[I] - x[I-1]); /* = F(i, J) */
//			F_v[I][j] = (rho[I  ][J-1]*(y[J] - y_v[j]) + rho[I][J]*(y_v[j] - y[J-1]))*v[I][j]/(y[J] - y[J-1]); /* = F(I, j) */
			F_u[i][J] = RHOL*u[i][J] + (RHOG*Alpha[I-1][J  ]*(x[I] - x_u[i]) + RHOG*Alpha[I][J]*(x_u[i] - x[I-1]))*USLIP/(x[I] - x[I-1]); /* = F(i, J) */
			F_v[I][j] = RHOL*v[I][j]; /* = F(I, j) */
		} /* for J */
	} /* for I */

} /* conv */

/* ################################################################# */
void ucoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b, int Istart, int Iend, int Jstart , int Jend,int *fixed)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the u equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold, mun, mus;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
			AREAe = AREAw;
			AREAs = x[I] - x[I-1];
			AREAn = AREAs;

			/* eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = ((F_u[i  ][J  ] + F_u[i-1][J  ])/2)*AREAw;
			Fe = ((F_u[i+1][J  ] + F_u[i  ][J  ])/2)*AREAe;
			Fs = ((F_v[I  ][j  ] + F_v[I-1][j  ])/2)*AREAs;
			Fn = ((F_v[I  ][j+1] + F_v[I-1][j+1])/2)*AREAn;

			/* eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b  */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = (mu[I-1][J]/(x_u[i  ] - x_u[i-1]))*AREAw;
			De = (mu[I  ][J]/(x_u[i+1] - x_u[i  ]))*AREAe;
			Ds = ((mu[I-1][J  ] + mu[I][J  ] + mu[I-1][J-1] + mu[I][J-1])/(4*(y[J  ] - y[J-1])))*AREAs;
			Dn = ((mu[I-1][J+1] + mu[I][J+1] + mu[I-1][J  ] + mu[I][J  ])/(4*(y[J+1] - y[J  ])))*AREAn;

			/* The source terms */

			mus = 0.25*(mu[I][J] + mu[I-1][J] + mu[I][J-1] + mu[I-1][J-1]);
			mun = 0.25*(mu[I][J] + mu[I-1][J] + mu[I][J+1] + mu[I-1][J+1]);
			SP[i][J] = 0.;
			Su[i][J] = (mu[I][J]*dudx[I][J] - mu[I-1][J]*dudx[I-1][J])/(x[I] - x[I-1]) + 
			           (mun*dvdx[i][j+1] - mus*dvdx[i][j])/(y_v[j+1] - y_v[j]);
//			Su[i][J] += -9.81*0.5*(rho[I-1][J] + rho[I][J]);
			Su[i][J] *= AREAw*AREAs;

			if (j > *(fixed) && j <= *(fixed+1) && ( (I > *(fixed+2) && I < *(fixed+3))|| (I > *(fixed+4) && I <= *(fixed+5))  )){
				SP[i][J] = -LARGE;
            }

			/* The coefficients (hybrid differencing sheme) */

			aW[i][J] = max3( Fw, Dw + Fw/2, 0.);
			aE[i][J] = max3(-Fe, De - Fe/2, 0.);
			aS[i][J] = max3( Fs, Ds + Fs/2, 0.);
			aN[i][J] = max3(-Fn, Dn - Fn/2, 0.);


			aPold    = 0.5*(rho_old[I-1][J] + rho_old[I][J])*AREAe*AREAn/Dt;

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[i][J] = aW[i][J] + aE[i][J] + aS[i][J] + aN[i][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Calculation of d[i][J] = d_u[i][J] defined in eq. 6.23 for use in the  */
			/* equation for pression correction (eq. 6.32). See subroutine pccoeff. */

			d_u[i][J] = AREAw*relax[NU]/aP[i][J];

			/* Putting the integrated pressure gradient into the source term b[i][J] */
			/* The reason is to get an equation on the generalised form  */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm.  */
			/* note: In reality b = a0p*fiP + Su = 0.  */

			b[i][J] = (p[I-1][J] - p[I][J])*AREAw + Su[I][J] + aPold*u_old[i][J];

			/* Introducing relaxation by eq. 6.36 . and putting also the last  */
			/* term on the right side into the source term b[i][J] */

			aP[i][J] /= relax[NU];
			b [i][J] += (1 - relax[NU])*aP[i][J]*u[i][J];

			/* now we have implemented eq. 6.36 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done  */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */

} /* ucoeff */

/* ################################################################# */
void vcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b, int Istart, int Iend, int Jstart , int Jend, int *fixed)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the v equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold, mue, muw;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y[J] - y[J-1]; /* See fig. 6.4 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i];
			AREAn = AREAs;

			/* eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = ((F_u[i  ][J] + F_u[i  ][J-1])/2)*AREAw;
			Fe = ((F_u[i+1][J] + F_u[i+1][J-1])/2)*AREAe;
			Fs = ((F_v[I  ][j] + F_v[I  ][j-1])/2)*AREAs;
			Fn = ((F_v[I  ][j] + F_v[I  ][j+1])/2)*AREAn;

			/* eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = ((mu[I-1][J-1] + mu[I  ][J-1] + mu[I-1][J] + mu[I  ][J])/(4*(x[I  ] - x[I-1])))*AREAw;
			De = ((mu[I  ][J-1] + mu[I+1][J-1] + mu[I  ][J] + mu[I+1][J])/(4*(x[I+1] - x[I  ])))*AREAe;
			Ds =  (mu[I][J-1]/(y_v[j  ] - y_v[j-1]))*AREAs;
			Dn =  (mu[I][J  ]/(y_v[j+1] - y_v[j  ]))*AREAn;

			/* The source terms */

			muw = 0.25*(mu[I][J] + mu[I-1][J] + mu[I][J-1] + mu[I-1][J-1]);
			mue = 0.25*(mu[I][J] + mu[I+1][J] + mu[I][J-1] + mu[I+1][J-1]);
			SP[I][j] = 0.;
			Su[I][j] = 0.;
			Su[I][j] = (mu[I][J]*dvdy[I][J] - mu[I][J-1]*dvdy[I][J-1])/(y[J] - y[J-1]) + 
			           (mue*dudy[i+1][j] - muw*dudy[i][j])/(x_u[i+1] - x_u[i]);
			Su[I][j] += -9.81*0.5*(rho[I][J-1] + rho[I][J]);
			Su[I][j] *= AREAw*AREAs;

			/* v can be fixed to zero by setting SP to a very large value */

            //https://www.horecasupply.nl/producten/glaswerk/bierglazen?prodId=20857

			if (j > *(fixed) && j <= *(fixed+1) && ( (I > *(fixed+2) && I < *(fixed+3))|| (I > *(fixed+4) && I <= *(fixed+5))  )){
				SP[I][j] = -LARGE;
            }

//			if (J == 2*NPJ/5 && i > 2*NPJ/3)
		//		SP[I][j] = -LARGE;

			/* The coefficients (hybrid differencing sheme) */
      		aN[I][j] = max3(-Fn, Dn - Fn/2, 0.);
            aS[I][j] = max3( Fs, Ds + Fs/2, 0.);
			aW[I][j] = max3( Fw, Dw + Fw/2, 0.);
			aE[I][j] = max3(-Fe, De - Fe/2, 0.);
			
			aPold    = 0.5*(rho_old[I][J-1] + rho_old[I][J])*AREAe*AREAn/Dt;

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[I][j] = aW[I][j] + aE[I][j] + aS[I][j] + aN[I][j] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Calculation of d[I][j] = d_v[I][j] defined in eq. 6.23 for use in the */
			/* equation for pression correction (eq. 6.32) (see subroutine pccoeff). */

			d_v[I][j] = AREAs*relax[NV]/aP[I][j];

			/* Putting the integrated pressure gradient into the source term b[I][j] */
			/* The reason is to get an equation on the generalised form */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm. */
			/* note: In reality b = a0p*fiP + Su = 0. */

			b[I][j] = (p[I][J-1] - p[I][J])*AREAs + Su[I][j] + aPold*v_old[I][j];

			/* Introducing relaxation by eq. 6.37 . and putting also the last */
			/* term on the right side into the source term b[i][J] */

			aP[I][j] /= relax[NV];
			b [I][j] += (1 - relax[NV])*aP[I][j]*v[I][j];

			/* now we have implemented eq. 6.37 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */

} /* vcoeff */

/* ################################################################# */
void pccoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b, int Istart, int Iend, int Jstart , int Jend)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the pc equation. ******/
	int    i, j, I, J;
	double AREAw, AREAe, AREAs, AREAn;
	double SSUM;

	SMAX = 0.;
	SSUM = 0.;
	SAVG = 0.;
	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The constant b´ in eq. 6.32 */

			b[I][J] = F_u[i][J]*AREAw - F_u[i+1][J]*AREAe + F_v[I][j]*AREAs - F_v[I][j+1]*AREAn
			          + (rho_old[I][J] - rho[I][J])*AREAe*AREAs/Dt;

			SP[I][J] = 0.;
			Su[I][J] = 0.;
			
			b[I][J] += Su[I][J];

		      SMAX     = max2(SMAX,fabs(b[I][J]));
		      SSUM    += fabs(b[I][J]);
			
			/* The coefficients */

			aE[I][J] = (rho[I  ][J  ] + rho[I+1][J  ])*d_u[i+1][J  ]*AREAe/2;
			aW[I][J] = (rho[I-1][J  ] + rho[I  ][J  ])*d_u[i  ][J  ]*AREAw/2;
			aN[I][J] = (rho[I  ][J  ] + rho[I  ][J+1])*d_v[I  ][j+1]*AREAn/2;
			aS[I][J] = (rho[I  ][J-1] + rho[I  ][J  ])*d_v[I  ][j  ]*AREAs/2;

			aP[I][J] = aE[I][J] + aW[I][J] + aN[I][J] + aS[I][J] - SP[I][J];

			pc[I][J] = 0.;

			/* note: At the points nearest the boundaries, some coefficients are */
			/* necessarily zero. For instance at I = 1 and J = 1, the coefficients */
			/* aS and aW are zero since they are on the outside of the calculation */
			/* domain. This is automatically satisfied through the initialisation */
			/* where d_u[i][J] and d_v[I][j] are set to zero at these points. */

			} /* for J */
		} /* for I */

		/* Average error in mass balance is summed error devided by */
		/* number of internal grid points */
	      SAVG = SSUM/((Iend - Istart)*(Jend - Jstart));

} /* pccoeff */



/* ################################################################# */
void storeresults(void)
/* ################################################################# */
{
/***** To newly calculated variables are stored in the arrays ******/
/***** for old variables, which can be used in the next timestep ******/
	int I, J, i, j;

	/* store u velocity */

	if (solve_fi[NU] == TRUE)
		for(i = Istart[NU]; i <= Iend[NU]; i++)
			for(J = Jstart[NU]; J <= Jend[NU]; J++)
				u_old[i][J] = u[i][J];

	/* store v velocity */

	if (solve_fi[NV] == TRUE)
		for(I = Istart[NV]; I <= Iend[NV]; I++)
			for(j = Jstart[NV]; j <= Jend[NV]; j++)
				v_old[I][j] = v[I][j];

	/* store pressure correction */

	if (solve_fi[NPC] == TRUE)
		for(I = Istart[NPC]; I <= Iend[NPC]; I++)
			for(J = Jstart[NPC]; J <= Jend[NPC]; J++)
				pc_old[I][J] = pc[I][J];

	/* store temperature */

	if (solve_fi[NT] == TRUE)
		for(I = Istart[NT]; I <= Iend[NT]; I++)
			for(J = Jstart[NT]; J <= Jend[NT]; J++)
				T_old[I][J] = T[I][J];

	/* store gas volume fraction */

	if (solve_fi[NALPHA] == TRUE)
		for(I = Istart[NALPHA]; I <= Iend[NALPHA]; I++)
			for(J = Jstart[NALPHA]; J <= Jend[NALPHA]; J++)
				Alpha_old[I][J] = Alpha[I][J];

} /* storeresults */


/* ################################################################# */
void calcmean(void)
/* ################################################################# */
{
/***** Purpose: Creating time averaged result table ******/
	int    I, J, i, j;
	double ugrid, vgrid;
	FILE   *fp;
	char   meandata[200];
	int    gnuplot_used = 1;
	int    Nzero = (int) 10./Dt - 1;

	if (NDt > Nzero) {

		for (I = 0; I <= NPI; I++) {
			i = I;
			for (J = 1; J <= NPJ; J++) {
				j = J;
				ugrid = 0.5*(u[i][J]+u[i+1][J  ]);
				vgrid = 0.5*(v[I][j]+v[I  ][j+1]);
				usum[I][J]     += ugrid;
				vsum[I][J]     += vgrid;
				psum[I][J]     += p[I][J];
				rhosum[I][J]   += rho[I][J];
				musum[I][J]    += mu[I][J];
				alphasum[I][J] += Alpha[I][J];

				umean[I][J]     = usum[I][J]/(NDt-Nzero);
				vmean[I][J]     = vsum[I][J]/(NDt-Nzero);
				pmean[I][J]     = psum[I][J]/(NDt-Nzero);
				rhomean[I][J]   = rhosum[I][J]/(NDt-Nzero);
				mumean[I][J]    = musum[I][J]/(NDt-Nzero);
				alphamean[I][J] = alphasum[I][J]/(NDt-Nzero);
			} /* for J */
		} /* for I */

		/* meandata = path outputdir + meandata.dat */
		sprintf (meandata, "%smeandata.dat", outputdir);

		if(NDt%DUMPFREQ == 0) {
//			if ((fp = fopen(meandata,"w"))==NULL) {
//				printf("could not open file: %s \n", meandata);
//				exit(1);
//			}
			while ((fp = fopen(meandata,"w"))==NULL);

			fprintf(fp, "#Time    %11.4e\n", current_time);
			fprintf(fp, "#x\t\ty\t\tugrid\t\tvgrid\t\tp\t\trho\t\tmu\t\tAlpha\n");

			for (I = 1; I <= NPI; I++) {
				i = I;
				for (J = 1; J <= NPJ; J++) {
					j = J;
					fprintf(fp, "%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
				      	       x[I], y[J], umean[I][J], vmean[I][J], pmean[I][J], rhomean[I][J], mumean[I][J], alphamean[I][J]);
//							 1     2     3      4      5        6          7         8
				} /* for J */
				if (gnuplot_used == 1) fprintf(fp, "\n");
			} /* for I */

			fclose(fp);
		} /* NDt%DUMPFREQ == 0 */
	} /* if NDt > Nzero */

} /* calcmean */

/* ################################################################# */
void velcorr(void)
/* ################################################################# */
{
/***** To correct the pressure and the velocities by eq. 6.24, 6.25 ******/
/*****  and a modified version of eq. 6.33. ******/
	int I, J, i, j;

	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;

			p[I][J] += relax[NPC]*pc[I][J]; /* equation 6.33 */

			/* Velocity correction */
			/* Note: the relaxation factors for u and v are included  */
			/* in the d_u and d_v terms (see page 146) */

			if (i != 1)
				u[i][J] += d_u[i][J]*(pc[I-1][J  ] - pc[I][J]); /* eq. 6.24 */

			if (j != 1)
				v[I][j] += d_v[I][j]*(pc[I  ][J-1] - pc[I][J]); /* eq. 6.25 */

		} /* for J */
	} /* for I */

} /* velcorr */

/* ################################################################# */
void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b, int Istart, int Iend, int Jstart , int Jend)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       Dw, De, Ds, Dn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*AREAw;
			Fe = F_u[i+1][J  ]*AREAe;
			Fs = F_v[I  ][j  ]*AREAs;
			Fn = F_v[I  ][j+1]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			/* The conductivity, Gamma, at the interface is calculated */
			/* with the use of a harmonic mean. */

			Dw = ((Gamma[I-1][J  ]*Gamma[I  ][J  ])/(Gamma[I-1][J  ]*(x[I  ] - x_u[i  ]) + Gamma[I  ][J  ]*(x_u[i  ] - x[I-1])))*AREAw;
			De = ((Gamma[I  ][J  ]*Gamma[I+1][J  ])/(Gamma[I  ][J  ]*(x[I+1] - x_u[i+1]) + Gamma[I+1][J  ]*(x_u[i+1] - x[I  ])))*AREAe;
			Ds = ((Gamma[I  ][J-1]*Gamma[I  ][J  ])/(Gamma[I  ][J-1]*(y[J  ] - y_v[j  ]) + Gamma[I  ][J  ]*(y_v[j  ] - y[J-1])))*AREAs;
			Dn = ((Gamma[I  ][J  ]*Gamma[I  ][J+1])/(Gamma[I  ][J  ]*(y[J+1] - y_v[j+1]) + Gamma[I  ][J+1]*(y_v[j+1] - y[J  ])))*AREAn;

			/* The source terms */

			SP[I][J] = 0.;
			Su[I][J] = 0.;

			/* The coefficients (hybrid differencing sheme) */

			aW[I][J] = max3( Fw, Dw + Fw/2, 0.);
			aE[I][J] = max3(-Fe, De - Fe/2, 0.);
			aS[I][J] = max3( Fs, Ds + Fs/2, 0.);
			aN[I][J] = max3(-Fn, Dn - Fn/2, 0.);
			aPold    = rho_old[I][J]*AREAe*AREAn/Dt;
//			if (J == Jstart) aS[I][J] = 0.;
//			if (J == Jend)   aN[I][J] = 0.;
//			if (I == Istart) aW[I][J] = 0.;

			if (i > 11 && i < 18 && J > 2*NPJ/5 && J < 3*NPJ/5){
				SP[i][J] = -LARGE;
				Su[i][J] = LARGE*373.;
			}

			/* eq. 8.31 with time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J] + aPold*T_old[I][J];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */

} /* Tcoeff */

/* ################################################################# */
void Alphacoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b, int Istart, int Iend, int Jstart , int Jend, int *fixed)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn, 
	       AREAw, AREAe, AREAs, AREAn,
	       aPold,
	       alphae, alphaw, alphan, alphas, alpha1, alpha2, alpha3;

	alphae = 0.;
	alphaw = 0.;
	alphan = 0.;
	alphas = 0.;
	
	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = u[i  ][J  ]*AREAw;
			Fe = u[i+1][J  ]*AREAe;
//			Fs = (v[I  ][j  ] + USLIP)*AREAs;
//			Fn = (v[I  ][j+1] + USLIP)*AREAn;
			Fs = (v[I  ][j  ] + USLIP*RHOL*0.5*(2. - Alpha[I][J] - Alpha[I][J-1])/(rho[I][J-1] + rho[I][J]))*AREAs;
			Fn = (v[I  ][j+1] + USLIP*RHOL*0.5*(2. - Alpha[I][J] - Alpha[I][J+1])/(rho[I][J+1] + rho[I][J]))*AREAn;

//			if (I == 1) Fw = (0.01 - (p[I] - p[I-1])/(x[I  ] - x[I-1])/Cdrag)*AREAw;

			/* The source terms */


            if (j > *(fixed) && j <= *(fixed+1) && ( (I > *(fixed+2) && I < *(fixed+3))|| (I > *(fixed+4) && I <= *(fixed+5))  )){
                SP[I][J] = -LARGE;
            }
            else{
                SP[I][J] = 0.;
            }

			Su[I][J] = 0.;

			/* The coefficients (upwind differencing sheme) */

            if (j > *(fixed) && j <= *(fixed+1) && ( (I > *(fixed+2) && I < *(fixed+3))|| (I > *(fixed+4) && I <= *(fixed+5))  )){
                aN[I][J]=0.;
                aS[I][J]=0.;
                aW[I][J]=0.;
                aE[I][J]=0.;
            }
            else{
			    aW[I][J] = max2( Fw, 0.);
			    aE[I][J] = max2(-Fe, 0.);
			    aS[I][J] = max2( Fs, 0.);
			    aN[I][J] = max2(-Fn, 0.);
            }

			aPold    = AREAe*AREAn/Dt;

			/* deferred correction alphaw (Barton minus Upwind) */
		if (I > Istart) {
			if (u[i][J] >= 0) {
				alpha1 = 0.5*(Alpha[I-1][J] - Alpha[I-2][J]);
				alpha2 = 0.5*(Alpha[I][J] - Alpha[I-1][J]);
				alpha3 = 0.;
				if (Alpha[I][J] <= Alpha[I-1][J])
					alphaw = min2(alpha3, max2(alpha1, alpha2));
				else	alphaw = max2(alpha3, min2(alpha1, alpha2));
			} else /* u[i][J] < 0 */ {
				alpha1 = 0.5*(Alpha[I][J] - Alpha[I+1][J]);
				alpha2 = 0.5*(Alpha[I-1][J] - Alpha[I][J]);
				alpha3 = 0.;
				if (Alpha[I][J] <= Alpha[I-1][J])
					alphaw = max2(alpha3, min2(alpha1, alpha2));
				else	alphaw = min2(alpha3, max2(alpha1, alpha2));
			} /* if */
		}

			/* deferred correction alphae (Barton minus Upwind) */
		if (I < Iend) {
			if (u[i+1][J] >= 0) {
				alpha1 = 0.5*(Alpha[I][J] - Alpha[I-1][J]);
				alpha2 = 0.5*(Alpha[I+1][J] - Alpha[I][J]);
				alpha3 = 0.;
				if (Alpha[I+1][J] <= Alpha[I][J])
					alphae = min2(alpha3, max2(alpha1, alpha2));
				else	alphae = max2(alpha3, min2(alpha1, alpha2));
			} else /* u[i+1][J] < 0 */ {
				alpha1 = 0.5*(Alpha[I+1][J] - Alpha[I+2][J]);
				alpha2 = 0.5*(Alpha[I][J] - Alpha[I+1][J]);
				alpha3 = 0.;
				if (Alpha[I+1][J] <= Alpha[I][J])
					alphae = max2(alpha3, min2(alpha1, alpha2));
				else	alphae = min2(alpha3, max2(alpha1, alpha2));
			} /* if */
		}
		
			/* deferred correction alphas (Barton minus Upwind) */
		if (J > Jstart) {
			if (v[I][j] >= 0) {
				alpha1 = 0.5*(Alpha[I][J-1] - Alpha[I][J-2]);
				alpha2 = 0.5*(Alpha[I][J] - Alpha[I][J-1]);
				alpha3 = 0.;
				if (Alpha[I][J] <= Alpha[I][J-1])
					alphas = min2(alpha3, max2(alpha1, alpha2));
				else	alphas = max2(alpha3, min2(alpha1, alpha2));
			} else /* v[I][j] < 0 */ {
				alpha1 = 0.5*(Alpha[I][J] - Alpha[I][J+1]);
				alpha2 = 0.5*(Alpha[I][J-1] - Alpha[I][J]);
				alpha3 = 0.;
				if (Alpha[I][J] <= Alpha[I][J-1])
					alphas = max2(alpha3, min2(alpha1, alpha2));
				else	alphas = min2(alpha3, max2(alpha1, alpha2));
			} /* if */
		}
		
			/* deferred correction alphan (Barton minus Upwind) */
		if (J < Jend) {
			if (v[I][j+1] >= 0) {
				alpha1 = 0.5*(Alpha[I][J] - Alpha[I][J-1]);
				alpha2 = 0.5*(Alpha[I][J+1] - Alpha[I][J]);
				alpha3 = 0.;
				if (Alpha[I][J+1] <= Alpha[I][J])
					alphan = min2(alpha3, max2(alpha1, alpha2));
				else	alphan = max2(alpha3, min2(alpha1, alpha2));
			} else /* v[I][j+1] < 0 */ {
				alpha1 = 0.5*(Alpha[I][J+1] - Alpha[I][J+2]);
				alpha2 = 0.5*(Alpha[I][J] - Alpha[I][J+1]);
				alpha3 = 0.;
				if (Alpha[I][J+1] <= Alpha[I][J])
					alphan = max2(alpha3, min2(alpha1, alpha2));
				else	alphan = min2(alpha3, max2(alpha1, alpha2));
			} /* if */
		}

			/* eq. 8.31 with time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J] + aPold;

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J] + aPold*Alpha_old[I][J];
			b[I][J] -= Fe*alphae - Fw*alphaw + Fn*alphan - Fs*alphas;

			aP[I][j] /= relax[NALPHA];
			b [I][j] += (1 - relax[NALPHA])*aP[I][j]*Alpha[I][j];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */

} /* Alphacoeff */

/* ################################################################# */
void density(void)
/* ################################################################# */
{
/***** Purpose: Calculate the density rho(I, J) in the fluid as a function *****/
/*****          of the ideal gas law. Note: rho at the walls are not needed *****/
/*****          in this case, and therefore not calculated. *****/
	int    I, J;

	for (I = 0; I <= NPI + 1; I++) {
		for (J = 0; J <= NPJ + 1; J++) {
			rho_old[I][J] = rho[I][J];
			rho[I][J] = (1. - Alpha[I][J])*RHOL + Alpha[I][J]*RHOG;
		} /* for J */
	} /* for I */

} /* density */

/* ################################################################# */
void viscosity(void)
/* ################################################################# */
{
/***** Purpose: Calculate the viscosity in the fluid as a function of *****/
/*****          temperature. Max error in the actual temp. interval is 0.5% *****/
	int    I, J, i, j;
	double crossterms, EijEij;

	for (I = 0; I <= NPI; I++)
		for (J = 1; J <= NPJ + 1; J++) {
//			mu[I][J] = (0.0395*T[I][J] + 6.58)*1.E-6;
			i = I;
			j = J;
			crossterms = 0.25*(dudy[i][j] + dudy[i+1][j] + dudy[i][j+1] + dudy[i+1][j+1] +
      			             dvdx[i][j] + dvdx[i+1][j] + dvdx[i][j+1] + dvdx[i+1][j+1]);
			EijEij  = 2.0*(sqr(dudx[I][J]) + sqr(dvdy[I][J])) + sqr(crossterms);
			mu[I][J] = 1.e-3 + sqr(CZERO*DELTA)*sqrt(EijEij)*RHOL;
			if (BIT == TRUE) mu[I][J] += CMUBIT * Alpha[I][J] * DBUB * USLIP;
	} /* for J */

} /* viscosity */

/* ################################################################# */
void conductivity(void)
/* ################################################################# */
{
/***** Purpose: Calculate the thermal conductivity in the fluid as a  ******/
/*****          function of temperature. Max error in the actual ******/
/*****          temperature interval is 0.9% ******/

	int    I, J;

	for (I = 0; I <= NPI; I++)
		for (J = 1; J <= NPJ; J++) {
			Gamma[I][J] = (6.1E-5*T[I][J] + 8.4E-3)/Cp[I][J];
			if (Gamma[I][J] < 0.){
				output();
				fprintf(stderr, "Error: Gamma[%d][%d] = %e\n", I, J, Gamma[I][J]);
				exit(1);
			} /* if */
		} /* for J */

} /* conductivity */


/* ################################################################# */
void output(void)
/* ################################################################# */
{
/***** Purpose: Creating result table ******/
	int    I, J, i, j;
	double ugrid, vgrid,stream;
	FILE   *fp, *str, *velu, *velv;

	fp = fopen("output.dat", "w");
	str = fopen("str.dat", "w");
	velu = fopen("velu.dat", "w");
	velv = fopen("velv.dat", "w");

	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			ugrid = 0.5*(u[i][J]+u[i+1][J  ]);
			vgrid = 0.5*(v[I][j]+v[I  ][j+1]);
			fprintf(fp, "%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
			             x[I], y[J], ugrid*XMAX/U_IN/NPI, vgrid*YMAX/U_IN/NPJ, p[I][J], T[I][J], rho[I][J], mu[I][J], Gamma[I][J], Alpha[I][J], u[i][J] + (p[I] - p[I-1])/(x[I] - x[I-1])/Cdrag);
//					 1     2     3                    4                    5        6        7          8         9            10           11
		} /* for J */
		fprintf(fp, "\n");
	} /* for I */

	fclose(fp);
	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			stream = -(u[i][J+1]-u[i][J])/(y[J+1]-y[J])
			         +(v[I+1][j]-v[I][j])/(x[I+1]-x[I]);
			fprintf(str, "%10.2e\t%10.2e\t%10.5e\n",
			             x_u[i], y_v[j], stream);
			fprintf(velu, "%10.2e\t%10.2e\t%10.5e\n",
			             x_u[i], y[J], u[i][J]);
			fprintf(velv, "%10.2e\t%10.2e\t%10.5e\n",
			             x[I], y_v[j], v[I][j]);
		} /* for J */
		fprintf(str, "\n");
		fprintf(velu, "\n");
		fprintf(velv, "\n");
	} /* for I */

	fclose(str);
	fclose(velu);
	fclose(velv);

} /* print */

/* ################################################################# */
void outputmanager(int number)
/* ################################################################# */
{

	int    I, J, i, j;
	double ugrid, vgrid;
	FILE *fp;
	char dumpname[200];
	int gnuplot_used = 1;

	/* dumpname = path outputdir + d12345.dat */
	sprintf (dumpname, "%sd%05i.dat", outputdir, number);

//	if ((fp = fopen(dumpname,"w"))==NULL) {
//		printf("could not open file: %s \n", dumpname);
//		exit(1);
//	}
	while ((fp = fopen(dumpname,"w"))==NULL);

	fprintf(fp, "#Time    %11.4e\n", current_time);
	fprintf(fp, "#x\t\ty\t\tugrid\t\tvgrid\t\tp\t\trho\t\tmu\t\tAlpha\n");
	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			ugrid = 0.5*(u[i][J]+u[i+1][J  ]);
//			vgrid = 0.5*(v[I][j]+v[I  ][j+1]);
			vgrid = 0.5*(v[I][j]+v[I  ][j+1]) - Alpha[I][J]*RHOG*USLIP/rho[I][J];
			fprintf(fp, "%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",
			             x[I], y[J], ugrid, vgrid, p[I][J], rho[I][J], mu[I][J], Alpha[I][J]);
//					 1     2     3      4      5        6          7         8
		} /* for J */
		if (gnuplot_used == 1) fprintf(fp, "\n");
	} /* for I */

	fclose(fp);

} /* outputmanager */

/* ################################################################# */
void readconstants(void)
/* ################################################################# */
{
/***** Purpose: reads the constants that are not defined as constants ******/
/***** in the constants.h file from the asm.inp file created ******/


	int    sc = 0;
	char   ch;
	FILE   *fp;

	fp = fopen("asm.inp", "r");
	if(!(fp)) {
		printf("Could not open asm.inp \n");
		exit(1);
	} /* if */
	
	fscanf(fp, "%*s");
	ch=getc(fp);
	while (ch!='\n') {
		if (!((ch==' ')&&(sc==0))) {
			outputdir[sc]=ch;
			sc++;
		}
		ch=getc(fp);
	}

	fscanf(fp,"%*s %i %*s %i", &NPI, &NPJ);
	fscanf(fp,"%*s %lf %*s %lf", &XMAX, &YMAX);
	fscanf(fp,"%*s %lf %*s %lf", &USLIP, &ALPHAIN);
	fscanf(fp,"%*s %lf %*s %lf", &IinLeft, &IinRight);
	fscanf(fp,"%*s %lf %*s %lf", &Dt, &TOTAL_TIME);
	fscanf(fp,"%*s %d", &TOTAL_DUMPS);
	fscanf(fp,"%*s %lf %*s %lf", &RHOG, &RHOL);
	fscanf(fp,"%*s %i", &NNOZ);
	fscanf(fp,"%*s %lf", &DNOZ);
    fscanf(fp,"%*s %lf", &V_IN);
	fscanf(fp,"%*s %i", &CONFIG);
	fclose(fp);
} /* readconstants */

/* ################################################################# */
void defaultconstants(void)
/* ################################################################# */
{
/***** Purpose: uses default constants instead of reading them from ******/
/***** the asm.inp file ******/

	sprintf (outputdir, "E:\\My Documents\\Onderwijs\\Icfd\\ASMtest\\calc1\\");
	NPI         = 50;
	NPJ         = 50;
	XMAX        = 0.5;
	YMAX        = 1.5;
	USLIP       = 0.2;
	ALPHAIN     = 0.05;
	IinLeft     = 0.66;
	IinRight    = 0.74;
	Dt          = 0.05;
	TOTAL_TIME  = 100.;
	TOTAL_DUMPS = 400;
	RHOG        = 1.2;
	RHOL        = 1000.;

} /* defaultconstants */

/* ################################################################# */
int *int_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type int */
	int *a;

	a = (int *) calloc(np, sizeof(int));

	return a;

} /* int_1D_array */

/* ################################################################# */
double *double_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type double */
	double *a;

	a = (double *) calloc(np, sizeof(double));

	return a;

} /* double_1D_array */

/* ################################################################# */
double **double_2D_matrix (int nm, int np)
/* ################################################################# */
{
/* create an 2D matrix with size [nm, np] of type double */
	int i;
	double **m;

	m = (double **) calloc(nm, sizeof(double *));
	for ( i = 0; i < nm; i++)
		m[i] = (double *) calloc(np, sizeof(double));

	return m;

} /* double_2D_matrix */

/* ################################################################# */
void memalloc(void)
/* ################################################################# */
{
	x    = double_1D_array(NPI + 2);
	x_u  = double_1D_array(NPI + 2);
	y    = double_1D_array(NPJ + 2);
	y_v  = double_1D_array(NPJ + 2);

	u      = double_2D_matrix(NPI + 2, NPJ + 2);
	v      = double_2D_matrix(NPI + 2, NPJ + 2);
	pc     = double_2D_matrix(NPI + 2, NPJ + 2);
	p      = double_2D_matrix(NPI + 2, NPJ + 2);
	T      = double_2D_matrix(NPI + 2, NPJ + 2);
	rho    = double_2D_matrix(NPI + 2, NPJ + 2);
	mu     = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma  = double_2D_matrix(NPI + 2, NPJ + 2);
	Cp     = double_2D_matrix(NPI + 2, NPJ + 2);
	Alpha  = double_2D_matrix(NPI + 2, NPJ + 2);

	u_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	v_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	pc_old = double_2D_matrix(NPI + 2, NPJ + 2);
	T_old  = double_2D_matrix(NPI + 2, NPJ + 2);
	rho_old= double_2D_matrix(NPI + 2, NPJ + 2);
	Alpha_old= double_2D_matrix(NPI + 2, NPJ + 2);

	usum  = double_2D_matrix(NPI + 2, NPJ + 2);
	vsum  = double_2D_matrix(NPI + 2, NPJ + 2);
	psum  = double_2D_matrix(NPI + 2, NPJ + 2);
	rhosum= double_2D_matrix(NPI + 2, NPJ + 2);
	musum = double_2D_matrix(NPI + 2, NPJ + 2);
	alphasum= double_2D_matrix(NPI + 2, NPJ + 2);

	umean  = double_2D_matrix(NPI + 2, NPJ + 2);
	vmean  = double_2D_matrix(NPI + 2, NPJ + 2);
	pmean  = double_2D_matrix(NPI + 2, NPJ + 2);
	rhomean= double_2D_matrix(NPI + 2, NPJ + 2);
	mumean = double_2D_matrix(NPI + 2, NPJ + 2);
	alphamean= double_2D_matrix(NPI + 2, NPJ + 2);

	dudx   = double_2D_matrix(NPI + 2, NPJ + 2);
	dudy   = double_2D_matrix(NPI + 2, NPJ + 2);
	dvdx   = double_2D_matrix(NPI + 2, NPJ + 2);
	dvdy   = double_2D_matrix(NPI + 2, NPJ + 2);

	aP     = double_2D_matrix(NPI + 2, NPJ + 2);
	aE     = double_2D_matrix(NPI + 2, NPJ + 2);
	aW     = double_2D_matrix(NPI + 2, NPJ + 2);
	aN     = double_2D_matrix(NPI + 2, NPJ + 2);
	aS     = double_2D_matrix(NPI + 2, NPJ + 2);
	b      = double_2D_matrix(NPI + 2, NPJ + 2);

	Istart = int_1D_array(NFIMAX + 1);
	Iend   = int_1D_array(NFIMAX + 1);
	Jstart = int_1D_array(NFIMAX + 1);
	Jend   = int_1D_array(NFIMAX + 1);

	SP     = double_2D_matrix(NPI + 2, NPJ + 2);
	Su     = double_2D_matrix(NPI + 2, NPJ + 2);

	F_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	F_v    = double_2D_matrix(NPI + 2, NPJ + 2);

	d_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	d_v    = double_2D_matrix(NPI + 2, NPJ + 2);

	relax  = double_1D_array(NFIMAX + 1);

	solve_fi = int_1D_array(NFIMAX + 1);
	print_fi = int_1D_array(NFIMAX + 1);


} /* memalloc */
