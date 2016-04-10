/*
 *  main.h
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 *	Variable Explanation:
 Nx Ny Nz:		Lattice size
 points:		Nx*Ny*Nz, number of bulk points (assuming no particle)
 npar:			number of particles (if 0, no extra surface point except wall will be added)
 nsurf:			number of surface points: 0 - no wall & no particle; 2*Nx*Ny - wall but no particle
 qpoints:		number of points having Q = points+nsurf: points - no wall & no particle;	Nx*Ny*(Nz+2) - wall but no particle
 newrun_on:		1 - read initial configurations (Q, Rho, u) from restart.dat; 0 - start new run
 wall_on:		1 - no-slip wall & anchoring wall will be present in z direction; 0 - Periodic b.c.
 rand_init:		whether to randomly initialize Q: -1 - from a best guess; 0 - random seed; >0: used rand_seed as seed
 t_print:		every that steps to print to the screen & check convergence
 t_write:		every that steps to write to files
 n_evol_Q:		control qdt for the time step of Q evolution: qdt = 1/n_evol_Q
 bulk0:			0 - no wall; Nx*Ny - wall
 type_*:		type of anchoring for *(top/bottom) wall: 0 - infinite anchoring; 1 - degenerate planar; 2 - finite nondegenerate
 W_*:			anchoring strength of *(top/bottom) wall
 n_*[]:			prefered orientation for *(top/bottom) wall
 
 neighbor[0-5]:	for bulk point is 6 neighbors in 6 directions:	+; neighbor is bulk (arm is 1);
																-: neighbor is surface (arm is 0.5)
				for surface point is 3 neighbors and 3 next neighbors in 3 directions
																+ +: arms 0.5 and 1.5
																- +: arms 1   and 1
																- -: arms 1   and 2
				(Exists if Q_on!=0)
 
 "direction":	0: -x; 1: +x; 2: -y; 3: +y; 4: -z; 5: +z
 
 nextf[0-14]:	the position in f after streaming; bounce back on wall/particle is accounted for.
				(Exists if flow_on!=0)
 
 info[0]:		type of Q-point: -ipar-2 - inside ipar'th particle; -1 - bulk point; >=0: the order in surf or vsurf
				(Always exists)
 
 surf[0-9]:		info of the surface: anchoring type, anchoring strength, normal (nx, ny, nz), prefered Q (5 elements).
				(Exists if Q_on!=0)
 
 vsurf[0-2]:	velocity of the surface point.
				(Exists if Q_on!=0 && flow_on!=0)
 
 W[i,j]:		dj_ui =  ( 0 3 6 )
						 ( 1 4 7 )
						 ( 2 5 8 )
				(Exists if Q_on!=0 && flow_on!=0)
 
 sigma_q[0-9]:	stress due to LC phase
				(Exists if Q_on!=0 && flow_on!=0)
 
 sigma_p[0-2]:	forcing in Navier-Stokes; derivative of sigma_q
				(Exists if Q_on!=0 && flow_on!=0)
 
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef double real;
#define third 0.3333333333333333
#define one24th 0.0416666666666667
#define one12th 0.0833333333333333
#define one30th 0.0333333333333333
#define two3rd  0.6666666666666667
#define four3rd 1.3333333333333333
#define eight3rd 2.666666666666667
#define PI 3.141592653589793

#define root 0

extern int Nx, Ny, Nz, points, npar, nsurf, qpoints;
extern int newrun_on, wall_x, wall_y, wall_z, flow_on, Q_on, debug_on, patch_on;
extern int rand_init, rand_seed;
extern int t_current, t_max, t_print, t_write;
extern int n_evol_Q,  uconverge, qconverge;
extern int bulk0;
extern int type_xlo, type_xhi, type_ylo, type_yhi, type_bot, type_top;
extern double dt, qdt, tau_f, itau_f, kappa, xi, xi1, xi2, Gamma_rot, rho, U_lc, S_lc, A_ldg, Temp, e[15][3], q_init;
extern double q_ch, twqL_ch;
//extern double ux_top, ux_bot, uy_top, uy_bot, xforce, yforce, zforce;
extern double ux_lo, uy_lo, uz_lo, ux_hi, uy_hi, uz_hi, xforce, yforce, zforce; 
extern double Q_tol, u_tol, k_eng, Q_diff, e_toto, e_tot, e_ld, e_el, e_ch, e_sf, Fld0;
extern double W_xlo, W_xhi, W_ylo, W_yhi, W_top, W_bot, n_xlo[3], n_xhi[3], n_ylo[3], n_yhi[3], n_top[3], n_bot[3];

extern int *neighb, *nextf, *info, *neighbsurf;
extern real *Q, *H, *surf, *Qsurf, *Hsurf;
extern real *Rho, *u, *W, *f, *p, *f2, *Cf, *sigma_q, *sigma_p;

extern int myid, numprocs;
extern MPI_Win winq, wins, winr, winu, winf, winf2, winp, winQsurf, winHsurf, winneighbsurf, winsurf, winnf;
extern MPI_Comm shmcomm;
extern int point, lpoint, qpoint, node, nodes;

void read_param();
void ntoq(double nx, double ny, double nz, double *q0, double *q1, double *q2, double *q3, double *q4, double *q5);
real randvec();
void cal_dQ();
void evol_Q();
void init_surf();
void build_neighbor();
double getS(double U);
double getF0(double U);
double trQQ(real *q);
double QQ(real *q, int k);
double OneThirdDelta(int i);
double One3rdDelta(int i, int j);
double Delta(int i, int j);
double QQQ(real *q);
void output3(int new);
void output1(int new, char d, int cx, int cy);
void output_time(double t_begin);
void lattice_vec();
void cal_rho(real *f0);
void cal_u(real *f0);
void cal_fequ(real *f0);
void cal_feqf_new(real *f0);
void cal_feqf_diff(real *fin, real *fout);
//void cal_sigma();
void cal_p();
void streaming(real *fin, real *fout);
int bounce(int i);
void build_stream();
void evol_f(real *fin, real *fout);
void cal_stress();
void cal_sigma_p();
void monitor();
void write_restart();
void read_restart();
void normalize(double *x, double *y, double *z);
void add_patch();
