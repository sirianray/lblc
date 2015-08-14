/*
 This is Lattice Boltzmann method for liquid crystal hydrodynamic system in a
 channel with particles
 
 newrun_on:			new run if !=0
 particle_on:		particle on if !=0
 wall_on:			periodic boundary condition in z if ==0; bounce back wall in z if !=0
 wall_anchoring_on:	apply nematic anchoring on two z walls if !=0 (still work even if wall_on==0)
 Q_on:				evolve Q if !=0
 flow_on:			evolve f's if !=0
 simul_evol_on:		simultaneously evolve Q and f's if !=0; alternatively evolve Q and f's if ==0
 debug_on:			debug mode on if !=0
 yforce:			body force in y direction (not applied to the particle interior)
 uy_top:			moving top wall at constant velocity
 uy_bottom:			moving bottom wall at constant velocity
 n_evol_Q:		How many times Q tensor evolve when velocity evolves once
 qdt:			delt_t for Q evolution = 1/(2*n_evol_Q)
 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef double real;
#define Nx 41
#define Ny 41
#define Nz 41
#define third   0.3333333333333333
#define one24th 0.0416666666666667
#define one12th 0.0833333333333333
#define one30th 0.0333333333333333
#define two3rd  0.6666666666666667
#define four3rd 1.3333333333333333
#define eight3rd 2.666666666666667

extern int e[15][3], t_print, t_max, t_current, n_evol_Q;
extern int newrun_on, particle_on, wall_on, wall_anchoring_on, flow_on, Q_on, simul_evol_on, debug_on;
extern real dt, tau_f, itau_f, kappa, T, xi, Gamma_rot, rho, yforce, U, A_ldg;
extern real uy_top, uy_bottom, S;
extern real H[Nx][Ny][Nz][3][3], convQ[Nx][Ny][Nz][5];
extern real sigma[Nx][Ny][Nz][3][3], u[Nx][Ny][Nz][3], Rho[Nx][Ny][Nz], W[Nx][Ny][Nz][3][3];
extern real f_eq[Nx][Ny][Nz][15], Cf[Nx][Ny][Nz][15], p[Nx][Ny][Nz][15];
//extern real Qo_top[5], Qo_bottom[5];
extern real f1[Nx][Ny][Nz][15], f2[Nx][Ny][Nz][15];
extern real Q1[Nx][Ny][Nz][3][3], Q2[Nx][Ny][Nz][3][3];
extern real k_eng, sigma_q[Nx][Ny][Nz][3][3], sigma_p[Nx][Ny][Nz][3], lap_Q[Nx][Ny][Nz];
extern real qthreshold, qdt;

// function declaration:

void read_param();
void lattice_vec();
void init(real f[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3]);
void evol_f(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3]);
void evol_p(real Q0[Nx][Ny][Nz][3][3]);
void streaming(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15]);
//void bc(real f1[Nx][Ny][Nz][15]);

//Modifications for finite anchoring
extern real Qsurf[Nx][Ny][2][5], Qsurf0[Nx][Ny][2][5];
extern real dQdxnu[Nx][Ny][2][5];
extern real W_top, W_bot, ntop[3], nbot[3], ninit[3];
extern int rand_init, rand_seed;
void surface_derivative(real Q0[Nx][Ny][Nz][3][3]);
real randvec();


void cal_fequ(real f_eq[Nx][Ny][Nz][15], real Rho0[Nx][Ny][Nz], real u0[Nx][Ny][Nz][3], real sigma0[Nx][Ny][Nz][3][3]);
// cal_feqf is the old subroutine, but has best performance for incompressibile flow
// cal_feqf_new is faster and supposed to be more accurate, 
// good for compressible flow but not the first choice for incompressible flow
//void cal_feqf(real f_eq[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real sigma0[Nx][Ny][Nz][3][3]);
void cal_feqf_new(real f_eq0[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real sigma[Nx][Ny][Nz][3][3]);
void cal_sigma(real Q0[Nx][Ny][Nz][3][3]);
void cal_p(real p[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real feq[Nx][Ny][Nz][15]);
void cal_W(real u0[Nx][Ny][Nz][3]);
void cal_dQ(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real convQ[Nx][Ny][Nz][5]);
void cal_stress_q(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real sigma_q[Nx][Ny][Nz][3][3]);
void cal_sigma_p(real Q0[Nx][Ny][Nz][3][3]);

void cal_u(real u0[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15]);
void cal_rho(real rho_local[Nx][Ny][Nz], real f[Nx][Ny][Nz][15]);
void cal_mom(real m[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15]);
											// calculate momentum from distribution function f's
void monitor(int flag, int print);			// flag>0: monitor flow; flag<0: monitor Q; flag=0: monitor both

// for data output
void print_Q1do(int flag);					// old version of print_Q1d, may be deleted in the future
void print_Q1d(int flag);					// 1d data of Q to file: q1, q2, ..., q5. They can be processed by measS.m
											// flag=1: append; flag=0: new file
void print_Q2d(int dim, int x, int flag);	// 2d data of Q to file: Q1, Q2, ..., Q5. They can be processed by measSQ.m
void print_u1d(FILE *FILE_rho, FILE *FILE_ux, FILE *FILE_uy, FILE *FILE_uz);
											// 1d data of hydrodynamics to files. Directly visulizable
void print_u2d(int dim, int x, int flag);	// 2d data of hydrodynamics to file: cuy, cuz, crho. Directly visulizable
void print_visual(int new);						// 3d data to file: grid, rho_3d, u_3d, Q_3d, stress_3d
											// They can be processed by iso from postprocess folder

// for restart
void write_restart(real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3]);
void read_restart (real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3]);

// for functionality
real getS();								// calculate order parameter S from U
real Delta(int i, int j);					// Delta function I
real OneThirdDelta(int i, int j);			// I/3
void clear2d(real A[3][3]);					// clear 2d array A[3][3] (reset to 0)
real QQ(real T[3][3], int i, int j);		// sum_k{T[i][k]*T[k][j]}
int bounce(int i);							// opposite direction of i on the lattice
void period_image(int *i, int *j, int *k);	// periodic image of point [i, j, k] in main box
