/*
	This is Lattice Boltzmann method for liquid crystal hydrodynamic system in a
channel

pbc_z:		    periodic boundary condition in z direction (z-wall existence)
wall_anchoring: apply anchoring on the two walls (anchoring can still work if pbc_z is on)
q_converge:     not evolve Q if it is nonzero (= 1 - q_run)
u_converge:     not evolve flow if it is nozero (= 1 - flow_on)
t_v_on:         time when Q tensor kicks in the evolution of flow
yforce:         body force in y direction (not applied to the particle interior)
uy_top:         moving top wall at constant velocity
uy_bottom:      moving bottom wall at constant velocity

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double real;
#define Nx 20
#define Ny 20
#define Nz 20

#define third   0.3333333333333333
#define one24th 0.0416666666666667
#define one12th 0.0833333333333333
#define one30th 0.0333333333333333
#define two3rd  0.6666666666666667
#define four3rd 1.3333333333333333
#define eight3rd 2.666666666666667

extern int e[15][3], t_print, t_max, movie;
extern int pbc_z, wall_anchoring, Q_converge, u_converge, debug_mode, simul_evol, Q_on, u_on;
extern int t_current, newrun, t_v_on, init_linear;
extern real dt, tau_f, tau_G, kappa, T, xi, Gamma_rot, rho, yforce, U, A_ldg;
extern real uy_top, uy_bottom;
extern real itau_f, itau_G;
extern real H[Nx][Ny][Nz][3][3], convQ[Nx][Ny][Nz][5];
extern real sigma[Nx][Ny][Nz][3][3], u[Nx][Ny][Nz][3], Po[Nx][Ny][Nz];
extern real W[Nx][Ny][Nz][3][3];
extern real S, Rho[Nx][Ny][Nz];
extern real f_eq[Nx][Ny][Nz][15], Cf[Nx][Ny][Nz][15], p[Nx][Ny][Nz][15];
extern real Qo_top[5], Qo_bottom[5];
extern real f1[Nx][Ny][Nz][15], f2[Nx][Ny][Nz][15];
extern real Q1[Nx][Ny][Nz][3][3], Q2[Nx][Ny][Nz][3][3];
extern real k_eng, sigma_q[Nx][Ny][Nz][3][3], tau[Nx][Ny][Nz][3][3], detau[Nx][Ny][Nz][3], lap_Q[Nx][Ny][Nz];
extern real mcheck[Nx][Ny][Nz][3];

//Modifications for finite anchoring
extern real Qsurf[Nx][Ny][2][5];
extern real dQdxnu[Nx][Ny][2][5];
extern real W_top, W_bot, ntop[3], nbot[3], ninit[3];
extern int rand_init, rand_seed;
void surface_derivative(real Q0[Nx][Ny][Nz][3][3]);
real randvec();

//Cavity modifications
extern int cav_height, cav_width, cav_flag[Nx][Ny][Nz], cav_on;
void cav_surf();
extern real Qtop[3][3], Qsides[3][3];
extern real ncavtop[3], ncavsides[3];

//Channel parameters
void read_param();
void lattice_vec();
void init(real f[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3]);
void evol_f(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15]);
void streaming(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15]);
//void bc(real f1[Nx][Ny][Nz][15]);

real getS();
real Delta(int i, int j);
real OneThirdDelta(int i, int j);
void clear2d(real A[3][3]);
void cal_fequ(real f_eq[Nx][Ny][Nz][15], real Rho0[Nx][Ny][Nz], real u0[Nx][Ny][Nz][3], real sigma0[Nx][Ny][Nz][3][3]);
// cal_feqf is the old subroutine, but has best performance for incompressibile flow
// cal_feqf_new is faster and supposed to be more accurate, good for compressible flow but not the first choice for incompressible flow
void cal_feqf(real f_eq[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real sigma0[Nx][Ny][Nz][3][3]);
void cal_feqf_new(real f_eq0[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real sigma[Nx][Ny][Nz][3][3]);
void cal_sigma(real Q0[Nx][Ny][Nz][3][3]);
void cal_p(real p[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real feq[Nx][Ny][Nz][15]);
void cal_W(real u0[Nx][Ny][Nz][3]);
void cal_dQ(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real Po[Nx][Ny][Nz], real convQ[Nx][Ny][Nz][5]);
void cal_stress_q(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real sigma_q[Nx][Ny][Nz][3][3]);

void cal_u(real u0[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15]);
void cal_rho(real rho_local[Nx][Ny][Nz], real f[Nx][Ny][Nz][15]);
void cal_mom(real m[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15]);
void monitor(int flag, int print);	// positive: monitor flow; negative: monitor Q; 0: monitor both
int bounce(int i);
real QQ(real T[3][3], int i, int j);
void print_Q(int flag);
void print_Q2(int flag);
void print_Q3(int dimension, int x, int flag);
void print_u(int dimension, int x, int flag);
void print_u2(FILE *FILE_rho, FILE *FILE_ux, FILE *FILE_uy, FILE *FILE_uz);
void write_restart(real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3]);
void read_restart (real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3]);

void cal_sigma_p(real Q0[Nx][Ny][Nz][3][3]);
void print_visual(int input);
void check(real Q0[Nx][Ny][Nz][3][3]);
