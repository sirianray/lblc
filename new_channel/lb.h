/*
	This is Lattice Boltzmann method for liquid crystal hydrodynamic system in a
channel

q_converge:     not evolve Q if it is nonzero (= 1 - q_run)
u_converge:     not evolve flow if it is nozero (= 1 - flow_on)
yforce:         body force in y direction (not applied to the particle interior)
uy_top:         moving top wall at constant velocity
uy_bottom:      moving bottom wall at constant velocity
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define third    0.3333333333333333
#define one24th  0.0416666666666667
#define one12th  0.0833333333333333
#define one30th  0.0333333333333333
#define two3rd   0.6666666666666667
#define four3rd  1.3333333333333333
#define eight3rd 2.6666666666666667

extern int e[15][3], t_print, t_max, movie, Nx, Ny, Nz, wall_degen, cav_degen;
extern int Q_converge, u_converge, simul_evol, Q_on, u_on, cav_inf;
extern int t_current, newrun, wall_inf, cav_degen;
extern double dt, tau_f, tau_G, kappa, T, xi, Gamma_rot, rho, yforce, U, A_ldg, W_cav;
extern double uy_top, uy_bottom;
extern double itau_f, itau_G;


//Modifications for 1D storage
extern double *Q1, *Q2, *u, *f1, *f2, *H, *convQ, *sigma, *W, *Rho, *dQdxnu_top, *dQdxnu_bot;
extern double *u, *dtau, *f_eq, *Cf, *p, *tau, *Qsurftop, *Qsurfbot, *sigma_q;
extern double Qo_top[5], Qo_bottom[5];
extern void alloc();
extern int *cav_flag;
void free_all();
void degenerate_anchoring();
void degenerate_cavity();

//Arrays
extern double S;
extern double k_eng;

//Modifications for finite anchoring
extern double W_wall, ntop[3], nbot[3];
extern int rand_init, rand_seed;
void surface_derivative(double *Q0);
double randvec();

//Cavity modifications
extern int cav_height, cav_width, cav_on;
void cav_surf();
extern double Qtop[5], Qsides[5];
extern double ncavtop[3], ncavsides[3];

//Channel parameters
void read_param();
void lattice_vec();
void init(double *f, double *Q0);
void evol_f(double *f1, double *f2);
void evol_Q(double *Qin, double *Qout);
void streaming(double *f1, double *f2);
void cal_dtau(double *Q0);

double getS();
double Delta(int i, int j);
double OneThirdDelta(int i, int j);
void cal_fequ(double *f, double *u0);
void cal_feqf_new();
void cal_p(double *f);
void cal_W();
void cal_dQ(double *Q0);
void cal_stress(double *Q0);
void cal_sigma();
void cal_dtau(double *Q0);
void derivative(int i, int j, int k, double *Q0, double *d2, double *ddQ);

void cal_u(double *f);
void cal_rho(double *f);
void cal_mom(double *m, double *f);
void monitor(int flag, int print);	// positive: monitor flow; negative: monitor Q; 0: monitor both
int bounce(int i);
double QQ(double M[3][3], int i, int j);
void write_restart(double *Q, double *Rho, double *u);
void read_restart (double *Q, double *Rho, double *u);

void print_visual(int input);
void check(double *Q0);
