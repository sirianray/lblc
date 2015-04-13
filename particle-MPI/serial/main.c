/*
 *  main.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */

#include "main.h"
#include "particle.h"

int Nx, Ny, Nz, points, npar, nsurf, qpoints;
int newrun_on, wall_on, flow_on, Q_on, debug_on;
int rand_init, rand_seed;
int t_current, t_max, t_print, t_write;
int n_evol_Q=1, uconverge=1, qconverge=1;
int bulk0=0, qlength;
int type_bot, type_top;
double dt=1.0, qdt, tau_f, itau_f, kappa, xi, xi1, xi2, Gamma_rot, rho, U_lc, S_lc, A_ldg, Temp=third, e[15][3];
double yforce, uy_top, uy_bottom;
double qthreshold=1e-15, uthreshold, k_eng=-1, Q_diff=100, e_toto=0., e_tot, e_ld, e_el, e_sf;
double W_top, W_bot, n_top[3], n_bot[3];

int *neighbor, *nextf, *info;
real *Q, *H, *surf, *vsurf;
real *Rho, *u, *W, *f, *p, *f2, *Cf, *sigma_q, *sigma_p;

int myid, numprocs;

int main(int argc, char *argv[]){
//	MPI_Init(&argc, &argv);
//	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
//	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	double e_toto, err;
	int converge=0, i, j, k, ii, id, id1;
	
	read_param();
	lattice_vec();
	allocate();
	p_allocate();
	
	if (flow_on!=0) {
		build_stream();
	}
	if (Q_on!=0) {
		build_neighbor();		
	}
	
	init_surf();
	p_init();
	p_iden();
	init();
	
	
	if (flow_on!=0) cal_fequ(f);	
	if(Q_on!=0) cal_dQ();
	
	if (Q_on!=0 && flow_on!=0 && npar>0) {
//		while (qconverge==0) {
//			t_current++;
//			cal_dQ();
//			evol_Q();
//			if (t_current%t_print==0 && myid==0) monitor();
//		}
		uconverge = 0;
		qconverge = 0;
	}
	
	if (Q_on!=0 && flow_on!=0) {
		cal_W();
		cal_stress();
		cal_sigma_p();
	}
	
	output1(1,'z',Nx/2,Ny/2);
	output3(1);
	monitor();
		
	while (t_current<t_max && uconverge*qconverge==0) {
		e_toto=e_tot;
		if (Q_on!=0 && qconverge==0) {
			if (flow_on!=0 && uconverge==0) cal_W();
			for (ii=0; ii<n_evol_Q; ii++) {
				cal_dQ();
				evol_Q();
			}			
		}

		if (flow_on!=0 && uconverge==0) {
			if (Q_on!=0 && qconverge==0) {
				cal_stress();
				cal_sigma_p();
			}
			evol_f(f,f2);
		}

		if (Q_on!=0 && qconverge==0) {
			if (flow_on!=0 && uconverge==0) cal_W();
			for (ii=0; ii<n_evol_Q; ii++) {
				cal_dQ();
				evol_Q();
			}			
		}

		if (flow_on!=0 && uconverge==0) {
			if (Q_on!=0 && qconverge==0) {
				cal_stress();
				cal_sigma_p();
			}
			evol_f(f2,f);
		}
		
		if (t_current%t_print==0 && myid==0) monitor();
		if (t_current%t_write==0 && myid==0) {
			output1(0,'z',Nx/2,Ny/2);
			output3(0);
		}
		
		t_current++;
	}
	
	
	output1(0,'z',Nx/2,Ny/2);
	output3(0);

//	printf("this is processor %d out of %d processors\n",myid,numprocs);
	printf("Nx=%d, Temp=%f\n",Nx,Temp);
	
	if (myid==0) write_restart();
	
	deallocate();
	p_deallocate();
//	MPI_Win_free(&win);
//	MPI_Comm_free(&shmcomm);
//	MPI_Finalize();
	
	return 0;
}