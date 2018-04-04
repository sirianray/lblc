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
int newrun_on, wall_x, wall_y, wall_z, flow_on, Q_on, debug_on, patch_on, phi_on;
int rand_init, rand_seed;
int t_current, t_max, t_print, t_write;
int n_evol_Q=1, uconverge=1, qconverge=1, pconverge=1, n_pre_evol=10000, n_evol_phi=1;
int bulk0=0;
int type_xlo, type_xhi, type_ylo, type_yhi, type_bot, type_top;
double dt=1.0, qdt, pdt, tau_f, itau_f, L1, L2, L3, L4, xi, xi1, xi2, Gamma_rot, rho, S_lc, S_lc2, A_ldg, Temp=third, e[15][3], q_init;
double q_ch, twqL_ch;
//double ux_top, ux_bot, uy_top, uy_bot, xforce, yforce, zforce;
double ux_lo, uy_lo, uz_lo, ux_hi, uy_hi, uz_hi, xforce, yforce, zforce;
double Q_tol=5e-28, u_tol=2e-15, k_eng=-1, Q_diff=100, e_toto=0., e_tot, e_ld, e_el, e_ch, e_sf, e_L1, e_L2, e_L3, e_L4, e_ly;
double W_xlo, W_xhi, W_ylo, W_yhi, W_top, W_bot, n_xlo[3], n_xhi[3], n_ylo[3], n_yhi[3], n_top[3], n_bot[3];
double K1, K2, K3, K4, K24;
double zeta;

// lyotropic related
double A_phi, ly_k, ly_Gamma, ly_e0, phi_diff=1., ly_an, ly_phi1, ly_phi2;
double ldg_a, ldg_b, ldg_c;
double ly_phi_tot, ly_phi_tot0;
real *ly_phi, *ly_mu, *ly_dphi;
MPI_Win winly_phi, winly_mu;

int *neighb, *nextf, *info, *neighbsurf;
real *Q, *H, *surf, *Qsurf, *Hsurf;
real *Rho, *u, *W, *f, *p, *f2, *Cf, *sigma_q, *sigma_p;

int myid=0, numprocs=1;
MPI_Win winq, wins, winr, winu, winf, winf2, winp, winQsurf, winHsurf, winneighbsurf, winsurf, winnf, wininfo, winneighb, winsigma_p;
MPI_Comm shmcomm;
int point, lpoint, qpoint, node, nodes;

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	double t_begin, zeta0;
	int i, j, k, ii, id, id1;

	t_begin = MPI_Wtime();		// record the begining CPU time
	
	read_param();
	lattice_vec();
	allocate();
	p_allocate();

	if (flow_on!=0) {
		build_stream();
	}

	if (Q_on!=0 || phi_on!=0) {
		build_neighbor();		
	}

	init_surf();
	add_patch();
	p_init();
	p_iden();
    p_set_neighb();
	init();

    zeta0 = zeta;
    zeta  = 0.;

	if (flow_on!=0) cal_fequ(f);	
	if (Q_on!=0) cal_dQ();
    if (phi_on!=0) cal_mu();

	if ((Q_on!=0 || phi_on!=0) && flow_on!=0 && newrun_on!=0) {
		while( (qconverge==0 || pconverge==0) && t_current<n_pre_evol) {
			t_current++;
            if (Q_on!=0) {
		        for (ii=0; ii<n_evol_Q; ii++) {
                    cal_dQ();
                    evol_Q();
                }
            }
            if (phi_on!=0) {
		        for (ii=0; ii<n_evol_phi; ii++) {
                    cal_mu();
                    evol_phi();
                }
            }
			if (t_current%t_print==0) monitor();
		}
		e_tot     =-1;
		k_eng     =-1;
		qconverge = 0;
        pconverge = 0;
		uconverge = 0;
		t_current =-1;
	}

	if (Q_on!=0 && flow_on!=0) {
		cal_W();
		cal_stress2();
		cal_sigma_p();
	}
    if (phi_on!=0 && flow_on!=0) {
        cal_mu();
        cal_stress_phi(Q_on);
    }
    MPI_Barrier(MPI_COMM_WORLD);

	output1(1,'z',Nx/2,Ny/2);
	output3(1);
	if(myid==0) printf("Q initialized\n");
	MPI_Barrier(MPI_COMM_WORLD);

    zeta = zeta0;

	if (t_current%t_print==0) monitor();
	while (t_current<t_max && uconverge*qconverge*pconverge==0) {
		e_toto=e_tot;
		if (Q_on!=0 && qconverge==0) {
			if (flow_on!=0 && uconverge==0) cal_W();
			for (ii=0; ii<n_evol_Q; ii++) {
				cal_dQ();
				evol_Q();
			}			
		}

        if (phi_on!=0 && pconverge==0) {
			if (flow_on!=0 && uconverge==0) cal_W();
			for (ii=0; ii<n_evol_phi; ii++) {
                cal_mu();
                evol_phi();
			}			
        }

		if (flow_on!=0 && uconverge==0) {
			if (Q_on!=0 && qconverge==0) {
				cal_stress2();
				cal_sigma_p();
			}
			if (phi_on!=0 && pconverge==0) {
                cal_mu();
				cal_stress_phi(Q_on);
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

        if (phi_on!=0 && pconverge==0) {
			if (flow_on!=0 && uconverge==0) cal_W();
			for (ii=0; ii<n_evol_phi; ii++) {
                cal_mu();
                evol_phi();
			}			
        }

		if (flow_on!=0 && uconverge==0) {
			if (Q_on!=0 && qconverge==0) {
				cal_stress2();
				cal_sigma_p();
			}
			if (phi_on!=0 && pconverge==0) {
                cal_mu();
				cal_stress_phi(Q_on);
			}
			evol_f(f2,f);
		}
		
		if (t_current%t_print==0) monitor();
		if (t_current%t_write==0 && t_current>0) {
			output1(0,'z',Nx/2,Ny/2);
			output3(0);
	        write_restart();
			fflush(stdout);
		}
		t_current++;
	}
	
	
	output_time(t_begin);
	output1(0,'z',Nx/2,Ny/2);
//	output3(1);
	
	write_restart();

	p_deallocate();
	deallocate();

	return 0;
}
