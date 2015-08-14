#include "lb.h"
#include "particle.h"

int e[15][3], t_print, t_max, t_current=0, n_evol_Q=1;
int newrun_on, particle_on, wall_on, wall_anchoring_on, flow_on, Q_on, simul_evol_on, debug_on;
real dt=1.0, tau_f, itau_f, kappa, T=third, xi, Gamma_rot, rho, yforce, U, A_ldg;
real uy_top, uy_bottom, S;
real H[Nx][Ny][Nz][3][3], convQ[Nx][Ny][Nz][5];
real sigma[Nx][Ny][Nz][3][3], u[Nx][Ny][Nz][3], Rho[Nx][Ny][Nz], W[Nx][Ny][Nz][3][3];
real f_eq[Nx][Ny][Nz][15], Cf[Nx][Ny][Nz][15], p[Nx][Ny][Nz][15];
//real Qo_top[5], Qo_bottom[5];
real f1[Nx][Ny][Nz][15], f2[Nx][Ny][Nz][15];
real Q1[Nx][Ny][Nz][3][3], Q2[Nx][Ny][Nz][3][3];
real k_eng=-1, sigma_q[Nx][Ny][Nz][3][3], sigma_p[Nx][Ny][Nz][3], lap_Q[Nx][Ny][Nz];
real qthreshold=5e-25, qdt;

//Modifications for finite anchoring
real Qsurf[Nx][Ny][2][5], Qsurf0[Nx][Ny][2][5];
real dQdxnu[Nx][Ny][2][5];
real W_top, W_bot, ntop[3], nbot[3], ninit[3];
int rand_init, rand_seed;

int main()
{
	int t, i, j, k, ii, jj, kk;
	FILE *out1, *out2, *out3, *out4;
	int time_converge, u_run, q_run;

	read_param();
	lattice_vec();
	if (particle_on!=0) {
		p_init();
		p_iden();
	}

	init(f1, Q1);
	u_run = flow_on;
	q_run = Q_on;
	
//	writing to data files
	out1=fopen("density","w");
	out2=fopen("ux","w");
	out3=fopen("uy","w");
	out4=fopen("uz","w");
	monitor(0,1);
	print_Q1d(0);
	print_u1d(out1,out2,out3,out4);
	print_visual(1);
	
//	evolution
	if (simul_evol_on==0) {	//	alternating evolution
		while ( (Q_on!=0 || flow_on!=0) && t_current<t_max ) {
			time_converge = 0;
			while ( u_run!=0 && flow_on!=0 && t_current<t_max) {
				t_current++;
				time_converge++;
				cal_sigma(Q1);
				evol_p(Q1);
				evol_f(f1,f2,Q1);
				cal_sigma(Q1);
				evol_p(Q1);
				evol_f(f2,f1,Q1);
				if(t_current%t_print==0){
					printf("t=%d:\n",t_current);
					monitor(1,1);
					print_u1d(out1, out2, out3, out4);
					print_Q1d(1);
					print_visual(0);
				}
				else {
					monitor(1,0);
				}
			}
			if ( u_run!=0 && (time_converge/t_print>1 || t_current<5)) {
				flow_on = 1;
				cal_W(u);
			}
			
			time_converge=0;
			while(q_run!=0 && Q_on!=0 && t_current<t_max){
				t_current++;
				evol_Q(Q1,Q2);
				evol_Q(Q2,Q1);
				time_converge++;
				if(t_current%t_print==0){
					printf("t=%d:\n",t_current);
					monitor(-1,1);
					print_u1d(out1, out2, out3, out4);
					print_Q1d(1);
					print_visual(0);
				}
				else {
					monitor(-1,0);
				}
				
			}
			if(q_run!=0 && time_converge/t_print>1){
				Q_on = 1;
				cal_stress_q(Q1,H,sigma_q);
				cal_sigma_p(Q1);
			}
		}
	} else {			// simultaneous evolution
//		equilibrate Q tensor before simultaneous evolution
		while ( (Q_on!=0 || flow_on!=0) && t_current<t_max ) {
			t_current++;
			cal_sigma(Q1);
			if (particle_on!=0) {
			  evol_p(Q1);
			}
			if (u_run!=0) {
			  evol_f(f1,f2,Q1);
			}
			if (q_run!=0) {
			  cal_W(u);
			for(i=0;i<n_evol_Q;i++){
			  evol_Q(Q1,Q2);
			  evol_Q(Q2,Q1);
                        }
			  if(u_run!=0) {
			    cal_stress_q(Q1,H,sigma_q);
	  		    cal_sigma_p(Q1);
        }
			}
			cal_sigma(Q2);
			if (particle_on!=0) {
			  evol_p(Q1);
			}
			if (u_run!=0) {
			  evol_f(f2,f1,Q2);
			}
			if (q_run!=0) {			
			  cal_W(u);
			for(i=0;i<n_evol_Q;i++){
                          evol_Q(Q1,Q2);
                          evol_Q(Q2,Q1);
                        }
			  if(u_run!=0) {
 			  cal_stress_q(Q1,H,sigma_q);
   			cal_sigma_p(Q1);
    	  }
			}
			
			if(t_current%t_print==0){
				printf("t=%d:\n",t_current);
				monitor(0,1);
				print_u1d(out1, out2, out3, out4);
				print_Q1d(1);
				print_visual(0);
			}
			else {
				monitor(0,0);
			}
		}
	}
	
//	print_visual();
	
	print_Q2d(0,Nx/2,0);
	print_Q2d(1,Ny/2,1);
//	print_Q2d(2,(Nz+1)/2,1);
	print_u2d(0, Nx/2, 0);
	write_restart(Q1,Rho,u);
	fclose(out1);
	fclose(out2);
	fclose(out3);
	fclose(out4);

	return 0;
}
