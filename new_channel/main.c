#include "lb.h"

int e[15][3], t_print, t_max, movie, Nx, Ny, Nz, wall_degen, cav_degen;
int  Q_on, u_converge, simul_evol, u_on, Q_converge, cav_inf;
int t_current=0, newrun, wall_inf;
double dt=1.0, tau_f, tau_G, kappa, T=third, xi, Gamma_rot, rho, yforce, U, A_ldg, *dQdxnu_top, *dQdxnu_bot;
double uy_top, uy_bottom;
double itau_f, itau_G;
double *H, *convQ, *sigma, *u, *W, *sigma_q, W_cav;
double *f_eq, *Cf, *p, *f1, *f2, *Q1, *Q2;
double Qo_top[5], Qo_bottom[5];
double k_eng=1, *tau, *dtau, *lap_Q;
double S, *Rho;

//Finite anchoring
double *Qsurfbot, *Qsurftop; 
double W_wall, nbot[3], ntop[3];
int rand_seed, rand_init;

//Cavity-Channel parameters
int cav_height, cav_width, *cav_flag, cav_on;
double Qtop[5], Qsides[5];
double ncavtop[3], ncavsides[3];

int main()
{
		  int t,i,j,k,ii,jj,kk;
		  int time_converge;

		  Q_converge = 0;
		  u_converge = 0;

		  read_param();
		  lattice_vec();
		  alloc();

		  //discover cavity points
		  cav_surf();

		  if(newrun) init(f1, Q1);

		  //equilibrate Q-field 
		  if(Q_on && rand_init) {
					 printf("\nBeginning initial equilibration\n");

					 while(!Q_converge){
								evol_Q(Q1,Q2);
								evol_Q(Q2,Q1);
								monitor(-1,0);
					 }

					 printf("Finished initial equilbration\n");
					 Q_converge = 0;
		  }

		  if(Q_on){
					 printf("Calculating initial stress\n");
					 //cal_stress(Q1);
					 cal_dtau(Q1);
		  }

		  //	evolution
		  if (!simul_evol) {	//	alternating evolution
					 printf("Alternating evolution\n\n");
					 while( ( (!Q_converge && Q_on) || (!u_converge && u_on) ) && t_current<t_max){

								time_converge = 0;

								while(u_on && !u_converge && t_current<t_max){

										  t_current++;

										  cal_sigma();
										  evol_f(f1,f2);
										  cal_sigma();
										  evol_f(f2,f1);

										  time_converge++;

										  if(t_current%t_print==0){

													 printf("\nt=%d:\n",t_current);
													 monitor(1,1);
													 if(movie) print_visual(0);
													 print_visual(1);
										  }
										  else {
													 monitor(1,0);
										  }
								}

								if(u_on) printf("\nvelocity converged\n");

								if(u_on && (time_converge>=4 || t_current<5)){
										  u_converge=0;
										  cal_W();
								}

								time_converge=0;

								while(Q_on && !Q_converge && t_current<t_max){

										  t_current++;

										  evol_Q(Q1,Q2);
										  evol_Q(Q2,Q1);

										  time_converge++;

										  if(t_current%t_print==0){
													 printf("\nt=%d:\n",t_current);
													 monitor(-1,1);
													 if(movie) print_visual(0);
										  }
										  else {
													 monitor(-1,0);
										  }
								}

								if(Q_on && time_converge>=4){
										  Q_converge=0;
										  cal_stress(Q1);
										  cal_dtau(Q1);
								}

								if(Q_on) printf("\nQ-field converged\n");
					 }
		  }
		  else {			// simulataneous evolution
					 printf("Simultaneous evolution\n\n");
					 t_current = 0;
					 Q_converge = 0;
					 u_converge = 0;
					 while( (Q_converge==0 || u_converge==0) && t_current<t_max){

								t_current++;
								cal_sigma();
								evol_f(f1,f2);
								if(Q_on){
										  cal_W();
										  evol_Q(Q1,Q2);
										  evol_Q(Q2,Q1);
										  cal_stress(Q1);
										  cal_dtau(Q1);
								}
								cal_sigma();
								evol_f(f2,f1);
								if(Q_on){
										  cal_W();
										  evol_Q(Q1,Q2);
										  evol_Q(Q2,Q1);
										  cal_stress(Q1);
										  cal_dtau(Q1);
								}
								if(t_current%t_print==0){
										  printf("\nt=%d:\n",t_current);
										  monitor(0,1);
										  if(movie) print_visual(0);
								}
								else {
										  monitor(0,0);
								}

					 }
		  }

		  evol_Q(Q1,Q2);

		  write_restart(Q1,Rho,u);
		  print_visual(1-movie);

		  free_all();

		  printf("FINISHED\n");

		  return 0;

}
