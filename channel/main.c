#include "lb.h"

int e[15][3], t_print, t_max, movie;
int pbc_z, wall_anchoring, Q_on, u_converge, debug_mode, simul_evol, u_on, Q_converge;
int t_current=0, newrun, t_v_on, init_linear;
real dt=1.0, tau_f, tau_G, kappa, T=third, xi, Gamma_rot, rho, yforce, U, A_ldg;
real uy_top, uy_bottom;
real itau_f, itau_G;
real H[Nx][Ny][Nz][3][3], convQ[Nx][Ny][Nz][5];
real sigma[Nx][Ny][Nz][3][3], u[Nx][Ny][Nz][3], Po[Nx][Ny][Nz];
real W[Nx][Ny][Nz][3][3];
real S, Rho[Nx][Ny][Nz];
real f_eq[Nx][Ny][Nz][15], Cf[Nx][Ny][Nz][15], p[Nx][Ny][Nz][15];
real Qo_top[5], Qo_bottom[5];
real f1[Nx][Ny][Nz][15], f2[Nx][Ny][Nz][15];
real Q1[Nx][Ny][Nz][3][3], Q2[Nx][Ny][Nz][3][3];
real k_eng=1, sigma_q[Nx][Ny][Nz][3][3], tau[Nx][Ny][Nz][3][3], detau[Nx][Ny][Nz][3], lap_Q[Nx][Ny][Nz];
real mcheck[Nx][Ny][Nz][3];

//Finite anchoring
real Qsurf[Nx][Ny][2][5]; // 0 = bottom, 1 = top
real dQdxnu[Nx][Ny][2][5];
real W_bot, W_top, nbot[3], ntop[3], ninit[3];
int rand_seed, rand_init;

//Cavity-Channel parameters
int cav_height, cav_width, cav_flag[Nx][Ny][Nz], cav_on;
real Qtop[3][3], Qsides[3][3];
real ncavtop[3], ncavsides[3];

int main()
{
		  int t,i,j,k,ii,jj,kk;
		  FILE *out1, *out2, *out3, *out4;
		  int time_converge;

		  Q_converge = 0;
		  u_converge = 0;

		  read_param();
		  lattice_vec();

		  //discover cavity points
		  cav_surf();

		  init(f1, Q1);

		  monitor(0,1);
		  print_visual(1);

		  //equilibrate Q-field 
		  if(rand_init && Q_on) {
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
					 cal_stress_q(Q1,H,sigma_q);
					 cal_detau(Q1);
		  }

		  if(u_on) printf("\nk_eng=%f\n",k_eng);

		  //	evolution
		  if (simul_evol==0) {	//	alternating evolution
					 while( ( (Q_converge==0 && Q_on) || (u_converge==0 && u_on) ) && t_current<t_max){

								time_converge = 0;

								while(u_on==1 && u_converge==0 && t_current<t_max){

										  t_current++;

										  cal_sigma(Q1);
										  evol_f(f1,f2);
										  cal_sigma(Q1);
										  evol_f(f2,f1);

										  time_converge++;

										  if(t_current%t_print==0){

													 printf("\nt=%d:\n",t_current);
													 monitor(1,1);
													 if(movie) print_visual(0);
										  }
										  else {
													 monitor(1,0);
										  }
								}

								if(u_on && (time_converge>=4 || t_current<5)){
										  u_converge=0;
										  cal_W(u);
								}

								time_converge=0;

								while(Q_on==1 && Q_converge==0 && t_current<t_max){

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
										  cal_stress_q(Q1,H,sigma_q);
										  cal_detau(Q1);
								}
					 }
		  }
		  else {			// simulataneous evolution
		  t_current = 0;
		  Q_converge = 0;
		  u_converge = 0;
					 while( (Q_converge==0 || u_converge==0) && t_current<t_max){
								t_current++;
								cal_sigma(Q1);
								monitor(0,0);
								printf("1\n");
								evol_f(f1,f2);
								monitor(0,0);
								printf("2\n");
								cal_W(u);
								print_visual(1);
								evol_Q(Q1,Q2);
								monitor(0,0);
								printf("3\n");
								cal_stress_q(Q2,H,sigma_q);
								cal_detau(Q2);
								cal_sigma(Q2);
								monitor(0,0);
								printf("4\n");

								evol_f(f2,f1);
								//print_u2(out1, out2, out3, out4);
								cal_W(u);
								evol_Q(Q2,Q1);
								cal_stress_q(Q1,H,sigma_q);
								cal_detau(Q1);
								cal_sigma(Q1);
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


		  if(debug_mode!=0){
					 printf("H\n");
					 for(i=(Nx+1)/2;i<(Nx+3)/2;i++){
								for(j=(Ny+1)/2;j<(Ny+3)/2;j++){
										  for(k=0;k<Nz;k++){
													 for(ii=0;ii<3;ii++){
																for(jj=0;jj<3;jj++)printf(" %13.8f",H[i][j][k][ii][jj]);
													 }
													 printf("\n");
										  }
								}
					 }
		  }


		  if(debug_mode!=0) check(Q1);

		  write_restart(Q1,Rho,u);
		  print_visual(1-movie);
					 
		  return 0;

}
