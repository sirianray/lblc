#include "lb.h"
#include <math.h>

void init(real f[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3])
{
		  int i, j, k, ii, jj, ip, jp, kp, im, jm, km;
		  real sita, sita0, lambda1, lambda2, nlocal[3], slocal, slocali;

		  //zero-out field variables
		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {
								for (k=0; k<Nz; k++) {
										  for(ii=0;ii<3;ii++){

													 u[i][j][k][ii] = 0;

													 for(jj=0;jj<3;jj++){
																Q0[i][j][k][ii][jj]=0;
													 }
										  }
								}
					 }
		  }

		  if(newrun==1){
					 for (i=0; i<Nx; i++) {
								for (j=0; j<Ny; j++) {
										  for (k=0; k<Nz; k++) {

													 Rho[i][j][k]=rho;

													 u[i][j][k][0]=0;

													 if(init_linear) {
																u[i][j][k][1] = uy_bottom+(uy_top-uy_bottom)*((double)k+0.5)/((double)(Nz));
													 } else {
																u[i][j][k][1] = 0.0;
													 }

													 u[i][j][k][2]=0;

													 if(cav_flag[i][j][k]==-1){
																u[i][j][k][0] = 0.0;
																u[i][j][k][1] = 0.0;
																u[i][j][k][2] = 0.0;
													 }


													 k_eng += 0.5*(u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2])*Rho[i][j][k];

													 if(Q_on && cav_flag[i][j][k]>-1){

																if(rand_init){
																		  real randx, randy, randz;
																		  randx = randvec();
																		  randy = randvec();
																		  randz = randvec();

																		  real mag = randx*randx + randy*randy + randz*randz;
																		  mag = sqrt(mag);

																		  randx /= mag;
																		  randy /= mag;
																		  randz /= mag;

																		  real Srand = 0.5*S;

																		  Q0[i][j][k][0][0] = Srand*(randx*randx - third);
																		  Q0[i][j][k][0][1] = Srand*(randx*randy);
																		  Q0[i][j][k][0][2] = Srand*(randx*randz);
																		  Q0[i][j][k][1][1] = Srand*(randy*randy - third);
																		  Q0[i][j][k][1][2] = Srand*(randy*randz);

																		  Q0[i][j][k][1][0] =  Q0[i][j][k][0][1];
																		  Q0[i][j][k][2][0] =  Q0[i][j][k][0][2];
																		  Q0[i][j][k][2][1] =  Q0[i][j][k][1][2];
																		  Q0[i][j][k][2][2] =  -Q0[i][j][k][0][0]-Q0[i][j][k][1][1];

																}
																else {
																		  Q0[i][j][k][0][0] =  S*(ninit[0]*ninit[0] - third);
																		  Q0[i][j][k][0][1] =  S*(ninit[0]*ninit[1]);
																		  Q0[i][j][k][0][2] =  S*(ninit[0]*ninit[2]);
																		  Q0[i][j][k][1][1] =  S*(ninit[1]*ninit[1] - third);
																		  Q0[i][j][k][1][2] =  S*(ninit[1]*ninit[2]);

																		  Q0[i][j][k][1][0] =  Q0[i][j][k][0][1];
																		  Q0[i][j][k][2][0] =  Q0[i][j][k][0][2];
																		  Q0[i][j][k][2][1] =  Q0[i][j][k][1][2];
																		  Q0[i][j][k][2][2] =  -Q0[i][j][k][0][0]-Q0[i][j][k][1][1];
																}
													 }
										  }
								}
					 }

					 //For finite anchoring
					 for(i=0;i<Nx;i++){
								for(j=0;j<Ny;j++){

										  Qsurf[i][j][0][0] =  S*(nbot[0]*nbot[0] - third);
										  Qsurf[i][j][0][1] =  S*(nbot[0]*nbot[1]);
										  Qsurf[i][j][0][2] =  S*(nbot[0]*nbot[2]);
										  Qsurf[i][j][0][3] =  S*(nbot[1]*nbot[1] - third);
										  Qsurf[i][j][0][4] =  S*(nbot[1]*nbot[2]);

										  Qsurf[i][j][1][0] =  S*(ntop[0]*ntop[0] - third);
										  Qsurf[i][j][1][1] =  S*(ntop[0]*ntop[1]);
										  Qsurf[i][j][1][2] =  S*(ntop[0]*ntop[2]);
										  Qsurf[i][j][1][3] =  S*(ntop[1]*ntop[1] - third);
										  Qsurf[i][j][1][4] =  S*(ntop[1]*ntop[2]);
								}
					 }
		  }
		  else  read_restart(Q0,Rho,u);

		  //Preferred anchoring
		  Qo_bottom[0] =  S*(nbot[0]*nbot[0] - third);
		  Qo_bottom[1] =  S*(nbot[0]*nbot[1]);
		  Qo_bottom[2] =  S*(nbot[0]*nbot[2]);
		  Qo_bottom[3] =  S*(nbot[1]*nbot[1] - third);
		  Qo_bottom[4] =  S*(nbot[1]*nbot[2]);

		  Qo_top[0] =  S*(ntop[0]*ntop[0] - third);
		  Qo_top[1] =  S*(ntop[0]*ntop[1]);
		  Qo_top[2] =  S*(ntop[0]*ntop[2]);
		  Qo_top[3] =  S*(ntop[1]*ntop[1] - third);
		  Qo_top[4] =  S*(ntop[1]*ntop[2]);

		  //populate Q for sides of cavity
		  for(i=0;i<3;i++){
					 for(j=0;j<3;j++){
								Qtop[i][j] = S*(ncavtop[i]*ncavtop[j] - OneThirdDelta(i,j));
								Qsides[i][j] = S*(ncavsides[i]*ncavsides[j] - OneThirdDelta(i,j));
					 }
		  }

		  cal_dQ(Q0, H, Po, convQ);

		  if(newrun==0) {
					 cal_stress_q(Q0, H, sigma_q);
					 cal_detau(Q0);
		  } else {
					 // initialize all stresses to 0
					 for (i=0; i<Nx; i++) {
								for (j=0; j<Ny; j++) {
										  for (k=0; k<Nz; k++) {
													 for(ii=0;ii<3;ii++) {
																for(jj=0;jj<3;jj++) sigma_q[i][j][k][ii][jj] = 0;
													 }
													 detau[i][j][k][0] = 0;
													 detau[i][j][k][1] = 0;
													 detau[i][j][k][2] = 0;
										  }
								}
					 }
		  }

		  cal_sigma(Q0);
		  cal_fequ(f,Rho,u,sigma);
		  cal_W(u);

}

real randvec(){
		  double r;

		  r = 2*( (double)rand()/(double)RAND_MAX - 1);
		  return r;
}

