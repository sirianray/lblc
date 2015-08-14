#include "lb.h"
#include <math.h>

void init(double *f, double *Q0)
{
		  int i, j, k, ii, jj, ip, jp, kp, im, jm, km, index;
		  int points = Nx*Ny*Nz;

		  //For hybrid initial condition
		  double theta0, theta, thetalocal, lambda1, lambda2, thetalocali, nlocal[3]={0};
		  theta0 = ntop[0]*nbot[0]+ntop[1]*nbot[1]+ntop[2]*nbot[2];
		  theta0 = theta0/sqrt(ntop[0]*ntop[0]+ntop[1]*ntop[1]+ntop[2]*ntop[2]);        
		  theta0 = theta0/sqrt(nbot[0]*nbot[0]+nbot[1]*nbot[1]+nbot[2]*nbot[2]);
		  theta0 = acos(theta0);

		  if(newrun==1){

					 for(i=0;i<(points);i++) Rho[i] = rho;

					 for (k=0; k<Nz; k++) {
								for (j=0; j<Ny; j++) {
										  for (i=0; i<Nx; i++) {

													 index = i + Nx*j + Ny*Nx*k;

													 if(cav_on && cav_flag[index]==-1){
																//for(ii=0;ii<3;ii++) u[3*index+ii] = 0;
																u[3*index+0] = 0;
																u[3*index+1] = uy_top*((double)k+0.5)/((double)Nz);
																u[3*index+2] = 0;
													 }


													 k_eng += 0.5*( u[3*index+0]*+u[3*index+0] + u[3*index+1]*+u[3*index+1] + u[3*index+2]*+u[3*index+2])*Rho[index];

													 if(Q_on || (cav_on && cav_flag[index]>=0)){

																if(rand_init){

																		  double randx, randy, randz;
																		  randx = randvec();
																		  randy = randvec();
																		  randz = randvec();

																		  double mag = randx*randx + randy*randy + randz*randz;
																		  mag = sqrt(mag);

																		  randx /= mag;
																		  randy /= mag;
																		  randz /= mag;

																		  double Srand = 0.5*S;

																		  Q0[5*index+0] = Srand*(randx*randx - third);
																		  Q0[5*index+1] = Srand*(randx*randy);
																		  Q0[5*index+2] = Srand*(randx*randz);
																		  Q0[5*index+3] = Srand*(randy*randy - third);
																		  Q0[5*index+4] = Srand*(randy*randz);

																}

																else {
																		  if (theta0<1e-2 && theta0>-1e-2) {
																					 lambda1 = 0.5;
																					 lambda2 = 0.5;
																		  } else {
																					 theta = theta0 * ((double)k+0.5)/((double)Nz);
																					 lambda1=(cos(theta)-cos(theta0)*cos(theta0-theta));
																					 lambda2=(cos(theta0-theta)-cos(theta)*cos(theta0));
																		  }
																		  nlocal[0] = lambda1 * nbot[0] + lambda2 * ntop[0];
																		  nlocal[1] = lambda1 * nbot[1] + lambda2 * ntop[1];
																		  nlocal[2] = lambda1 * nbot[2] + lambda2 * ntop[2];
																		  thetalocal = nlocal[0]*nlocal[0] + nlocal[1]*nlocal[1] + nlocal[2]*nlocal[2];
																		  thetalocali= 1.0/thetalocal;

																		  Q0[5*index+0] =  S*(nlocal[0]*nlocal[0]*thetalocali - third);
																		  Q0[5*index+1] =  S*(nlocal[0]*nlocal[1]*thetalocali);
																		  Q0[5*index+2] =  S*(nlocal[0]*nlocal[2]*thetalocali);
																		  Q0[5*index+3] =  S*(nlocal[1]*nlocal[1]*thetalocali - third);
																		  Q0[5*index+4] =  S*(nlocal[1]*nlocal[2]*thetalocali);
																}

													 }
										  }
								}
					 }

					 //For finite anchoring
					 for(i=0;i<Nx;i++){
								for(j=0;j<Ny;j++){

										  index = i + j*Nx;

										  if(!rand_init){

													 Qsurfbot[5*index+0] =  S*(nbot[0]*nbot[0] - third);
													 Qsurfbot[5*index+1] =  S*(nbot[0]*nbot[1]);
													 Qsurfbot[5*index+2] =  S*(nbot[0]*nbot[2]);
													 Qsurfbot[5*index+3] =  S*(nbot[1]*nbot[1] - third);
													 Qsurfbot[5*index+4] =  S*(nbot[1]*nbot[2]);

													 Qsurftop[5*index+0] =  S*(ntop[0]*ntop[0] - third);
													 Qsurftop[5*index+1] =  S*(ntop[0]*ntop[1]);
													 Qsurftop[5*index+2] =  S*(ntop[0]*ntop[2]);
													 Qsurftop[5*index+3] =  S*(ntop[1]*ntop[1] - third);
													 Qsurftop[5*index+4] =  S*(ntop[1]*ntop[2]);
										  }

										  if(rand_init){

													 double randx, randy, randz;
													 randx = randvec();
													 randy = randvec();
													 randz = randvec();

													 double mag = randx*randx + randy*randy + randz*randz;
													 mag = sqrt(mag);

													 randx /= mag;
													 randy /= mag;
													 randz /= mag;

													 double Srand = 0.5*S;

													 Qsurfbot[5*index+0] = Srand*(randx*randx - third);
													 Qsurfbot[5*index+1] = Srand*(randx*randy);
													 Qsurfbot[5*index+2] = Srand*(randx*randz);
													 Qsurfbot[5*index+3] = Srand*(randy*randy - third);
													 Qsurfbot[5*index+4] = Srand*(randy*randz);

													 randx = randvec();
													 randy = randvec();
													 randz = randvec();

													 mag = sqrt(mag);

													 randx /= mag;
													 randy /= mag;
													 randz /= mag;

													 Qsurftop[5*index+0] = Srand*(randx*randx - third);
													 Qsurftop[5*index+1] = Srand*(randx*randy);
													 Qsurftop[5*index+2] = Srand*(randx*randz);
													 Qsurftop[5*index+3] = Srand*(randy*randy - third);
													 Qsurftop[5*index+4] = Srand*(randy*randz);

										  }
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
		  Qtop[0] = S*(ncavtop[0]*ncavtop[0] - third);
		  Qtop[1] = S*(ncavtop[0]*ncavtop[1]);
		  Qtop[2] = S*(ncavtop[0]*ncavtop[2]);
		  Qtop[3] = S*(ncavtop[1]*ncavtop[1] - third);
		  Qtop[4] = S*(ncavtop[1]*ncavtop[2]);

		  Qsides[0] = S*(ncavsides[0]*ncavsides[0] - third);
		  Qsides[1] = S*(ncavsides[0]*ncavsides[1]);
		  Qsides[2] = S*(ncavsides[0]*ncavsides[2]);
		  Qsides[3] = S*(ncavsides[1]*ncavsides[1] - third);
		  Qsides[4] = S*(ncavsides[1]*ncavsides[2]);

		  cal_dQ(Q0);

		 // if(!newrun){
		 //  		 cal_stress(Q0);
		 //  		 cal_dtau(Q0);
		 // }

		  cal_sigma();
		  cal_fequ(f,u);
		  cal_W();

}

double randvec(){
		  double r;

		  r = 2*( (double)rand()/(double)RAND_MAX - 1);
		  return r;
}

