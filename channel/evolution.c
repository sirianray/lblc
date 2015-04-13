#include "lb.h" 
void evol_f(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15])
{
		  int i,j,k,ii,jj;
		  int ip, jp, kp;
		  real Cf1, Cf2, CG1, CG2;
		  real temp;

		  cal_feqf_new(f_eq, f1, sigma);
		  cal_p(p, f1, f_eq);

		  streaming(f1,f2);

		  for(i=0;i<Nx;i++){
					 for(j=0;j<Ny;j++){
								for(k=0;k<Nz;k++){

										  if( (cav_on && cav_flag[i][j][k]==0) || !cav_on) { //Evolve bulk points first

													 for(ii=0;ii<15;ii++){

																Cf[i][j][k][ii] = - itau_f * (f1[i][j][k][ii] - f_eq[i][j][k][ii]) + p[i][j][k][ii];

																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if (kp>=0 && kp<Nz) { // if adjacent site is within bulk

																		  if (ip >= Nx) {
																					 ip -= Nx;
																		  } else if (ip<0) {
																					 ip += Nx;
																		  }

																		  if (jp >= Ny) {
																					 jp -= Ny;
																		  } else if (jp<0) {
																					 jp += Ny;
																		  }

																		  jj = ii;
																}
																else { // if adjacent site passes boundary 

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																}


																f2[ip][jp][kp][jj] += dt * Cf[i][j][k][ii];
																//					f1[i][j][k][ii] = f1[i][j][k][ii] + dt * Cf[i][j][k][ii];

													 }
										  }

										  if(cav_on && cav_flag[i][j][k]==1){ //Evolve cavity points

													 for(ii=0;ii<15;ii++){

																Cf[i][j][k][ii] = - itau_f * (f1[i][j][k][ii] - f_eq[i][j][k][ii]) + p[i][j][k][ii];

																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if(ip <0)  ip += Nx;
																if(ip>=Nx) ip -= Nx;

																if(jp <0)  jp += Ny;
																if(jp>=Ny) jp -= Ny;

																jj=ii;

																if(kp>=0){
																		  if(cav_flag[ip][jp][kp]==-1){

																					 ip = i;
																					 jp = j;
																					 kp = k;

																					 jj = bounce(ii);
																		  }
																}

																if(kp<0){

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																		 
																}

																f2[ip][jp][kp][jj] += dt * Cf[i][j][k][ii];
																//					f1[i][j][k][ii] = f1[i][j][k][ii] + dt * Cf[i][j][k][ii];
													 }
										  }
								}
					 }
		  }


		  cal_rho(Rho,f2);
		  cal_u(u,f2);
		  cal_sigma(Q1);
		  cal_feqf_new(f_eq,f2,sigma);
		  cal_p(p,f2,f_eq);

		  for(k=0;k<Nz;k++){
					 for(i=0;i<Nx;i++){
								for(j=0;j<Ny;j++){

										  if( (cav_on && cav_flag[i][j][k]==0) || !cav_on){

													 for(ii=0;ii<15;ii++) {
																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if (kp>=0 && kp <Nz) {
																		  if (ip >= Nx) {
																					 ip -= Nx;
																		  } else if (ip<0) {
																					 ip += Nx;
																		  }
																		  if (jp >= Ny) {
																					 jp -= Ny;
																		  } else if (jp<0) {
																					 jp += Ny;
																		  }

																		  jj = ii;

																} else {

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																}

																Cf1 = Cf[i][j][k][ii];

																Cf2 = - itau_f * (f2[ip][jp][kp][jj] - f_eq[ip][jp][kp][jj]) + p[ip][jp][kp][jj];
																//					Cf2 = - (f2[i][j][k][ii] - f_eq[i][j][k][ii]);
																//					Cf2 = Cf1;
																
																f1[i][j][k][ii] = f1[i][j][k][ii] + dt * 0.5 * (Cf1 + Cf2);
													 }
										  }

										  if(cav_flag[i][j][k]==1 && cav_on){

													 for(ii=0;ii<15;ii++) {

																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if(ip <0)  ip += Nx;
																if(ip>=Nx) ip -= Nx;

																if(jp <0)  jp += Ny;
																if(jp>=Ny) jp -= Ny;

																jj=ii;

																if(kp>=0){
																		  if(cav_flag[ip][jp][kp]==-1){

																					 ip = i;
																					 jp = j;
																					 kp = k;

																					 jj = bounce(ii);
																		  }
																}

																if(kp<0){

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																}

																Cf1 = Cf[i][j][k][ii];

																Cf2 = - itau_f * (f2[ip][jp][kp][jj] - f_eq[ip][jp][kp][jj]) + p[ip][jp][kp][jj];
																//					Cf2 = - (f2[i][j][k][ii] - f_eq[i][j][k][ii]);
																//					Cf2 = Cf1;
																
																f1[i][j][k][ii] = f1[i][j][k][ii] + dt * 0.5 * (Cf1 + Cf2);
													 }
										  }

								}
					 }
		  }

		  streaming(f1,f2);

		  cal_rho(Rho,f2);
		  cal_u(u, f2);
		  //	cal_W(u);
		  //	cal_sigma();
}

void streaming(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15])
{
		  int  i, j, k, ii, iip, jj, kk, in, jn, kn;
		  real Bx2;
		  real dx, dy, dz, dr2i, rrad;
		  int  ipar;


		  for (ii=0; ii<15; ii++) {

					 if(ii<7) {
								Bx2 = two3rd;
					 } else {
								Bx2 = one12th;
					 }

					 iip = bounce(ii);

					 for (k=0; k<Nz; k++) {
								for (i=0; i<Nx; i++) {
										  for (j=0; j<Ny; j++) {

													 in = i + e[ii][0];
													 jn = j + e[ii][1];
													 kn = k + e[ii][2];

													 if(in==-1) in=Nx-1;
													 if(in==Nx) in=0;

													 if(jn==-1) jn=Ny-1;
													 if(jn==Ny) jn=0;

													 if(!cav_on || cav_flag[i][j][k]==0){

																if (kn>=0 && kn<Nz) { //if adjacent point is within bulk

																		  f2[in][jn][kn][ii]=f1[i][j][k][ii];

																} else { //if at top or bottom wall

																		  if (pbc_z==1){

																					 if(kn==-1) kn=Nz-1;
																					 if(kn==Nz) kn=0;

																					 f2[in][jn][kn][ii]=f1[i][j][k][ii];
																		  } 

																		  else { 

																					 f2[i][j][k][iip]=f1[i][j][k][ii];

																					 if (k==0) {
																								f2[i][j][k][iip] += -Bx2*Rho[i][j][k]*(uy_bottom*(double)e[ii][1]);
																					 } 
																					 else if(k==Nz-1) {
																								f2[i][j][k][iip] += -Bx2*Rho[i][j][k]*(uy_top*(double)e[ii][1]);
																					 }
																		  } 
																}
													 }

													 if(cav_on && cav_flag[i][j][k]==1) {

																if(kn>=0){
																		  if(cav_flag[in][jn][kn]==-1) f2[i][j][k][iip]=f1[i][j][k][ii]; //bounce back
																		  if(cav_flag[in][jn][kn]>= 0) f2[in][jn][kn][ii]=f1[i][j][k][ii]; //adjacent point is bulk
																}

																if(kn<0) f2[i][j][k][iip]=f1[i][j][k][ii]; //bottom corner points bounce back 

													 }
										  }
								}
					 }
		  }

}
