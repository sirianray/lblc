#include "lb.h" 
void evol_f(double *fin, double *fout)
{
		  int i,j,k,ii,jj, index, next;
		  int ip, jp, kp;
		  double Cfin, Cfout, CG1, CG2;
		  double temp;

		  cal_feqf_new(fin);
		  cal_p(fin);

		  streaming(fin,fout);

		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Nx + k*Nx*Ny;

										  if( (cav_on && cav_flag[index]==0) || !cav_on) { //Evolve bulk points first

													 for(ii=0;ii<15;ii++){

																Cf[15*index+ii] = - itau_f * (fin[15*index+ii] - f_eq[15*index+ii]) + p[15*index+ii];

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

																next = ip + jp*Nx + kp*Nx*Ny; 

																fout[15*next+jj] += dt * Cf[15*index+ii];

													 }
										  }

										  if(cav_on && cav_flag[index]==1){ //Evolve cavity points

													 for(ii=0;ii<15;ii++){

																Cf[15*index+ii] = - itau_f * (fin[15*index+ii] - f_eq[15*index+ii]) + p[15*index+ii];

																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if(ip < 0)  ip += Nx;
																if(ip>=Nx) ip -= Nx;

																if(jp < 0)  jp += Ny;
																if(jp>=Ny) jp -= Ny;

																jj=ii;

																next = ip + jp*Nx + kp*Nx*Ny;
																
																if(cav_flag[next]==-1 || kp<0){

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																}


																next = ip + jp*Nx + kp*Nx*Ny;
																fout[15*next+jj] += dt * Cf[15*index+ii];
													 }
										  }
								}
					 }
		  }


		  cal_rho(fout);
		  cal_u(fout);
		  cal_sigma();
		  cal_feqf_new(fout);
		  cal_p(fout);

		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Nx + k*Nx*Ny;

										  if( (cav_on && cav_flag[index]==0) || !cav_on){

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

																next = ip + jp*Nx + kp*Nx*Ny;

																Cfin = Cf[15*index+ii];

																Cfout = - itau_f * (fout[15*next+jj] - f_eq[15*next+jj]) + p[15*next+jj];
																fin[15*index+ii] = fin[15*index+ii] + dt * 0.5 * (Cfin + Cfout);
													 }
										  }

										  if(cav_flag[index]==1 && cav_on){

													 for(ii=0;ii<15;ii++) {

																ip = i + e[ii][0];
																jp = j + e[ii][1];
																kp = k + e[ii][2];

																if(ip <0)  ip += Nx;
																if(ip>=Nx) ip -= Nx;

																if(jp <0)  jp += Ny;
																if(jp>=Ny) jp -= Ny;

																jj=ii;

																next = ip + jp*Nx + kp*Nx*Ny;

																if(kp<0 || cav_flag[next]==-1){

																		  ip = i;
																		  jp = j;
																		  kp = k;

																		  jj = bounce(ii);
																}

																next = ip + jp*Nx + kp*Nx*Ny;

																Cfin = Cf[15*index+ii];

																Cfout = - itau_f * (fout[15*next+jj] - f_eq[15*next+jj]) + p[15*next+jj];
																
																fin[15*index+ii] = fin[15*index+ii] + dt * 0.5 * (Cfin + Cfout);
													 }
										  }

								}
					 }
		  }

		  streaming(fin,fout);

		  cal_rho(fout);
		  cal_u(fout);
}

void streaming(double *fin, double *fout)
{
		  int  i, j, k, ii, iip, jj, kk, in, jn, kn;
		  int index, next;
		  double Bx2;
		  int  ipar;


		  for (ii=0; ii<15; ii++) {

					 if(ii<7) {
								Bx2 = two3rd;
					 }else {
								Bx2 = one12th;
					 }

					 iip = bounce(ii);

					 for (k=0; k<Nz; k++) {
								for (j=0; j<Ny; j++) {
										  for (i=0; i<Nx; i++) {

													 index = i + j*Nx + k*Nx*Ny;

													 in = i + e[ii][0];
													 jn = j + e[ii][1];
													 kn = k + e[ii][2];

													 if(in==-1) in=Nx-1;
													 if(in==Nx) in=0;

													 if(jn==-1) jn=Ny-1;
													 if(jn==Ny) jn=0;

													 next = in + jn*Nx + kn*Nx*Ny;

													 if(!cav_on || cav_flag[index]==0){

																if (kn>=0 && kn<Nz) { //if adjacent point is within bulk

																		  fout[15*next+ii]=fin[15*index+ii];

																} else { //if at top or bottom wall

																		  fout[15*index+iip]=fin[15*index+ii];

																		  if(k==0)     fout[15*index+iip] += -Bx2*Rho[index]*(uy_bottom*(double)e[ii][1]);
																		  if(k==Nz-1) fout[15*index+iip] += -Bx2*Rho[index]*(uy_top*(double)e[ii][1]);
																}
													 }

													 if(cav_on && cav_flag[index]==1) {

																if(kn>=0){
																		  if(cav_flag[next]==-1) fout[15*index+iip]=fin[15*index+ii]; //bounce back
																		  if(cav_flag[next]>= 0) fout[15*next+ii]=fin[15*index+ii]; //adjacent point is bulk
																}

																if(kn<0) fout[15*index+iip]=fin[15*index+ii]; //bottom corner points bounce back 

													 }
										  }
								}
					 }
		  }

}
