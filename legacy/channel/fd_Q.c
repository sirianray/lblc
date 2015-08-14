#include "lb.h"

void evol_Q(real Q1[Nx][Ny][Nz][3][3], real Q2[Nx][Ny][Nz][3][3])
{
		  int i, j, k, ii, jj, kk;
		  real D[3][3], Omega[3][3], S[3][3], trQW, temp;

		  cal_dQ(Q1, H, Po, convQ);

		  //Bulk evolution
		  for(i=0;i<Nx;i++){
					 for(j=0;j<Ny;j++){
								for(k=0;k<Nz;k++){

										  if((cav_flag[i][j][k]>=0 && cav_on) || !cav_on) {

													 trQW=0.0;
													 for(ii=0;ii<3;ii++){
																for(jj=0;jj<3;jj++){
																		  D[ii][jj] = 0.5*(W[i][j][k][ii][jj] + W[i][j][k][jj][ii]);
																		  Omega[ii][jj] = 0.5*(W[i][j][k][ii][jj] - W[i][j][k][jj][ii]);
																		  trQW = trQW + Q1[i][j][k][ii][jj]*W[i][j][k][jj][ii];
																}
													 }

													 for(ii=0;ii<3;ii++){
																for(jj=0;jj<3;jj++){
																		  temp = 0;
																		  for(kk=0;kk<3;kk++){
																					 temp += (xi*D[ii][kk] + Omega[ii][kk])*(Q1[i][j][k][kk][jj] + OneThirdDelta(kk,jj)) + (Q1[i][j][k][ii][kk] + OneThirdDelta(ii,kk))*(xi*D[kk][jj] - Omega[kk][jj]);
																		  }
																		  temp += - 2.0*xi*(Q1[i][j][k][ii][jj] + OneThirdDelta(ii,jj) )*trQW;
																		  S[ii][jj] = temp;
																}
													 }

													 Q2[i][j][k][0][0] = Q1[i][j][k][0][0] + dt * (Gamma_rot * H[i][j][k][0][0] + S[0][0] - convQ[i][j][k][0]);
													 Q2[i][j][k][0][1] = Q1[i][j][k][0][1] + dt * (Gamma_rot * H[i][j][k][0][1] + S[0][1] - convQ[i][j][k][1]);
													 Q2[i][j][k][0][2] = Q1[i][j][k][0][2] + dt * (Gamma_rot * H[i][j][k][0][2] + S[0][2] - convQ[i][j][k][2]);
													 Q2[i][j][k][1][1] = Q1[i][j][k][1][1] + dt * (Gamma_rot * H[i][j][k][1][1] + S[1][1] - convQ[i][j][k][3]);
													 Q2[i][j][k][1][2] = Q1[i][j][k][1][2] + dt * (Gamma_rot * H[i][j][k][1][2] + S[1][2] - convQ[i][j][k][4]);

													 Q2[i][j][k][1][0] = Q2[i][j][k][0][1];
													 Q2[i][j][k][2][0] = Q2[i][j][k][0][2];
													 Q2[i][j][k][2][1] = Q2[i][j][k][1][2];
													 Q2[i][j][k][2][2] =-Q2[i][j][k][0][0] - Q2[i][j][k][1][1];

										  }
								}
					 }
		  }

		  //Evoluation for top and bottom surface
		  if(!cav_on){
					 int n=0;
					 if(W_bot>=0 || W_top>=0){
								surface_derivative(Q1);

								for(i=0;i<Nx;i++){
										  for(j=0;j<Ny;j++){
													 for(n=0;n<5;n++){

																if(W_bot>=0) Qsurf[i][j][0][n] = Qsurf[i][j][0][n] - dt*Gamma_rot * (-1 * kappa * dQdxnu[i][j][0][n] + W_bot*( Qsurf[i][j][0][n] - Qo_bottom[n]) );
																if(W_top>=0) Qsurf[i][j][1][n] = Qsurf[i][j][1][n] - dt*Gamma_rot * (kappa * dQdxnu[i][j][1][n] + W_top*( Qsurf[i][j][1][n] - Qo_top[n]));

													 }
										  }
								}
					 }
		  }
}

void surface_derivative(real Q0[Nx][Ny][Nz][3][3]){

		  //bottom surface
		  int i, j, k, n;

		  real Qvec1[5], Qvec0[5];
		  real QvecNz1[5], QvecNz2[5];

		  for(i=0;i<5;i++){
					 Qvec1[i] = 0;
					 Qvec0[i] = 0;
					 QvecNz1[i] = 0;
					 QvecNz2[i] = 0;
		  }

		  //		  real dx1 = 0.5, dx2 = 1.0, idx = dx1*dx2*(dx1+dx2);
		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {

								Qvec0[0] = Q0[i][j][0][0][0];
								Qvec0[1] = Q0[i][j][0][0][1];
								Qvec0[2] = Q0[i][j][0][0][2];
								Qvec0[3] = Q0[i][j][0][1][1];
								Qvec0[4] = Q0[i][j][0][1][2];

								Qvec1[0] = Q0[i][j][1][0][0];
								Qvec1[1] = Q0[i][j][1][0][1];
								Qvec1[2] = Q0[i][j][1][0][2];
								Qvec1[3] = Q0[i][j][1][1][1];
								Qvec1[4] = Q0[i][j][1][1][2];

								for (n=0; n<5; n++) {
										  //Bottom Surface
										  //										  dQdxnu[i][j][0][n] = idx*(-dx1*dx1*Qvec1[n] + (dx1+dx2)*(dx1+dx2)*Qvec0[n] - (dx2*dx2+2*dx1*dx2)*Qsurf[i][j][0][n]);
										  dQdxnu[i][j][0][n] = -eight3rd*Qsurf[i][j][0][n] + 3.0*Qvec0[n] - third*Qvec1[n];

								}
					 }
		  }


		  //		  dx1 = 1.0, dx2 = 0.5, idx = dx1*dx2*(dx1+dx2);
		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {

								QvecNz1[0] = Q0[i][j][Nz-1][0][0];
								QvecNz1[1] = Q0[i][j][Nz-1][0][1];
								QvecNz1[2] = Q0[i][j][Nz-1][0][2];
								QvecNz1[3] = Q0[i][j][Nz-1][1][1];
								QvecNz1[4] = Q0[i][j][Nz-1][1][2];

								QvecNz2[0] = Q0[i][j][Nz-2][0][0];
								QvecNz2[1] = Q0[i][j][Nz-2][0][1];
								QvecNz2[2] = Q0[i][j][Nz-2][0][2];
								QvecNz2[3] = Q0[i][j][Nz-2][1][1];
								QvecNz2[4] = Q0[i][j][Nz-2][1][2];

								for (n=0; n<5; n++) {
										  //Top Surface
										  //										  dQdxnu[i][j][1][n] = idx*( (dx1*dx1+2*dx1*dx2)*Qsurf[i][j][1][n] - (dx1+dx2)*(dx1+dx2)*QvecNz1[n] + dx2*dx2 * QvecNz2[n]);
										  dQdxnu[i][j][1][n] = eight3rd*Qsurf[i][j][1][n] - 3.0*QvecNz1[n] + third*QvecNz2[n];
								}
					 }
		  }

}

void cal_dQ(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real Po[Nx][Ny][Nz], real convQ[Nx][Ny][Nz][5])
{
		  int i, j, k, ii=0, jj=0, kk=0, nn=0;
		  int ip, im, jp, jm, kp, km;

		  real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
		  t10=t11=t12=t20=t21=t22=t30=t31=t32=t40=t41=t42=t50=t51=t52=0;

		  real temp=0, Qu[5]={0}, Qd[5]={0};

		  real d2x0, d2x1, d2x2, d2x3, d2x4, d2y0, d2y1, d2y2, d2y3, d2y4, d2z0, d2z1, d2z2, d2z3, d2z4;
		  d2x0=d2x1=d2x2=d2x3=d2x4=d2y0=d2y1=d2y2=d2y3=d2y4=d2z0=d2z1=d2z2=d2z3=d2z4=0;

		  real dx, dy, dz, dr2i;
		  dx=dy=dz=dr2i=0;

		  real a0, a1, a2, a3, a4;
		  a0=a1=a2=a3=a4=0;

		  real trQ2=0, trH3rd=0, d2Q[3][3]={0}, dQ2=0;
		  real ux=0, uy=0, uz=0;

		  int halfNy = 0.5*Ny;
		  int h = cav_height;
		  real Qxu[5]={0}, Qxd[5]={0}, Qyu[5]={0}, Qyd[5]={0}, Qzu[5]={0}, Qzd[5]={0};

		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {
								for (k=0; k<Nz; k++) {

										  //Q at current point
										  a0 = Q0[i][j][k][0][0];
										  a1 = Q0[i][j][k][0][1];
										  a2 = Q0[i][j][k][0][2];
										  a3 = Q0[i][j][k][1][1];
										  a4 = Q0[i][j][k][1][2];

										  im = i-1;
										  jm = j-1;
										  km = k-1;

										  ip = i+1;
										  jp = j+1;
										  kp = k+1;

										  //PBC
										  if(i==0) im = Nx-1;
										  if(j==0) jm = Ny-1;
										  if(k==0) km = Nz-1;

										  if(i==Nx-1) ip = 0;
										  if(j==Ny-1) jp = 0;
										  if(k==Nz-1) kp = 0;

										  //Bulk points away from cavity
										  if( (cav_flag[i][j][k]==0 && cav_on==1) || cav_on==0){

													 //Bottom Wall
													 if(k==0 && (wall_anchoring && W_bot<0) ){
																Qd[0] = Qo_bottom[0];
																Qd[1] = Qo_bottom[1];
																Qd[2] = Qo_bottom[2];
																Qd[3] = Qo_bottom[3];
																Qd[4] = Qo_bottom[4];
													 } 
													 else if(k==0 && (wall_anchoring && W_bot>=0) ){
																Qd[0] = Qsurf[i][j][0][0];
																Qd[1] = Qsurf[i][j][0][1];
																Qd[2] = Qsurf[i][j][0][2];
																Qd[3] = Qsurf[i][j][0][3];
																Qd[4] = Qsurf[i][j][0][4];
													 } 
													 else {
																Qd[0] = Q0[i][j][km][0][0];
																Qd[1] = Q0[i][j][km][0][1];
																Qd[2] = Q0[i][j][km][0][2];
																Qd[3] = Q0[i][j][km][1][1];
																Qd[4] = Q0[i][j][km][1][2];
													 }

													 //Top Wall
													 if(k==Nz-1 && (wall_anchoring && W_top<0) ) {
																Qu[0] = Qo_top[0];
																Qu[1] = Qo_top[1];
																Qu[2] = Qo_top[2];
																Qu[3] = Qo_top[3];
																Qu[4] = Qo_top[4];

													 } 
													 else if(k==Nz-1 && (wall_anchoring !=0 && W_top>=0) ) {
																Qu[0] = Qsurf[i][j][1][0];
																Qu[1] = Qsurf[i][j][1][1];
																Qu[2] = Qsurf[i][j][1][2];
																Qu[3] = Qsurf[i][j][1][3];
																Qu[4] = Qsurf[i][j][1][4];
													 } 
													 else {
																Qu[0] = Q0[i][j][kp][0][0];
																Qu[1] = Q0[i][j][kp][0][1];
																Qu[2] = Q0[i][j][kp][0][2];
																Qu[3] = Q0[i][j][kp][1][1];
																Qu[4] = Q0[i][j][kp][1][2];
													 }

													 //x-first and second derivatives
													 d2x0 = Q0[ip][j][k][0][0] + Q0[im][j][k][0][0] - 2.0*a0;
													 d2x1 = Q0[ip][j][k][0][1] + Q0[im][j][k][0][1] - 2.0*a1;
													 d2x2 = Q0[ip][j][k][0][2] + Q0[im][j][k][0][2] - 2.0*a2;
													 d2x3 = Q0[ip][j][k][1][1] + Q0[im][j][k][1][1] - 2.0*a3;
													 d2x4 = Q0[ip][j][k][1][2] + Q0[im][j][k][1][2] - 2.0*a4;
													 t10  = 0.5 * (Q0[ip][j][k][0][0]-Q0[im][j][k][0][0]);
													 t20  = 0.5 * (Q0[ip][j][k][0][1]-Q0[im][j][k][0][1]);
													 t30  = 0.5 * (Q0[ip][j][k][0][2]-Q0[im][j][k][0][2]);
													 t40  = 0.5 * (Q0[ip][j][k][1][1]-Q0[im][j][k][1][1]);
													 t50  = 0.5 * (Q0[ip][j][k][1][2]-Q0[im][j][k][1][2]);

													 //y-first and second derivatives
													 d2y0 = Q0[i][jm][k][0][0]+Q0[i][jp][k][0][0]-2.0*a0;
													 d2y1 = Q0[i][jm][k][0][1]+Q0[i][jp][k][0][1]-2.0*a1;
													 d2y2 = Q0[i][jm][k][0][2]+Q0[i][jp][k][0][2]-2.0*a2;
													 d2y3 = Q0[i][jm][k][1][1]+Q0[i][jp][k][1][1]-2.0*a3;
													 d2y4 = Q0[i][jm][k][1][2]+Q0[i][jp][k][1][2]-2.0*a4;
													 t11  = 0.5 * (Q0[i][jp][k][0][0]-Q0[i][jm][k][0][0]);
													 t21  = 0.5 * (Q0[i][jp][k][0][1]-Q0[i][jm][k][0][1]);
													 t31  = 0.5 * (Q0[i][jp][k][0][2]-Q0[i][jm][k][0][2]);
													 t41  = 0.5 * (Q0[i][jp][k][1][1]-Q0[i][jm][k][1][1]);
													 t51  = 0.5 * (Q0[i][jp][k][1][2]-Q0[i][jm][k][1][2]);

													 if(k==0 && wall_anchoring !=0 ){
																d2z0 = eight3rd*Qd[0] + four3rd*Qu[0] - 4.0 * a0;
																d2z1 = eight3rd*Qd[1] + four3rd*Qu[1] - 4.0 * a1;
																d2z2 = eight3rd*Qd[2] + four3rd*Qu[2] - 4.0 * a2;
																d2z3 = eight3rd*Qd[3] + four3rd*Qu[3] - 4.0 * a3;
																d2z4 = eight3rd*Qd[4] + four3rd*Qu[4] - 4.0 * a4;
																t12  =-four3rd*Qd[0] + third*Qu[0] +       a0;
																t22  =-four3rd*Qd[1] + third*Qu[1] +       a1;
																t32  =-four3rd*Qd[2] + third*Qu[2] +       a2;
																t42  =-four3rd*Qd[3] + third*Qu[3] +       a3;
																t52  =-four3rd*Qd[4] + third*Qu[4] +       a4;
													 } 

													 else if(k==Nz-1 && wall_anchoring!=0){
																d2z0 = eight3rd*Qu[0] + four3rd*Qd[0] - 4.0 * a0;
																d2z1 = eight3rd*Qu[1] + four3rd*Qd[1] - 4.0 * a1;
																d2z2 = eight3rd*Qu[2] + four3rd*Qd[2] - 4.0 * a2;
																d2z3 = eight3rd*Qu[3] + four3rd*Qd[3] - 4.0 * a3;
																d2z4 = eight3rd*Qu[4] + four3rd*Qd[4] - 4.0 * a4;
																t12  = four3rd*Qu[0] - third*Qd[0] -       a0;
																t22  = four3rd*Qu[1] - third*Qd[1] -       a1;
																t32  = four3rd*Qu[2] - third*Qd[2] -       a2;
																t42  = four3rd*Qu[3] - third*Qd[3] -       a3;
																t52  = four3rd*Qu[4] - third*Qd[4] -       a4;
													 }
													 else {
																d2z0 = Qd[0] + Qu[0] - 2.0*a0;
																d2z1 = Qd[1] + Qu[1] - 2.0*a1;
																d2z2 = Qd[2] + Qu[2] - 2.0*a2;
																d2z3 = Qd[3] + Qu[3] - 2.0*a3;
																d2z4 = Qd[4] + Qu[4] - 2.0*a4;
																t12  = 0.5 * ( Qu[0] - Qd[0] );
																t22  = 0.5 * ( Qu[1] - Qd[1] );
																t32  = 0.5 * ( Qu[2] - Qd[2] );
																t42  = 0.5 * ( Qu[3] - Qd[3] );
																t52  = 0.5 * ( Qu[4] - Qd[4] );
													 }

										  }

										  //Derivatives for Cavity
										  if(cav_on){
													 if(cav_flag[i][j][k]==1){

																//if on left side of cavity
																if(j<halfNy && k<h){

																		  Qxu[0] = Q0[ip][j][k][0][0];
																		  Qxu[1] = Q0[ip][j][k][0][1];
																		  Qxu[2] = Q0[ip][j][k][0][2];
																		  Qxu[3] = Q0[ip][j][k][1][1];
																		  Qxu[4] = Q0[ip][j][k][1][2];

																		  Qxd[0] = Q0[im][j][k][0][0];
																		  Qxd[1] = Q0[im][j][k][0][1];
																		  Qxd[2] = Q0[im][j][k][0][2];
																		  Qxd[3] = Q0[im][j][k][1][1];
																		  Qxd[4] = Q0[im][j][k][1][2];

																		  Qyu[0] = Qsides[0][0];
																		  Qyu[1] = Qsides[0][1];
																		  Qyu[2] = Qsides[0][2];
																		  Qyu[3] = Qsides[1][1];
																		  Qyu[4] = Qsides[1][2];

																		  Qyd[0] = Q0[i][jm][k][0][0];
																		  Qyd[1] = Q0[i][jm][k][0][1];
																		  Qyd[2] = Q0[i][jm][k][0][2];
																		  Qyd[3] = Q0[i][jm][k][1][1];
																		  Qyd[4] = Q0[i][jm][k][1][2];


																		  Qzu[0] = Q0[i][j][kp][0][0];
																		  Qzu[1] = Q0[i][j][kp][0][1];
																		  Qzu[2] = Q0[i][j][kp][0][2];
																		  Qzu[3] = Q0[i][j][kp][1][1];
																		  Qzu[4] = Q0[i][j][kp][1][2];

																		  if(k>0){
																					 Qzd[0] = Q0[i][j][km][0][0];
																					 Qzd[1] = Q0[i][j][km][0][1];
																					 Qzd[2] = Q0[i][j][km][0][2];
																					 Qzd[3] = Q0[i][j][km][1][1];
																					 Qzd[4] = Q0[i][j][km][1][2];

																		  }
																		  else{
																					 Qzd[0] = Qo_bottom[0];
																					 Qzd[1] = Qo_bottom[1];
																					 Qzd[2] = Qo_bottom[2];
																					 Qzd[3] = Qo_bottom[3];
																					 Qzd[4] = Qo_bottom[4];

																		  }

																		  //x-first and second derivatives
																		  d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
																		  d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
																		  d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
																		  d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
																		  d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
																		  t10  = 0.5 * (Qxu[0]-Qxd[0]);
																		  t20  = 0.5 * (Qxu[1]-Qxd[1]);
																		  t30  = 0.5 * (Qxu[2]-Qxd[2]);
																		  t40  = 0.5 * (Qxu[3]-Qxd[3]);
																		  t50  = 0.5 * (Qxu[4]-Qxd[4]);

																		  //y-first and second derivatives
																		  d2y0 = eight3rd*Qyu[0]+four3rd*Qyd[0]-4.0*a0;
																		  d2y1 = eight3rd*Qyu[1]+four3rd*Qyd[1]-4.0*a1;
																		  d2y2 = eight3rd*Qyu[2]+four3rd*Qyd[2]-4.0*a2;
																		  d2y3 = eight3rd*Qyu[3]+four3rd*Qyd[3]-4.0*a3;
																		  d2y4 = eight3rd*Qyu[4]+four3rd*Qyd[4]-4.0*a4;
																		  t11  = four3rd*Qyu[0]-third*Qyd[0]-a0;
																		  t21  = four3rd*Qyu[1]-third*Qyd[1]-a1;
																		  t31  = four3rd*Qyu[2]-third*Qyd[2]-a2;
																		  t41  = four3rd*Qyu[3]-third*Qyd[3]-a3;
																		  t51  = four3rd*Qyu[4]-third*Qyd[4]-a4;


																		  //z-first and second derivatives
																		  if(k>0){
																					 d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																					 d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																					 d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																					 d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																					 d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																					 t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																					 t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																					 t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																					 t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																					 t52  = 0.5 * ( Qzu[4] - Qzd[4] );
																		  }
																		  else{
																					 d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																					 d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																					 d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																					 d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																					 d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;
																					 t12  = -four3rd*Qzd[0] + third*Qzu[0]   +     a0;
																					 t22  = -four3rd*Qzd[1] + third*Qzu[1]   +     a1;
																					 t32  = -four3rd*Qzd[2] + third*Qzu[2]   +     a2;
																					 t42  = -four3rd*Qzd[3] + third*Qzu[3]   +     a3;
																					 t52  = -four3rd*Qzd[4] + third*Qzu[4]   +     a4;
																		  }
																}

																//if on right side of cavity
																if(j>halfNy && k<h){

																		  Qxu[0] = Q0[ip][j][k][0][0];
																		  Qxu[1] = Q0[ip][j][k][0][1];
																		  Qxu[2] = Q0[ip][j][k][0][2];
																		  Qxu[3] = Q0[ip][j][k][1][1];
																		  Qxu[4] = Q0[ip][j][k][1][2];

																		  Qxd[0] = Q0[im][j][k][0][0];
																		  Qxd[1] = Q0[im][j][k][0][1];
																		  Qxd[2] = Q0[im][j][k][0][2];
																		  Qxd[3] = Q0[im][j][k][1][1];
																		  Qxd[4] = Q0[im][j][k][1][2];

																		  Qyu[0] = Q0[i][jp][k][0][0];
																		  Qyu[1] = Q0[i][jp][k][0][1];
																		  Qyu[2] = Q0[i][jp][k][0][2];
																		  Qyu[3] = Q0[i][jp][k][1][1];
																		  Qyu[4] = Q0[i][jp][k][1][2];

																		  Qyd[0] = Qsides[0][0];
																		  Qyd[1] = Qsides[0][1];
																		  Qyd[2] = Qsides[0][2];
																		  Qyd[3] = Qsides[1][1];
																		  Qyd[4] = Qsides[1][2];

																		  Qzu[0] = Q0[i][j][kp][0][0];
																		  Qzu[1] = Q0[i][j][kp][0][1];
																		  Qzu[2] = Q0[i][j][kp][0][2];
																		  Qzu[3] = Q0[i][j][kp][1][1];
																		  Qzu[4] = Q0[i][j][kp][1][2];

																		  if(k>0){
																					 Qzd[0] = Q0[i][j][km][0][0];
																					 Qzd[1] = Q0[i][j][km][0][1];
																					 Qzd[2] = Q0[i][j][km][0][2];
																					 Qzd[3] = Q0[i][j][km][1][1];
																					 Qzd[4] = Q0[i][j][km][1][2];
																		  }
																		  else{
																					 Qzd[0] = Qo_bottom[0];
																					 Qzd[1] = Qo_bottom[1];
																					 Qzd[2] = Qo_bottom[2];
																					 Qzd[3] = Qo_bottom[3];
																					 Qzd[4] = Qo_bottom[4];

																		  }

																		  //x-first and second derivatives
																		  d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
																		  d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
																		  d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
																		  d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
																		  d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
																		  t10  = 0.5 * (Qxu[0]-Qxd[0]);
																		  t20  = 0.5 * (Qxu[1]-Qxd[1]);
																		  t30  = 0.5 * (Qxu[2]-Qxd[2]);
																		  t40  = 0.5 * (Qxu[3]-Qxd[3]);
																		  t50  = 0.5 * (Qxu[4]-Qxd[4]);

																		  //y-first and second derivatives
																		  d2y0 = eight3rd*Qyd[0]+four3rd*Qyu[0]-4.0*a0;
																		  d2y1 = eight3rd*Qyd[1]+four3rd*Qyu[1]-4.0*a1;
																		  d2y2 = eight3rd*Qyd[2]+four3rd*Qyu[2]-4.0*a2;
																		  d2y3 = eight3rd*Qyd[3]+four3rd*Qyu[3]-4.0*a3;
																		  d2y4 = eight3rd*Qyd[4]+four3rd*Qyu[4]-4.0*a4;
																		  t11  = -four3rd*Qyd[0]+third*Qyu[0]+a0;
																		  t21  = -four3rd*Qyd[1]+third*Qyu[1]+a1;
																		  t31  = -four3rd*Qyd[2]+third*Qyu[2]+a2;
																		  t41  = -four3rd*Qyd[3]+third*Qyu[3]+a3;
																		  t51  = -four3rd*Qyd[4]+third*Qyu[4]+a4;


																		  //z-first and second derivatives
																		  if(k>0){
																					 d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																					 d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																					 d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																					 d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																					 d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																					 t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																					 t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																					 t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																					 t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																					 t52  = 0.5 * ( Qzu[4] - Qzd[4] );
																		  }
																		  else{
																					 d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																					 d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																					 d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																					 d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																					 d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;

																					 t12  = -four3rd*Qzd[0]+third*Qzu[0]+a0;
																					 t22  = -four3rd*Qzd[1]+third*Qzu[1]+a1;
																					 t32  = -four3rd*Qzd[2]+third*Qzu[2]+a2;
																					 t42  = -four3rd*Qzd[3]+third*Qzu[3]+a3;
																					 t52  = -four3rd*Qzd[4]+third*Qzu[4]+a4;
																		  }

																}

																//If on top of cavity
																if(k==h){

																		  Qxu[0] = Q0[ip][j][k][0][0];
																		  Qxu[1] = Q0[ip][j][k][0][1];
																		  Qxu[2] = Q0[ip][j][k][0][2];
																		  Qxu[3] = Q0[ip][j][k][1][1];
																		  Qxu[4] = Q0[ip][j][k][1][2];

																		  Qxd[0] = Q0[im][j][k][0][0];
																		  Qxd[1] = Q0[im][j][k][0][1];
																		  Qxd[2] = Q0[im][j][k][0][2];
																		  Qxd[3] = Q0[im][j][k][1][1];
																		  Qxd[4] = Q0[im][j][k][1][2];

																		  Qyu[0] = Q0[i][jp][k][0][0];
																		  Qyu[1] = Q0[i][jp][k][0][1];
																		  Qyu[2] = Q0[i][jp][k][0][2];
																		  Qyu[3] = Q0[i][jp][k][1][1];
																		  Qyu[4] = Q0[i][jp][k][1][2];

																		  Qyd[0] = Q0[i][jm][k][0][0];
																		  Qyd[1] = Q0[i][jm][k][0][1];
																		  Qyd[2] = Q0[i][jm][k][0][2];
																		  Qyd[3] = Q0[i][jm][k][1][1];
																		  Qyd[4] = Q0[i][jm][k][1][2];

																		  Qzu[0] = Q0[i][j][kp][0][0];
																		  Qzu[1] = Q0[i][j][kp][0][1];
																		  Qzu[2] = Q0[i][j][kp][0][2];
																		  Qzu[3] = Q0[i][j][kp][1][1];
																		  Qzu[4] = Q0[i][j][kp][1][2];

																		  if(cav_flag[i][j][km]==-1){
																					 Qzd[0] = Qtop[0][0];
																					 Qzd[1] = Qtop[0][1];
																					 Qzd[2] = Qtop[0][2];
																					 Qzd[3] = Qtop[1][1];
																					 Qzd[4] = Qtop[1][2];
																		  }
																		  else {
																					 Qzd[0] = Q0[i][j][km][0][0];
																					 Qzd[1] = Q0[i][j][km][0][1];
																					 Qzd[2] = Q0[i][j][km][0][2];
																					 Qzd[3] = Q0[i][j][km][1][1];
																					 Qzd[4] = Q0[i][j][km][1][2];
																		  }


																		  //x-first and second derivatives
																		  d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
																		  d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
																		  d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
																		  d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
																		  d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
																		  t10  = 0.5 * (Qxu[0]-Qxd[0]);
																		  t20  = 0.5 * (Qxu[1]-Qxd[1]);
																		  t30  = 0.5 * (Qxu[2]-Qxd[2]);
																		  t40  = 0.5 * (Qxu[3]-Qxd[3]);
																		  t50  = 0.5 * (Qxu[4]-Qxd[4]);

																		  //y-first and second derivatives
																		  d2y0 = Qyd[0]+Qyd[0]-2.0*a0;
																		  d2y1 = Qyd[1]+Qyd[1]-2.0*a1;
																		  d2y2 = Qyd[2]+Qyd[2]-2.0*a2;
																		  d2y3 = Qyd[3]+Qyd[3]-2.0*a3;
																		  d2y4 = Qyd[4]+Qyd[4]-2.0*a4;
																		  t11  = 0.5 * (Qyu[0]-Qyd[0]);
																		  t21  = 0.5 * (Qyu[1]-Qyd[1]);
																		  t31  = 0.5 * (Qyu[2]-Qyd[2]);
																		  t41  = 0.5 * (Qyu[3]-Qyd[3]);
																		  t51  = 0.5 * (Qyu[4]-Qyd[4]);

																		  //z-first and second derivatives
																		  if(cav_flag[i][j][km]>=0){
																					 d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																					 d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																					 d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																					 d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																					 d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																					 t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																					 t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																					 t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																					 t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																					 t52  = 0.5 * ( Qzu[4] - Qzd[4] );
																		  }
																		  else {
																					 d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																					 d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																					 d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																					 d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																					 d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;
																					 t12  = -four3rd*Qzd[0] + third*Qzu[0] + a0;
																					 t22  = -four3rd*Qzd[1] + third*Qzu[1] + a1;
																					 t32  = -four3rd*Qzd[2] + third*Qzu[2] + a2;
																					 t42  = -four3rd*Qzd[3] + third*Qzu[3] + a3;
																					 t52  = -four3rd*Qzd[4] + third*Qzu[4] + a4;
																		  }

																}
													 }
										  }

										  //lapacian Q
										  d2Q[0][0] = d2x0 + d2y0 + d2z0;
										  d2Q[0][1] = d2x1 + d2y1 + d2z1;
										  d2Q[0][2] = d2x2 + d2y2 + d2z2;
										  d2Q[1][1] = d2x3 + d2y3 + d2z3;
										  d2Q[1][2] = d2x4 + d2y4 + d2z4;
										  d2Q[1][0] = d2Q[0][1];
										  d2Q[2][0] = d2Q[0][2];
										  d2Q[2][1] = d2Q[1][2];
										  d2Q[2][2] =-d2Q[0][0] - d2Q[1][1];

										  //lap_Q[i][j][k]=2.0*(d2Q[0][1]+d2Q[0][2]+d2Q[1][2]);

										  //	calculate H				
										  trQ2 = 2.0*(a0*a0+a1*a1+a2*a2+a3*a3+a4*a4+a0*a3);

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																temp = Q0[i][j][k][ii][jj];
																H0[i][j][k][ii][jj] = -1*(A_ldg*(1.0-third*U)*temp - A_ldg*U*( QQ(Q0[i][j][k],ii,jj) - trQ2*(temp + OneThirdDelta(ii,jj))) - kappa * d2Q[ii][jj]);
													 }
										  }

										  //check trace of H
										  if(debug_mode!=0){
													 trH3rd = third*(H0[i][j][k][0][0]+H0[i][j][k][1][1]+H0[i][j][k][2][2]);
													 if (trH3rd>1e-15 || trH3rd<-1e-15)printf("trH3rd nonzero: %d %d %d: %le\n",i, j, k,trH3rd);
										  }

										  //	calculate convective Q
										  ux = u[i][j][k][0];
										  uy = u[i][j][k][1];
										  uz = u[i][j][k][2];
										  convQ[i][j][k][0] = ux * t10 + uy * t11 + uz * t12;
										  convQ[i][j][k][1] = ux * t20 + uy * t21 + uz * t22;
										  convQ[i][j][k][2] = ux * t30 + uy * t31 + uz * t32;
										  convQ[i][j][k][3] = ux * t40 + uy * t41 + uz * t42;
										  convQ[i][j][k][4] = ux * t50 + uy * t51 + uz * t52;
								}
					 }
		  }
}


// calculate Q tensor related sigma
void cal_stress_q(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real sigma_q[Nx][Ny][Nz][3][3])
{

		  int i, j, k, ii=0, jj=0, kk=0, mm;
		  int ip, im, jp, jm, kp, km;
		  int ipar;
		  real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
		  real d2x0, d2x1, d2x2, d2x3, d2x4, d2y0, d2y1, d2y2, d2y3, d2y4, d2z0, d2z1, d2z2, d2z3, d2z4;
		  real temp, temp2, Qu[5], Qd[5];
		  real dx, dy, dz, dr2i;
		  real a0, a1, a2, a3, a4;
		  real trQ2, trH3rd, d2Q[3][3], dQ[3][3], dQ2;
		  real sum;

		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {
								for (k=0; k<Nz; k++) {

										  a0 = Q0[i][j][k][0][0];
										  a1 = Q0[i][j][k][0][1];
										  a2 = Q0[i][j][k][0][2];
										  a3 = Q0[i][j][k][1][1];
										  a4 = Q0[i][j][k][1][2];

										  im = i-1;
										  jm = j-1;
										  ip = i+1;
										  jp = j+1;
										  km = k-1;
										  kp = k+1;

										  if(i==0) im = Nx-1;
										  if(j==0) jm = Ny-1;
										  if(k==0) km = Nz-1;

										  if(i==Nx-1) ip = 0;
										  if(j==Ny-1) jp = 0;
										  if(k==Nz-1) kp = 0;

										  //Bulk points away from cavity
										  if(cav_flag[i][j][k]==0 || cav_on!=0){

													 if(k==0 && (wall_anchoring!=0 && W_bot<0) ){
																Qd[0] = Qo_bottom[0];
																Qd[1] = Qo_bottom[1];
																Qd[2] = Qo_bottom[2];
																Qd[3] = Qo_bottom[3];
																Qd[4] = Qo_bottom[4];
													 }
													 else if(k==0 && (wall_anchoring!=0 && W_bot<=1) ){
																Qd[0] = Qsurf[i][j][0][0];
																Qd[1] = Qsurf[i][j][0][1];
																Qd[2] = Qsurf[i][j][0][2];
																Qd[3] = Qsurf[i][j][0][3];
																Qd[4] = Qsurf[i][j][0][4];
													 }
													 else {
																Qd[0] = Q0[i][j][km][0][0];
																Qd[1] = Q0[i][j][km][0][1];
																Qd[2] = Q0[i][j][km][0][2];
																Qd[3] = Q0[i][j][km][1][1];
																Qd[4] = Q0[i][j][km][1][2];
													 }

													 if(k==Nz-1 && (wall_anchoring !=0 && W_top<0) ) {
																Qu[0] = Qo_top[0];
																Qu[1] = Qo_top[1];
																Qu[2] = Qo_top[2];
																Qu[3] = Qo_top[3];
																Qu[4] = Qo_top[4];
													 } 
													 else if(k==Nz-1 && (wall_anchoring !=0 && W_top>=0) ) {
																Qu[0] = Qsurf[i][j][1][0];
																Qu[1] = Qsurf[i][j][1][1];
																Qu[2] = Qsurf[i][j][1][2];
																Qu[3] = Qsurf[i][j][1][3];
																Qu[4] = Qsurf[i][j][1][4];
													 } 
													 else {
																Qu[0] = Q0[i][j][kp][0][0];
																Qu[1] = Q0[i][j][kp][0][1];
																Qu[2] = Q0[i][j][kp][0][2];
																Qu[3] = Q0[i][j][kp][1][1];
																Qu[4] = Q0[i][j][kp][1][2];
													 }

													 //x-first and second derivatives
													 d2x0 = Q0[ip][j][k][0][0] + Q0[im][j][k][0][0] - 2.0*a0;
													 d2x1 = Q0[ip][j][k][0][1] + Q0[im][j][k][0][1] - 2.0*a1;
													 d2x2 = Q0[ip][j][k][0][2] + Q0[im][j][k][0][2] - 2.0*a2;
													 d2x3 = Q0[ip][j][k][1][1] + Q0[im][j][k][1][1] - 2.0*a3;
													 d2x4 = Q0[ip][j][k][1][2] + Q0[im][j][k][1][2] - 2.0*a4;
													 t10  = 0.5 * (Q0[ip][j][k][0][0]-Q0[im][j][k][0][0]);
													 t20  = 0.5 * (Q0[ip][j][k][0][1]-Q0[im][j][k][0][1]);
													 t30  = 0.5 * (Q0[ip][j][k][0][2]-Q0[im][j][k][0][2]);
													 t40  = 0.5 * (Q0[ip][j][k][1][1]-Q0[im][j][k][1][1]);
													 t50  = 0.5 * (Q0[ip][j][k][1][2]-Q0[im][j][k][1][2]);

													 //y-first and second derivatives
													 d2y0 = Q0[i][jm][k][0][0]+Q0[i][jp][k][0][0]-2.0*a0;
													 d2y1 = Q0[i][jm][k][0][1]+Q0[i][jp][k][0][1]-2.0*a1;
													 d2y2 = Q0[i][jm][k][0][2]+Q0[i][jp][k][0][2]-2.0*a2;
													 d2y3 = Q0[i][jm][k][1][1]+Q0[i][jp][k][1][1]-2.0*a3;
													 d2y4 = Q0[i][jm][k][1][2]+Q0[i][jp][k][1][2]-2.0*a4;
													 t11  = 0.5 * (Q0[i][jp][k][0][0]-Q0[i][jm][k][0][0]);
													 t21  = 0.5 * (Q0[i][jp][k][0][1]-Q0[i][jm][k][0][1]);
													 t31  = 0.5 * (Q0[i][jp][k][0][2]-Q0[i][jm][k][0][2]);
													 t41  = 0.5 * (Q0[i][jp][k][1][1]-Q0[i][jm][k][1][1]);
													 t51  = 0.5 * (Q0[i][jp][k][1][2]-Q0[i][jm][k][1][2]);

													 if(k==0 && wall_anchoring !=0 ){
																d2z0 = eight3rd*Qd[0] + four3rd*Qu[0] - 4.0 * a0;
																d2z1 = eight3rd*Qd[1] + four3rd*Qu[1] - 4.0 * a1;
																d2z2 = eight3rd*Qd[2] + four3rd*Qu[2] - 4.0 * a2;
																d2z3 = eight3rd*Qd[3] + four3rd*Qu[3] - 4.0 * a3;
																d2z4 = eight3rd*Qd[4] + four3rd*Qu[4] - 4.0 * a4;
																t12  =-four3rd*Qd[0] + third*Qu[0] +       a0;
																t22  =-four3rd*Qd[1] + third*Qu[1] +       a1;
																t32  =-four3rd*Qd[2] + third*Qu[2] +       a2;
																t42  =-four3rd*Qd[3] + third*Qu[3] +       a3;
																t52  =-four3rd*Qd[4] + third*Qu[4] +       a4;
													 } 

													 else if(k==Nz-1 && wall_anchoring!=0){
																d2z0 = eight3rd*Qu[0] + four3rd*Qd[0] - 4.0 * a0;
																d2z1 = eight3rd*Qu[1] + four3rd*Qd[1] - 4.0 * a1;
																d2z2 = eight3rd*Qu[2] + four3rd*Qd[2] - 4.0 * a2;
																d2z3 = eight3rd*Qu[3] + four3rd*Qd[3] - 4.0 * a3;
																d2z4 = eight3rd*Qu[4] + four3rd*Qd[4] - 4.0 * a4;
																t12  = four3rd*Qu[0] - third*Qd[0] -       a0;
																t22  = four3rd*Qu[1] - third*Qd[1] -       a1;
																t32  = four3rd*Qu[2] - third*Qd[2] -       a2;
																t42  = four3rd*Qu[3] - third*Qd[3] -       a3;
																t52  = four3rd*Qu[4] - third*Qd[4] -       a4;
													 }
													 else {
																d2z0 = Qd[0] + Qu[0] - 2.0*a0;
																d2z1 = Qd[1] + Qu[1] - 2.0*a1;
																d2z2 = Qd[2] + Qu[2] - 2.0*a2;
																d2z3 = Qd[3] + Qu[3] - 2.0*a3;
																d2z4 = Qd[4] + Qu[4] - 2.0*a4;
																t12  = 0.5 * ( Qu[0] - Qd[0] );
																t22  = 0.5 * ( Qu[1] - Qd[1] );
																t32  = 0.5 * ( Qu[2] - Qd[2] );
																t42  = 0.5 * ( Qu[3] - Qd[3] );
																t52  = 0.5 * ( Qu[4] - Qd[4] );
													 }
										  }

										  if(cav_flag[i][j][k]==1 && cav_on==1){

													 //cavity Q's for derivatives (u=up, d=down)
													 int halfNy = 0.5*Ny;
													 int h = cav_height;
													 real Qxu[5]={0}, Qxd[5]={0}, Qyu[5]={0}, Qyd[5]={0}, Qzu[5]={0}, Qzd[5]={0};

													 //if on left side of cavity
													 if(j<halfNy && k<h){
																Qxu[0] = Q0[ip][j][k][0][0];
																Qxu[1] = Q0[ip][j][k][0][1];
																Qxu[2] = Q0[ip][j][k][0][2];
																Qxu[3] = Q0[ip][j][k][1][1];
																Qxu[4] = Q0[ip][j][k][1][2];

																Qxd[0] = Q0[im][j][k][0][0];
																Qxd[1] = Q0[im][j][k][0][1];
																Qxd[2] = Q0[im][j][k][0][2];
																Qxd[3] = Q0[im][j][k][1][1];
																Qxd[4] = Q0[im][j][k][1][2];

																Qyd[0] = Q0[i][jm][k][0][0];
																Qyd[1] = Q0[i][jm][k][0][1];
																Qyd[2] = Q0[i][jm][k][0][2];
																Qyd[3] = Q0[i][jm][k][1][1];
																Qyd[4] = Q0[i][jm][k][1][2];

																Qyu[0] = Qsides[0][0];
																Qyu[1] = Qsides[0][1];
																Qyu[2] = Qsides[0][2];
																Qyu[3] = Qsides[1][1];
																Qyu[4] = Qsides[1][2];

																Qzu[0] = Q0[i][j][kp][0][0];
																Qzu[1] = Q0[i][j][kp][0][1];
																Qzu[2] = Q0[i][j][kp][0][2];
																Qzu[3] = Q0[i][j][kp][1][1];
																Qzu[4] = Q0[i][j][kp][1][2];

																if(k>0){
																		  Qzd[0] = Q0[i][j][km][0][0];
																		  Qzd[1] = Q0[i][j][km][0][1];
																		  Qzd[2] = Q0[i][j][km][0][2];
																		  Qzd[3] = Q0[i][j][km][1][1];
																		  Qzd[4] = Q0[i][j][km][1][2];

																}
																else{
																		  Qzd[0] = Qsurf[i][j][0][0];
																		  Qzd[1] = Qsurf[i][j][0][1];
																		  Qzd[2] = Qsurf[i][j][0][2];
																		  Qzd[3] = Qsurf[i][j][0][3];
																		  Qzd[4] = Qsurf[i][j][0][4];

																}

																//x-first and second derivatives
																d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
																d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
																d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
																d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
																d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
																t10  = 0.5 * (Qxu[0]-Qxd[0]);
																t20  = 0.5 * (Qxu[1]-Qxd[1]);
																t30  = 0.5 * (Qxu[2]-Qxd[2]);
																t40  = 0.5 * (Qxu[3]-Qxd[3]);
																t50  = 0.5 * (Qxu[4]-Qxd[4]);

																//y-first and second derivatives
																d2y0 = eight3rd*Qyu[0]+four3rd*Qyd[0]-4.0*a0;
																d2y1 = eight3rd*Qyu[1]+four3rd*Qyd[1]-4.0*a1;
																d2y2 = eight3rd*Qyu[2]+four3rd*Qyd[2]-4.0*a2;
																d2y3 = eight3rd*Qyu[3]+four3rd*Qyd[3]-4.0*a3;
																d2y4 = eight3rd*Qyu[4]+four3rd*Qyd[4]-4.0*a4;
																t11  = four3rd*Qyu[0]-third*Qyd[0]-a0;
																t21  = four3rd*Qyu[1]-third*Qyd[1]-a1;
																t31  = four3rd*Qyu[2]-third*Qyd[2]-a2;
																t41  = four3rd*Qyu[3]-third*Qyd[3]-a3;
																t51  = four3rd*Qyu[4]-third*Qyd[4]-a4;


																//z-first and second derivatives
																if(k>0){
																		  d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																		  d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																		  d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																		  d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																		  d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																		  t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																		  t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																		  t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																		  t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																		  t52  = 0.5 * ( Qzu[4] - Qzd[4] );
																}
																else{
																		  d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																		  d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																		  d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																		  d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																		  d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;
																		  t12  = -four3rd*Qzd[0]+third*Qzu[0]+a0;
																		  t22  = -four3rd*Qzd[1]+third*Qzu[1]+a1;
																		  t32  = -four3rd*Qzd[2]+third*Qzu[2]+a2;
																		  t42  = -four3rd*Qzd[3]+third*Qzu[3]+a3;
																		  t52  = -four3rd*Qzd[4]+third*Qzu[4]+a4;
																}
													 }

													 //if on right side of cavity
													 if(j>halfNy && k<h){
																Qxu[0] = Q0[ip][j][k][0][0];
																Qxu[1] = Q0[ip][j][k][0][1];
																Qxu[2] = Q0[ip][j][k][0][2];
																Qxu[3] = Q0[ip][j][k][1][1];
																Qxu[4] = Q0[ip][j][k][1][2];

																Qxd[0] = Q0[im][j][k][0][0];
																Qxd[1] = Q0[im][j][k][0][1];
																Qxd[2] = Q0[im][j][k][0][2];
																Qxd[3] = Q0[im][j][k][1][1];
																Qxd[4] = Q0[im][j][k][1][2];

																Qyu[0] = Q0[i][jp][k][0][0];
																Qyu[1] = Q0[i][jp][k][0][1];
																Qyu[2] = Q0[i][jp][k][0][2];
																Qyu[3] = Q0[i][jp][k][1][1];
																Qyu[4] = Q0[i][jp][k][1][2];

																Qyd[0] = Qsides[0][0];
																Qyd[1] = Qsides[0][1];
																Qyd[2] = Qsides[0][2];
																Qyd[3] = Qsides[1][1];
																Qyd[4] = Qsides[1][2];

																Qzu[0] = Q0[i][j][kp][0][0];
																Qzu[1] = Q0[i][j][kp][0][1];
																Qzu[2] = Q0[i][j][kp][0][2];
																Qzu[3] = Q0[i][j][kp][1][1];
																Qzu[4] = Q0[i][j][kp][1][2];

																if(k>0){
																		  Qzd[0] = Q0[i][j][km][0][0];
																		  Qzd[1] = Q0[i][j][km][0][1];
																		  Qzd[2] = Q0[i][j][km][0][2];
																		  Qzd[3] = Q0[i][j][km][1][1];
																		  Qzd[4] = Q0[i][j][km][1][2];
																}
																else{
																		  Qzd[0] = Qsurf[i][j][0][0];
																		  Qzd[1] = Qsurf[i][j][0][1];
																		  Qzd[2] = Qsurf[i][j][0][2];
																		  Qzd[3] = Qsurf[i][j][0][3];
																		  Qzd[4] = Qsurf[i][j][0][4];

																}

																//x-first and second derivatives
																d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
																d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
																d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
																d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
																d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
																t10  = 0.5 * (Qxu[0]-Qxd[0]);
																t20  = 0.5 * (Qxu[1]-Qxd[1]);
																t30  = 0.5 * (Qxu[2]-Qxd[2]);
																t40  = 0.5 * (Qxu[3]-Qxd[3]);
																t50  = 0.5 * (Qxu[4]-Qxd[4]);

																//y-first and second derivatives
																d2y0 = eight3rd*Qyd[0]+four3rd*Qyu[0]-4.0*a0;
																d2y1 = eight3rd*Qyd[1]+four3rd*Qyu[1]-4.0*a1;
																d2y2 = eight3rd*Qyd[2]+four3rd*Qyu[2]-4.0*a2;
																d2y3 = eight3rd*Qyd[3]+four3rd*Qyu[3]-4.0*a3;
																d2y4 = eight3rd*Qyd[4]+four3rd*Qyu[4]-4.0*a4;
																t11  = -four3rd*Qyd[0]+third*Qyu[0]+a0;
																t21  = -four3rd*Qyd[1]+third*Qyu[1]+a1;
																t31  = -four3rd*Qyd[2]+third*Qyu[2]+a2;
																t41  = -four3rd*Qyd[3]+third*Qyu[3]+a3;
																t51  = -four3rd*Qyd[4]+third*Qyu[4]+a4;


																//z-first and second derivatives
																if(k>0){
																		  d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																		  d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																		  d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																		  d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																		  d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																		  t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																		  t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																		  t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																		  t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																		  t52  = 0.5 * ( Qzu[4] - Qzd[4] );
																}
																else{
																		  d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																		  d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																		  d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																		  d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																		  d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;
																		  t12  = -four3rd*Qzd[0]+third*Qzu[0]+a0;
																		  t22  = -four3rd*Qzd[1]+third*Qzu[1]+a1;
																		  t32  = -four3rd*Qzd[2]+third*Qzu[2]+a2;
																		  t42  = -four3rd*Qzd[3]+third*Qzu[3]+a3;
																		  t52  = -four3rd*Qzd[4]+third*Qzu[4]+a4;
																}

													 }

													 //If on top of cavity
													 if(k==h){
																Qxu[0] = Q0[ip][j][k][0][0];
																Qxu[1] = Q0[ip][j][k][0][1];
																Qxu[2] = Q0[ip][j][k][0][2];
																Qxu[3] = Q0[ip][j][k][1][1];
																Qxu[4] = Q0[ip][j][k][1][2];

																Qxd[0] = Q0[im][j][k][0][0];
																Qxd[1] = Q0[im][j][k][0][1];
																Qxd[2] = Q0[im][j][k][0][2];
																Qxd[3] = Q0[im][j][k][1][1];
																Qxd[4] = Q0[im][j][k][1][2];

																Qyu[0] = Q0[i][jp][k][0][0];
																Qyu[1] = Q0[i][jp][k][0][1];
																Qyu[2] = Q0[i][jp][k][0][2];
																Qyu[3] = Q0[i][jp][k][1][1];
																Qyu[4] = Q0[i][jp][k][1][2];

																Qyd[0] = Q0[i][jm][k][0][0];
																Qyd[1] = Q0[i][jm][k][0][1];
																Qyd[2] = Q0[i][jm][k][0][2];
																Qyd[3] = Q0[i][jm][k][1][1];
																Qyd[4] = Q0[i][jm][k][1][2];

																Qzu[0] = Q0[i][j][kp][0][0];
																Qzu[1] = Q0[i][j][kp][0][1];
																Qzu[2] = Q0[i][j][kp][0][2];
																Qzu[3] = Q0[i][j][kp][1][1];
																Qzu[4] = Q0[i][j][kp][1][2];

																if(cav_flag[i][j][km]==-1){
																		  Qzd[0] = Qtop[0][0];
																		  Qzd[1] = Qtop[0][1];
																		  Qzd[2] = Qtop[0][2];
																		  Qzd[3] = Qtop[1][1];
																		  Qzd[4] = Qtop[1][2];
																}
																else {
																		  Qzd[0] = Qtop[0][0];
																		  Qzd[1] = Qtop[0][1];
																		  Qzd[2] = Qtop[0][2];
																		  Qzd[3] = Qtop[1][1];
																		  Qzd[4] = Qtop[1][2];
																}
													 }

													 //x-first and second derivatives
													 d2x0 = Qxu[0] + Qxd[0] - 2.0*a0;
													 d2x1 = Qxu[1] + Qxd[1] - 2.0*a1;
													 d2x2 = Qxu[2] + Qxd[2] - 2.0*a2;
													 d2x3 = Qxu[3] + Qxd[3] - 2.0*a3;
													 d2x4 = Qxu[4] + Qxd[4] - 2.0*a4;
													 t10  = 0.5 * (Qxu[0]-Qxd[0]);
													 t20  = 0.5 * (Qxu[1]-Qxd[1]);
													 t30  = 0.5 * (Qxu[2]-Qxd[2]);
													 t40  = 0.5 * (Qxu[3]-Qxd[3]);
													 t50  = 0.5 * (Qxu[4]-Qxd[4]);

													 //y-first and second derivatives
													 d2y0 = Qyd[0]+Qyd[0]-2.0*a0;
													 d2y1 = Qyd[1]+Qyd[1]-2.0*a1;
													 d2y2 = Qyd[2]+Qyd[2]-2.0*a2;
													 d2y3 = Qyd[3]+Qyd[3]-2.0*a3;
													 d2y4 = Qyd[4]+Qyd[4]-2.0*a4;
													 t11  = 0.5 * (Qyu[0]-Qyd[0]);
													 t21  = 0.5 * (Qyu[1]-Qyd[1]);
													 t31  = 0.5 * (Qyu[2]-Qyd[2]);
													 t41  = 0.5 * (Qyu[3]-Qyd[3]);
													 t51  = 0.5 * (Qyu[4]-Qyd[4]);

													 //z-first and second derivatives
													 if(cav_flag[i][j][km]==1){
																d2z0 = Qzd[0] + Qzu[0] - 2.0*a0;
																d2z1 = Qzd[1] + Qzu[1] - 2.0*a1;
																d2z2 = Qzd[2] + Qzu[2] - 2.0*a2;
																d2z3 = Qzd[3] + Qzu[3] - 2.0*a3;
																d2z4 = Qzd[4] + Qzu[4] - 2.0*a4;
																t12  = 0.5 * ( Qzu[0] - Qzd[0] );
																t22  = 0.5 * ( Qzu[1] - Qzd[1] );
																t32  = 0.5 * ( Qzu[2] - Qzd[2] );
																t42  = 0.5 * ( Qzu[3] - Qzd[3] );
																t52  = 0.5 * ( Qzu[4] - Qzd[4] );
													 }
													 else {
																d2z0 = eight3rd*Qzd[0] + four3rd*Qzu[0] - 4.0*a0;
																d2z1 = eight3rd*Qzd[1] + four3rd*Qzu[1] - 4.0*a1;
																d2z2 = eight3rd*Qzd[2] + four3rd*Qzu[2] - 4.0*a2;
																d2z3 = eight3rd*Qzd[3] + four3rd*Qzu[3] - 4.0*a3;
																d2z4 = eight3rd*Qzd[4] + four3rd*Qzu[4] - 4.0*a4;
																t12  = -four3rd*Qzd[0] + third*Qzu[0] + a0;
																t22  = -four3rd*Qzd[1] + third*Qzu[1] + a1;
																t32  = -four3rd*Qzd[2] + third*Qzu[2] + a2;
																t42  = -four3rd*Qzd[3] + third*Qzu[3] + a3;
																t52  = -four3rd*Qzd[4] + third*Qzu[4] + a4;
													 }

										  }

										  dQ[0][0]=2.*(t10*t10+t20*t20+t30*t30+t40*t40+t50*t50+t10*t40);
										  dQ[1][1]=2.*(t11*t11+t21*t21+t31*t31+t41*t41+t51*t51+t11*t41);
										  dQ[2][2]=2.*(t12*t12+t22*t22+t32*t32+t42*t42+t52*t52+t12*t42);
										  dQ[0][1]=2.*(t10*t11+t20*t21+t30*t31+t40*t41+t50*t51)+t10*t41+t11*t40;
										  dQ[0][2]=2.*(t10*t12+t20*t22+t30*t32+t40*t42+t50*t52)+t10*t42+t12*t40;
										  dQ[1][2]=2.*(t11*t12+t21*t22+t31*t32+t41*t42+t51*t52)+t11*t42+t12*t41;
										  dQ[1][0]=dQ[0][1];
										  dQ[2][0]=dQ[0][2];
										  dQ[2][1]=dQ[1][2];

										  dQ2 = dQ[0][0]+dQ[1][1]+dQ[2][2];
										  //										  lap_Q[i][j][k] = dQ2;

										  //lapacian Q
										  d2Q[0][0] = d2x0 + d2y0 + d2z0;
										  d2Q[0][1] = d2x1 + d2y1 + d2z1;
										  d2Q[0][2] = d2x2 + d2y2 + d2z2;
										  d2Q[1][1] = d2x3 + d2y3 + d2z3;
										  d2Q[1][2] = d2x4 + d2y4 + d2z4;
										  d2Q[1][0] = d2Q[0][1];
										  d2Q[2][0] = d2Q[0][2];
										  d2Q[2][1] = d2Q[1][2];
										  d2Q[2][2] =-d2Q[0][0] - d2Q[1][1];

										  //											lap_Q[i][j][k]=2.0*(d2Q[0][1]+d2Q[0][2]+d2Q[1][2]);

										  //	calculate H				
										  trQ2 = 2.0*(a0*a0+a1*a1+a2*a2+a3*a3+a4*a4+a0*a3);

										  int kk,nn;
										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																temp = Q0[i][j][k][ii][jj];
																H0[i][j][k][ii][jj] = -1*(A_ldg*(1.0-third*U)*temp - A_ldg*U*( QQ(Q0[i][j][k],ii,jj) - trQ2*(temp + OneThirdDelta(ii,jj))) - kappa * d2Q[ii][jj]);
													 }
										  }

										  //check trace of H
										  if(debug_mode!=0){
													 trH3rd = third*(H0[i][j][k][0][0]+H0[i][j][k][1][1]+H0[i][j][k][2][2]);
													 if (trH3rd>1e-15 || trH3rd<-1e-15)printf("trH3rd nonzero: %d %d %d: %le\n",i, j, k,trH3rd);
										  }

										  for(sum=0.0,ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																sum = sum + Q0[i][j][k][ii][jj]*H0[i][j][k][ii][jj];
													 }
										  }

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																temp = 0.5*kappa*dQ2*Delta(ii,jj) + 2*xi*(Q0[i][j][k][ii][jj] + OneThirdDelta(ii,jj))*sum;
																temp2= 0;
//																temp = 2.0*xi*(Q0[i][j][k][ii][jj] + OneThirdDelta(ii,jj))*sum;

																for(mm=0;mm<3;mm++){
																		  temp += - xi*H0[i][j][k][ii][mm]*(Q0[i][j][k][mm][jj] + OneThirdDelta(mm,jj)) - xi*(Q0[i][j][k][ii][mm] + OneThirdDelta(ii,mm))*H0[i][j][k][mm][jj] ;
																}

																for(mm=0;mm<3;mm++){
																		  temp2+=Q0[i][j][k][ii][mm]*H0[i][j][k][mm][jj]-H0[i][j][k][ii][mm]*Q0[i][j][k][mm][jj];
																}

																temp += -kappa * dQ[ii][jj];
																if (t_current>t_v_on){
																		  sigma_q[i][j][k][ii][jj] = temp;
																			tau[i][j][k][ii][jj]   = temp2;
																}

																else {
																		  sigma_q[i][j][k][ii][jj] = 0;	// decouple velocity
																}
													 }
										  }
								}
					 }
		  }
}

void cal_detau(real Q0[Nx][Ny][Nz][3][3])
{
		  int i, j, k, ii=0, jj=0, kk=0;
		  int ip, im, jp, jm, kp, km;
		  int ipar, isurf, jsurf;
		  real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
		  real temp, Qu[5], Qd[5];
		  real d2x0, d2x1, d2x2, d2x3, d2x4, d2y0, d2y1, d2y2, d2y3, d2y4, d2z0, d2z1, d2z2, d2z3, d2z4;
		  real dx, dy, dz, dr2i;
		  real a0, a1, a2, a3, a4;
		  real trQ2, trH3rd, d2Q[3][3], dQ2;
		  real d2Q1, d2Q2, d2Q3, d2Q4, d2Q5;
		  real ux, uy, uz;	

		  int halfNy = 0.5*Ny;
		  int h = cav_height;

		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {
								for (k=0; k<Nz; k++) {

										  im = i - 1;
										  jm = j - 1;
										  km = k - 1;

										  ip = i + 1;
										  jp = j + 1;
										  kp = k + 1;

										  if(i==0) im = Nx-1;
										  if(j==0) jm = Ny-1;
										  if(k==0) km = Nz-1;

										  if(i==Nx-1) ip = 0;
										  if(j==Ny-1) jp = 0;
										  if(k==Nz-1) kp = 0;

										  if( (cav_on && cav_flag[i][j][k]==0) || !cav_on){

													 //x- and y-components
													 detau[i][j][k][0] = 0.5 * ( tau[ip][j][k][0][0] - tau[im][j][k][0][0] + tau[i][jp][k][0][1] - tau[i][jm][k][0][1] );
													 detau[i][j][k][1] = 0.5 * ( tau[ip][j][k][1][0] - tau[im][j][k][1][0] + tau[i][jp][k][1][1] - tau[i][jm][k][1][1] );

													 //z-components
													 detau[i][j][k][2] = 0.5 * ( tau[ip][j][k][2][0] - tau[im][j][k][2][0] + tau[i][jp][k][2][1] - tau[i][jm][k][2][1] );

													 if ( k==0 && pbc_z==0) {

																detau[i][j][k][0] += -0.5*tau[i][j][2][0][2] + 2.0*tau[i][j][1][0][2] - 1.5*tau[i][j][0][0][2];
																//					detau[i][j][k][0] += tau[i][j][1][0][2] - tau[i][j][0][0][2];
																detau[i][j][k][1] += -0.5*tau[i][j][2][1][2] + 2.0*tau[i][j][1][1][2] - 1.5*tau[i][j][0][1][2];
																detau[i][j][k][2] += -0.5*tau[i][j][2][2][2] + 2.0*tau[i][j][1][2][2] - 1.5*tau[i][j][0][2][2];

													 } else if ( k==Nz-1 && pbc_z==0) {

																detau[i][j][k][0] += 0.5*tau[i][j][Nz-3][0][2] - 2.0*tau[i][j][Nz-2][0][2] + 1.5*tau[i][j][Nz-1][0][2];
																//					detau[i][j][k][0] += tau[i][j][Nz-1][0][2] - tau[i][j][Nz-2][0][2];
																detau[i][j][k][1] += 0.5*tau[i][j][Nz-3][1][2] - 2.0*tau[i][j][Nz-2][1][2] + 1.5*tau[i][j][Nz-1][1][2];
																detau[i][j][k][2] += 0.5*tau[i][j][Nz-3][2][2] - 2.0*tau[i][j][Nz-2][2][2] + 1.5*tau[i][j][Nz-1][2][2];

													 } else {

																detau[i][j][k][0] += 0.5 * ( tau[i][j][kp][0][2] - tau[i][j][km][0][2] );
																detau[i][j][k][1] += 0.5 * ( tau[i][j][kp][1][2] - tau[i][j][km][1][2] );
																detau[i][j][k][2] += 0.5 * ( tau[i][j][kp][2][2] - tau[i][j][km][2][2] );

													 }
										  }

										  if( cav_on && cav_flag[i][j][k]==1 ) {


													 //left side of feature
													 if(j<halfNy && k<h){

																//x-comp
																detau[i][j][k][0] = 0.5 * ( tau[ip][j][k][0][0] - tau[im][j][k][0][0] );
																detau[i][j][k][0] += 0.5*tau[i][j-2][k][0][1] - 2*tau[i][j-1][k][0][1]+1.5*tau[i][j][k][0][1];

																//y-comp
																detau[i][j][k][1] = 0.5 * ( tau[ip][j][k][1][0] - tau[im][j][k][1][0] );
																detau[i][j][k][1] += 0.5*tau[i][j-2][k][1][1] - 2*tau[i][j-1][k][1][1]+1.5*tau[i][j][k][1][1];
																
																//z-comp
																detau[i][j][k][2]  =  0.5 * ( tau[ip][j][k][2][0] - tau[im][j][k][2][0] );
																detau[i][j][k][2] +=  0.5*tau[i][j-2][k][2][1] - 2.0*tau[i][j-1][k][2][1] + 1.5*tau[i][j][k][2][1];

																//modifications for z derivative 
																if(k>0){

																		  detau[i][j][k][0] += 0.5 * ( tau[i][j][kp][0][2] - tau[i][j][km][0][2] );
																		  detau[i][j][k][1] += 0.5 * ( tau[i][j][kp][1][2] - tau[i][j][km][1][2] );
																		  detau[i][j][k][2] += 0.5 * ( tau[i][j][kp][2][2] - tau[i][j][km][2][2] );
																} 
																else{

																		  detau[i][j][k][0] += -0.5*tau[i][j][2][0][2] + 2.0*tau[i][j][1][0][2] - 1.5*tau[i][j][0][0][2];
																		  detau[i][j][k][1] += -0.5*tau[i][j][2][1][2] + 2.0*tau[i][j][1][1][2] - 1.5*tau[i][j][0][1][2];
																		  detau[i][j][k][2] += -0.5*tau[i][j][2][2][2] + 2.0*tau[i][j][1][2][2] - 1.5*tau[i][j][0][2][2];

																}

													 }


													 //right side of feature
													 if(j>halfNy && k<h){

																//x-comp
																detau[i][j][k][0] = 0.5 * ( tau[ip][j][k][0][0] - tau[im][j][k][0][0] );
																detau[i][j][k][0] += -0.5*tau[i][j+2][k][0][1] + 2*tau[i][j+1][k][0][1] - 1.5*tau[i][j][k][0][1];

																//y-comp
																detau[i][j][k][1] = 0.5 * ( tau[ip][j][k][1][0] - tau[im][j][k][1][0] );
																detau[i][j][k][1] += -0.5*tau[i][j+2][k][1][1] + 2*tau[i][j+1][k][1][1] - 1.5*tau[i][j][k][1][1];
																
																//z-comp
																detau[i][j][k][2]  =  0.5 * ( tau[ip][j][k][2][0] - tau[im][j][k][2][0] );
																detau[i][j][k][2] +=  -0.5*tau[i][j+2][k][2][1] + 2.0*tau[i][j+1][k][2][1] - 1.5*tau[i][j][k][2][1];

																//modifications for z derivative 
																if(k>0){

																		  detau[i][j][k][0] += 0.5 * ( tau[i][j][kp][0][2] - tau[i][j][km][0][2] );
																		  detau[i][j][k][1] += 0.5 * ( tau[i][j][kp][1][2] - tau[i][j][km][1][2] );
																		  detau[i][j][k][2] += 0.5 * ( tau[i][j][kp][2][2] - tau[i][j][km][2][2] );
																} 
																else{

																		  detau[i][j][k][0] += -0.5*tau[i][j][2][0][2] + 2.0*tau[i][j][1][0][2] - 1.5*tau[i][j][0][0][2];
																		  detau[i][j][k][1] += -0.5*tau[i][j][2][1][2] + 2.0*tau[i][j][1][1][2] - 1.5*tau[i][j][0][1][2];
																		  detau[i][j][k][2] += -0.5*tau[i][j][2][2][2] + 2.0*tau[i][j][1][2][2] - 1.5*tau[i][j][0][2][2];

																}

													 }

													 //top of feature
													 if(k==h){

																detau[i][j][k][0] = 0.5 * ( tau[ip][j][k][0][0] - tau[im][j][k][0][0] + tau[i][jp][k][0][1] - tau[i][jm][k][0][1] );
																detau[i][j][k][1] = 0.5 * ( tau[ip][j][k][1][0] - tau[im][j][k][1][0] + tau[i][jp][k][1][1] - tau[i][jm][k][1][1] );
																detau[i][j][k][2] = 0.5 * ( tau[ip][j][k][2][0] - tau[im][j][k][2][0] + tau[i][jp][k][2][1] - tau[i][jm][k][2][1] );
																
																if(cav_flag[i][j][km]>-1){

																		  detau[i][j][k][0] += 0.5 * ( tau[i][j][kp][0][2] - tau[i][j][km][0][2] );
																		  detau[i][j][k][1] += 0.5 * ( tau[i][j][kp][1][2] - tau[i][j][km][1][2] );
																		  detau[i][j][k][2] += 0.5 * ( tau[i][j][kp][2][2] - tau[i][j][km][2][2] );
																}
																else{

																		  detau[i][j][k][0] += -0.5*tau[i][j][k+2][0][2] + 2.0*tau[i][j][k+1][0][2] - 1.5*tau[i][j][k][0][2];
																		  detau[i][j][k][1] += -0.5*tau[i][j][k+2][1][2] + 2.0*tau[i][j][k+1][1][2] - 1.5*tau[i][j][k][1][2];
																		  detau[i][j][k][2] += -0.5*tau[i][j][k+2][2][2] + 2.0*tau[i][j][k+1][2][2] - 1.5*tau[i][j][k][2][2];
																}

													 }


										  }


								}
					 }
		  }


}
