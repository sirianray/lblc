#include "lb.h"

void evol_Q(double *Qin, double *Qout){

		  int i, j, k, ii, jj, kk, index;
		  double D[3][3], Omega[3][3], S[3][3], trQW, temp;
		  double Wmat[3][3], Qmat[3][3], Svec[5] = {0};

		  cal_dQ(Qin);

		  //Bulk evolution
		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Nx + k*Nx*Ny;

										  if( (cav_flag[index]>=0 && cav_on) || !cav_on) {

													 kk=0;
													 for(ii=0;ii<2;ii++){
																for(jj=ii;jj<3;jj++){
																		  Qmat[ii][jj] = Qin[5*index+kk];
																		  ++kk;
																}
													 }

													 Qmat[1][0] = Qmat[0][1];
													 Qmat[2][0] = Qmat[0][2];
													 Qmat[2][1] = Qmat[1][2];
													 Qmat[2][2] = -Qmat[0][0]-Qmat[1][1];

													 if(u_on){

																trQW=0.0;


																kk=0;
																for(ii=0;ii<3;ii++){
																		  for(jj=0;jj<3;jj++){
																					 Wmat[ii][jj] = W[9*index+kk];
																					 ++kk;
																		  }
																}

																for(ii=0;ii<3;ii++){
																		  for(jj=0;jj<3;jj++){
																					 D[ii][jj] = 0.5*(Wmat[ii][jj] + Wmat[jj][ii]);
																					 Omega[ii][jj] = 0.5*(Wmat[ii][jj] - Wmat[jj][ii]);
																					 trQW = trQW + Qmat[ii][jj]*Wmat[jj][ii];
																		  }
																}

																for(ii=0;ii<3;ii++){
																		  for(jj=0;jj<3;jj++){

																					 temp = 0;

																					 for(kk=0;kk<3;kk++){
																								temp += (xi*D[ii][kk] + Omega[ii][kk])*(Qmat[kk][jj] + OneThirdDelta(kk,jj)) + (Qmat[ii][kk] + OneThirdDelta(ii,kk))*(xi*D[kk][jj] - Omega[kk][jj]);
																					 }

																					 temp += - 2.0*xi*(Qmat[ii][jj] + OneThirdDelta(ii,jj) )*trQW;
																					 S[ii][jj] = temp;
																		  }
																}

																kk=0;
																for(ii=0;ii<2;ii++){
																		  for(jj=ii;jj<3;jj++){
																					 Svec[kk] = S[ii][jj];
																					 ++kk;
																		  }
																}

													 }


													 for(ii=0;ii<5;ii++) {
																//if(H[5*index+ii]>1.0) printf("PROBLEM IN H-FIELD!\n");
																Qout[5*index+ii] = Qin[5*index+ii] + 0.5*(Gamma_rot * H[5*index+ii] + Svec[ii] - convQ[5*index+ii]);
													 }
										  }
								}
					 }
		  }

		  //Evoluation for top and bottom surface
		  if(!wall_inf){

					 surface_derivative(Qin);

					 if(wall_degen) degenerate_anchoring();
					 else {
								for(j=0;j<Ny;j++){
										  for(i=0;i<Nx;i++){

													 index = i + j*Ny;

													 //Using outward normals
													 for(ii=0;ii<5;ii++){
																Qsurfbot[5*index+ii] += - dt*Gamma_rot * (-1 * kappa * dQdxnu_bot[5*index+ii] + W_wall*( Qsurfbot[5*index+ii] - Qo_bottom[ii]) );
													 }

													 for(ii=0;ii<5;ii++){
																Qsurftop[5*index+ii] += - dt*Gamma_rot * (kappa * dQdxnu_top[5*index+ii] + W_wall*( Qsurftop[5*index+ii] - Qo_top[ii]));
													 }
										  }
								}
					 }
		  }
}

void degenerate_anchoring(){

		  double Pmat[3][3], Qtilde[3][3], Qtildevec[5], Qperp[3][3], Qperpvec[5], Q[3][3], Qmatbot[3][3], Qmattop[3][3], nu[3], nu_bot[3] = {0, 0, -1}, nu_top[3] = {0, 0, +1}, nuQnu, temp;
		  int surf_points = 2*Nx*Ny, ii,jj,nn,i,k,kk;

		  for(k=0;k<2;k++){

					 for(ii=0;ii<3;ii++){
								if(k==0) nu[ii] = nu_bot[ii]; 
								if(k==1) nu[ii] = nu_top[ii]; 

								for(i=0;i<surf_points;i++){

										  for(ii=0;ii<2;ii++){
													 for(jj=ii;jj<3;jj++){
																if(k==0) Q[ii][jj]=Qsurfbot[5*i+k];
																if(k==1) Q[ii][jj]=Qsurftop[5*i+k];
																++k;
													 }
										  }

										  Q[1][0] = Q[0][1];
										  Q[2][0] = Q[0][2];
										  Q[2][1] = Q[1][2];
										  Q[2][2] = -Q[0][0]-Q[1][1];

										  nuQnu = 0;
										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																Qtilde[ii][jj] = Q[ii][jj] + OneThirdDelta(ii,jj)*S;
																Pmat[ii][jj] = Delta(ii,jj) - nu[ii]*nu[jj];
																nuQnu += nu[ii]*Qtilde[ii][jj]*nu[jj];
													 }
										  }

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																for(kk=0;kk<3;kk++){
																		  for(nn=0;nn<3;nn++){
																					 Qperp[ii][jj] += Pmat[ii][kk]*Q[kk][nn]*Pmat[nn][jj];
																		  }
																}
													 }
										  }

										  for(ii=0;ii<2;ii++){
													 for(jj=ii;jj<3;jj++){
																Qtildevec[kk] = Qtilde[ii][jj]; 
																Qperpvec[kk] = Qperp[ii][jj]; 
																++kk;
													 }
										  }

										  if(k==0){
													 for(ii=0;ii<5;ii++){
																if(ii==0 || ii==3) temp = nuQnu;
																else temp = 0;
																Qsurfbot[5*i+ii] = Qsurfbot[5*i+ii] - dt*Gamma_rot * (-1 * kappa * dQdxnu_bot[5*i+ii] + 2*W_wall*( (Qtildevec[ii]-Qperpvec[ii]) - temp));
													 }
										  }

										  if(k==1){
													 for(ii=0;ii<5;ii++){
																if(ii==0 || ii==3) temp = nuQnu;
																else temp = 0;
																Qsurftop[5*i+ii] = Qsurftop[5*i+ii] - dt*Gamma_rot * (     kappa * dQdxnu_top[5*i+ii] + 2*W_wall*( (Qtildevec[ii] - Qperpvec[ii]) - temp));
													 }
										  }
								}
					 }
		  }
}

void degenerate_cavity(){
}

void surface_derivative(double *Q0){

		  //bottom surface
		  int i, j, k, ii, index;

		  double Qvec1[5]={0}, Qvec0[5]={0};
		  double QvecNz1[5]={0}, QvecNz2[5]={0};

		  //Bottom Surface
		  //dx1 = 0.5, dx2 = 1.0, idx = dx1*dx2*(dx1+dx2);
		  for (j=0; j<Ny; j++) {
					 for (i=0; i<Nx; i++) {

								index = i + j*Nx;
								for(ii=0;ii<5;ii++) Qvec0[ii] = Q0[5*index+ii];

								k = 1;
								index = i + j*Nx + k*Nx*Ny;
								for(ii=0;ii<5;ii++) Qvec1[ii] = Q0[5*index+ii];

								index = i + j*Nx;
								for(ii=0; ii<5; ii++)  dQdxnu_bot[5*index+ii] = -eight3rd*Qsurfbot[5*index+ii] + 3.0*Qvec0[ii] - third*Qvec1[ii];
					 }
		  }

		  //Top Surface
		  //dx1 = 1.0, dx2 = 0.5, idx = dx1*dx2*(dx1+dx2);
		  for (j=0; j<Ny; j++) {
					 for (i=0; i<Nx; i++) {

								k = Nz-1;
								index = i + j*Nx + k*Nx*Ny;
								for(ii=0;ii<5;ii++) QvecNz1[ii] = Q0[5*index+ii];

								k = Nz-2;
								index = i + j*Nx + k*Nx*Ny;
								for(ii=0;ii<5;ii++) QvecNz2[ii] = Q0[5*index+ii];

								index = i + j*Nx;
								for (ii=0; ii<5; ii++)  dQdxnu_top[5*index+ii] = eight3rd*Qsurftop[5*index+ii] - 3.0*QvecNz1[ii] + third*QvecNz2[ii];
					 }
		  }

}

void cal_dQ(double *Q0)
{
		  int i, j, k, ii=0, jj=0;
		  int ip, im, jp, jm, kp, km;

		  double temp=0, trQ2=0, trh3rd=0, d2Q[3][3]={0}, uvec[3]={0};

		  int halfNy = 0.5*Ny, h = cav_height;
		  double Hmat[3][3]={0}, Qmat[3][3]={0}, Qijk[5] = {0}, d2[15]={0}, ddQ[15]={0};

		  int index, surf_index;
		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Nx + k*Nx*Ny;
										  surf_index = i + j*Nx;

										  derivative(i,j,k,Q0,d2,ddQ);

										  //Q at current point
										  for(ii=0;ii<5;ii++) Qijk[ii] = Q0[5*index+ii];

										  //lapacian Q

										  int n, m;
										  ii=0;
										  for(n=0;n<2;n++){
													 for(m=n;m<3;m++){
																d2Q[n][m] = d2[5*0+ii]+d2[5*1+ii]+d2[5*2+ii];
																++ii;
													 }
										  }

										  d2Q[1][0] = d2Q[0][1];
										  d2Q[2][0] = d2Q[0][2];
										  d2Q[2][1] = d2Q[1][2];
										  d2Q[2][2] =-d2Q[0][0] - d2Q[1][1];

										  // calculate H				
										  trQ2=0;
										  for(ii=0;ii<5;ii++){
													 trQ2 += 2.0*(Qijk[ii]*Qijk[ii]);
										  }
										  trQ2 += 2.0*Qijk[0]*Qijk[3];

										  //Populate Q-matrix 
										  ii=0;
										  for(n=0;n<2;n++){
													 for(m=n;m<3;m++){
																Qmat[n][m] = Q0[5*index+ii];
																++ii;
													 }
										  }
										  Qmat[1][0] = Qmat[0][1];
										  Qmat[2][0] = Qmat[0][2];
										  Qmat[2][1] = Qmat[1][2];
										  Qmat[2][2] = -Qmat[0][0]-Qmat[1][1];

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																temp = Qmat[ii][jj];
																Hmat[ii][jj] = -1*(A_ldg*(1.0-third*U)*temp - A_ldg*U*( QQ(Qmat,ii,jj) - trQ2*(temp + OneThirdDelta(ii,jj))) - kappa * d2Q[ii][jj]);
													 }
										  }

										  ii=0;
										  for(n=0;n<2;n++){
													 for(m=n;m<3;m++){
																H[5*index+ii] = Hmat[n][m];
																++ii;
													 }
										  }


										  //	calculate convective Q
										  if(u_on){
													 for(ii=0;ii<3;ii++) uvec[ii] = u[3*index+ii];
													 for(ii=0;ii<5;ii++)  convQ[5*index+ii] = uvec[0]*ddQ[5*0+ii] + uvec[1]*ddQ[5*1+ii] + uvec[2]*ddQ[5*2+ii];
										  }
								}
					 }
		  }
}

// calculate Q tensor related sigma
void cal_stress(double *Q0)
{

		  int i, j, k, ii=0, jj=0, kk=0, mm, index, surf_index;
		  int ip, im, jp, jm, kp, km, m,n;
		  double temp, temp2;
		  double trQ2, trH3rd, d2Q[3][3], dQ2;
		  double sum;

		  double Hmat[3][3]={0},Qmat[3][3]={0}, Qijk[5] = {0}, d2[15]={0}, ddQvec[15]={0}, ddQ[3][5]={0};
		  double sigma_mat[3][3]={0}, tau_mat[3][3]={0};

		  int halfNy = 0.5*Ny;
		  int h = cav_height;

		  for (k=0;k<Nz;k++) {
					 for (j=0;j<Ny;j++) {
								for (i=0;i<Nx;i++) {

										  index = i + j*Nx + k*Nx*Ny;
										  surf_index = i + j*Nx;

										  for(ii=0;ii<5;ii++) Qijk[ii] = Q0[5*index+ii];

										  derivative(i,j,k,Q0,d2,ddQvec);

										  for(m=0;m<3;m++){
													 for(n=0;n<5;n++){
																ddQ[m][n] = ddQvec[5*m+n];
													 }
										  }

										  double dQ[3][3]={0};
										  for(m=0;m<2;m++){
													 for(n=m;n<2;n++){

																for(ii=0;ii<5;ii++) dQ[m][n] += 2*(ddQ[m][ii]*ddQ[n][ii]);

																dQ[m][n] += ddQ[m][0]*ddQ[n][3]+ddQ[m][3]*ddQ[n][0];

													 }
										  }

										  dQ[1][0]=dQ[0][1];
										  dQ[2][0]=dQ[0][2];
										  dQ[2][1]=dQ[1][2];

										  dQ2=0;
										  for(ii=0;ii<3;ii++) dQ2 += dQ[ii][ii];

										  ii=0;
										  for(n=0;n<2;n++){
													 for(m=n;m<3;m++){
																d2Q[n][m] = d2[5*0+ii]+d2[5*1+ii]+d2[5*2+ii];
																++ii;
													 }
										  }

										  d2Q[1][0] = d2Q[0][1];
										  d2Q[2][0] = d2Q[0][2];
										  d2Q[2][1] = d2Q[1][2];
										  d2Q[2][2] =-d2Q[0][0] - d2Q[1][1];

										  //	calculate H				
										  trQ2=0;
										  for(ii=0;ii<5;ii++){
													 trQ2 += 2.0*(Qijk[ii]*Qijk[ii]);
										  }
										  trQ2 += 2.0*Qijk[0]*Qijk[3];

										  kk=0;
										  for(ii=0;ii<2;ii++){
													 for(jj=ii;jj<3;jj++){
																Qmat[ii][jj] = Qijk[kk];
																++kk;
													 }
										  }
										  Qmat[1][0] = Qmat[0][1];
										  Qmat[2][0] = Qmat[0][2];
										  Qmat[2][1] = Qmat[1][2];
										  Qmat[2][2] = -Qmat[0][0]-Qmat[1][1];

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																temp = Qmat[ii][jj];
																Hmat[ii][jj] = -1*(A_ldg*(1.0-third*U)*temp - A_ldg*U*( QQ(Qmat,ii,jj) - trQ2*(temp + OneThirdDelta(ii,jj))) - kappa * d2Q[ii][jj]);
													 }
										  }

										  sum = 0;
										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																sum = sum + Qmat[ii][jj]*Hmat[ii][jj];
													 }
										  }

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){

																//Populate Sigma (Q-dependent terms only!)
																temp = 0.5*kappa*dQ2*Delta(ii,jj) + 2*xi*(Qmat[ii][jj] + OneThirdDelta(ii,jj))*sum;
																for(mm=0;mm<3;mm++){
																		  temp += - xi*Hmat[ii][mm]*(Qmat[mm][jj] + OneThirdDelta(mm,jj)) - xi*(Qmat[ii][mm] + OneThirdDelta(ii,mm))*Hmat[mm][jj] ;
																}
																temp += -kappa * dQ[ii][jj];
																sigma_mat[ii][jj] = temp;

																//populate tau
																temp2= 0;
																for(mm=0;mm<3;mm++){
																		  temp2 += Qmat[ii][mm]*Hmat[mm][jj] - Hmat[ii][mm]*Qmat[mm][jj];
																}
																tau_mat[ii][jj] = temp2;
													 }
										  }
										  

										  //Transform matrices to vectors 
										  kk=0;
										  for(ii=0;ii<2;ii++){
													 for(jj=ii;jj<3;jj++){
																H[5*index+kk] = Hmat[ii][jj];
																++kk;
													 }
										  }


										//  kk=0;
										//  for(ii=0;ii<3;ii++){
										//			 for(jj=ii;jj<3;jj++){
										//						sigma_q[6*index+kk] =  sigma_mat[ii][jj];
										//						++kk;
										//			 }
										//  }

										  kk=0;
										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																tau[9*index+kk] = tau_mat[ii][jj] + sigma_mat[ii][jj];
																++kk;
													 }
										  }
										  //printf("\n");
								}
					 }
		  }
		  //exit(1);
}

void cal_dtau(double *Q0)
{
		  int i, j, k, ii=0, jj=0, kk=0, index;
		  int ip, im, jp, jm, kp, km;

		  int halfNy = 0.5*Ny, jm2;
		  int h = cav_height;
		  int index2, index1, index0, indexNz3, indexNz2, indexNz1, jp2, kp2;

		  for (k=0; k<Nz; k++) {
					 for (j=0; j<Ny; j++) {
								for (i=0; i<Nx; i++) {

										  index = i + j*Nx + k*Nx*Ny;

										  im = (i-1) + j*Nx + k*Nx*Ny;
										  jm = i + (j-1)*Nx + k*Nx*Ny;
										  km = i + j*Nx + (k-1)*Nx*Ny;

										  ip = (i+1) + j*Nx + k*Nx*Ny;
										  jp = i + (j+1)*Nx + k*Nx*Ny;
										  kp = i + j*Nx + (k+1)*Nx*Ny;

										  //pbc
										  if(i==0) im = (Nx-1) + j*Nx + k*Nx*Ny;
										  if(j==0) jm = i + (Ny-1)*Nx + k*Nx*Ny;

										  if(i==Nx-1) ip = (0) + j*Nx + k*Nx*Ny;
										  if(j==Ny-1) jp = i + (0)*Nx + k*Nx*Ny;

										  if( (cav_on && cav_flag[index]==0) || !cav_on){

													 //[0 - 00   1 - 01   2 - 02 
													 // 3 - 10   4 - 11   5 - 12
													 // 6 - 20   7 - 21   8 - 22]

													 //x- and y-components
													 dtau[3*index+0] = 0.5 * ( tau[9*ip+0] - tau[9*im+0] + tau[9*jp+1] - tau[9*jm+1] );
													 dtau[3*index+1] = 0.5 * ( tau[9*ip+3] - tau[9*im+3] + tau[9*jp+4] - tau[9*jm+4] );
													 dtau[3*index+2] = 0.5 * ( tau[9*ip+6] - tau[9*im+6] + tau[9*jp+7] - tau[9*jm+7] );

													 if (k==0) {

																index2 = i + j*Nx + 2*Nx*Ny;
																index1 = i + j*Nx + 1*Nx*Ny;
																index0 = i + j*Nx + 0*Nx*Ny;

																dtau[3*index+0] += -0.5*tau[9*index2+2] + 2.0*tau[9*index1+2] - 1.5*tau[9*index0+2];
																dtau[3*index+1] += -0.5*tau[9*index2+5] + 2.0*tau[9*index1+5] - 1.5*tau[9*index0+5];
																dtau[3*index+2] += -0.5*tau[9*index2+8] + 2.0*tau[9*index1+8] - 1.5*tau[9*index0+8];

													 }else if (k==Nz-1) {

																indexNz3 = i + j*Nx + (Nz-3)*Nx*Ny;
																indexNz2 = i + j*Nx + (Nz-2)*Nx*Ny;
																indexNz1 = i + j*Nx + (Nz-1)*Nx*Ny;

																dtau[3*index+0] += 0.5*tau[9*indexNz3+2] - 2.0*tau[9*indexNz2+2] + 1.5*tau[9*indexNz1+2];
																dtau[3*index+1] += 0.5*tau[9*indexNz3+5] - 2.0*tau[9*indexNz2+5] + 1.5*tau[9*indexNz1+5];
																dtau[3*index+2] += 0.5*tau[9*indexNz3+8] - 2.0*tau[9*indexNz2+8] + 1.5*tau[9*indexNz1+8];

													 }else{

																dtau[3*index+0] += 0.5 * ( tau[9*kp+2] - tau[9*km+2] );
																dtau[3*index+1] += 0.5 * ( tau[9*kp+5] - tau[9*km+5] );
																dtau[3*index+2] += 0.5 * ( tau[9*kp+8] - tau[9*km+8] );

													 }
										  }

										  if( cav_on && cav_flag[index]) {


										   		 //left side of feature
										   		 if(j<halfNy && k<h){

										   					jm2 = i + (j-2)*Nx + k*Nx*Ny;

										   					//x-comp
										   					dtau[3*index+0] = 0.5 * ( tau[9*ip+0] - tau[9*im+0] );
										   					dtau[3*index+0] += 0.5*tau[9*jm2+1] - 2*tau[9*jm+1]+1.5*tau[9*index+1];

										   					//y-comp
										   					dtau[3*index+1] = 0.5 * ( tau[3*ip+3] - tau[3*im+3]);
										   					dtau[3*index+1] += 0.5*tau[9*jm2+4] - 2*tau[9*jm+4]+1.5*tau[9*index+4];

										   					//z-comp
										   					dtau[3*index+2]  =  0.5 * ( tau[3*ip+6] - tau[3*im+6] );
										   					dtau[3*index+2] += 0.5*tau[9*jm2+7] - 2*tau[9*jm+7]+1.5*tau[9*index+7];

										   					//modifications for z derivative 
										   					if(k>0){

										   							  dtau[3*index+0] += 0.5 * ( tau[9*kp+2] - tau[9*km+2] );
										   							  dtau[3*index+1] += 0.5 * ( tau[9*kp+5] - tau[9*km+5] );
										   							  dtau[3*index+2] += 0.5 * ( tau[9*kp+8] - tau[9*km+8] );
										   					} 
										   					else{
										   							  index2 = i + j*Nx + (2)*Nx*Ny;
										   							  index1 = i + j*Nx + (1)*Nx*Ny;
										   							  index0 = i + j*Nx + (0)*Nx*Ny;

										   							  dtau[3*index+0] += -0.5*tau[9*index2+2] + 2.0*tau[9*index1+2] - 1.5*tau[9*index0+2];
										   							  dtau[3*index+1] += -0.5*tau[9*index2+5] + 2.0*tau[9*index1+5] - 1.5*tau[9*index0+5];
										   							  dtau[3*index+2] += -0.5*tau[9*index2+8] + 2.0*tau[9*index1+8] - 1.5*tau[9*index0+8];

										   					}

										   		 }

										   		 //right side of feature
										   		 if(j>halfNy && k<h){

										   					jp2 = i + (j+2)*Nx + k*Nx*Ny;

										   					//x-comp
										   					dtau[3*index+0]  = 0.5 * ( tau[9*ip+0] - tau[9*im+0] );
										   					dtau[3*index+0] += -0.5*tau[9*jp2+1] + 2*tau[9*jp+1] - 1.5*tau[9*index+1];

										   					//y-comp
										   					dtau[3*index+1]  = 0.5 * ( tau[9*ip+3] - tau[9*im+3] );
										   					dtau[3*index+1] += -0.5*tau[9*jp2+4] + 2*tau[9*jp+4] - 1.5*tau[9*index+4];

										   					//z-comp
										   					dtau[3*index+2]  = 0.5 * ( tau[9*ip+6] - tau[9*im+6] );
										   					dtau[3*index+2] += -0.5*tau[9*jp2+7] + 2*tau[9*jp+7] - 1.5*tau[9*index+7];

										   					//modifications for z derivative 
										   					if(k>0){

										   							  dtau[3*index+0] += 0.5 * ( tau[9*kp+2] - tau[9*km+2] );
										   							  dtau[3*index+1] += 0.5 * ( tau[9*kp+5] - tau[9*km+5] );
										   							  dtau[3*index+2] += 0.5 * ( tau[9*kp+8] - tau[9*km+8] );
										   					} 
										   					else{
										   							  index2 = i + j*Nx + (2)*Nx*Ny;
										   							  index1 = i + j*Nx + (1)*Nx*Ny;
										   							  index0 = i + j*Nx + (0)*Nx*Ny;

										   							  dtau[3*index+0] += -0.5*tau[9*index2+2] + 2.0*tau[9*index1+2] - 1.5*tau[9*index0+2];
										   							  dtau[3*index+1] += -0.5*tau[9*index2+5] + 2.0*tau[9*index1+5] - 1.5*tau[9*index0+5];
										   							  dtau[3*index+2] += -0.5*tau[9*index2+8] + 2.0*tau[9*index1+8] - 1.5*tau[9*index0+8];

										   					}

										   		 }

										   		 //top of feature
										   		 if(k==h){

										   					dtau[3*index+0] = 0.5 * ( tau[9*ip+0] - tau[9*im+0] + tau[9*jp+1] - tau[9*jm+1] );
										   					dtau[3*index+1] = 0.5 * ( tau[9*ip+3] - tau[9*im+3] + tau[9*jp+4] - tau[9*jm+4] );
										   					dtau[3*index+2] = 0.5 * ( tau[9*ip+6] - tau[9*im+6] + tau[9*jp+7] - tau[9*jm+7] );

										   					if(cav_flag[i+j*Nx+(k-1)*Nx*Ny]>-1){

										   							  dtau[3*index+0] += 0.5 * ( tau[9*kp+2] - tau[9*km+2] );
										   							  dtau[3*index+1] += 0.5 * ( tau[9*kp+5] - tau[9*km+5] );
										   							  dtau[3*index+2] += 0.5 * ( tau[9*kp+8] - tau[9*km+8] );
										   					}
										   					else{
										   							  kp2 = i + j*Nx +(k+2)*Nx*Ny;

										   							  dtau[3*index+0] += -0.5*tau[9*kp2+2] + 2.0*tau[9*kp+2] - 1.5*tau[9*index+2];
										   							  dtau[3*index+1] += -0.5*tau[9*kp2+5] + 2.0*tau[9*kp+5] - 1.5*tau[9*index+5];
										   							  dtau[3*index+2] += -0.5*tau[9*kp2+8] + 2.0*tau[9*kp+8] - 1.5*tau[9*index+8];
										   					}
										   		 }
										  }
								}
					 }
		  }
}

void derivative(int i, int j, int k, double *Q0, double *d2, double *ddQ){

		  int ii,jj;
		  int index = i + j*Nx + k*Nx*Ny;
		  int  im = (i-1) + j*Nx + k*Nx*Ny;
		  int  jm = i + (j-1)*Nx + k*Nx*Ny;
		  int  km = i + j*Nx + (k-1)*Nx*Ny;

		  int  ip = (i+1) + j*Nx + k*Nx*Ny;
		  int  jp = i + (j+1)*Nx + k*Nx*Ny;
		  int  kp = i + j*Nx + (k+1)*Nx*Ny;
		  double Qup[3][5], Qdown[3][5];

		  //Q at current point
		  double Qijk[5]={0};
		  for(ii=0;ii<5;ii++) Qijk[ii] = Q0[5*index+ii];

		  //pbc for x and y
		  if(i==0) im = (Nx-1) + j*Nx + k*Nx*Ny;
		  if(j==0) jm = i + (Ny-1)*Nx + k*Nx*Ny;

		  if(i==Nx-1) ip = (0) + j*Nx + k*Nx*Ny;
		  if(j==Ny-1) jp = i + (0)*Nx + k*Nx*Ny;

		  //derivative coefficients
		  double c[3][3], c2[3][3];

		  if(!cav_on || (cav_on && cav_flag[index]>-1) ){

					 //x coefficients will never change
					 for(ii=0;ii<5;ii++){
								Qup[0][ii] = Q0[5*ip+ii];
								Qdown[0][ii] = Q0[5*im+ii];
					 }

					 c[0][0] =  1;
					 c[0][1] = -2;
					 c[0][2] =  1;

					 c2[0][0] = 0.5;
					 c2[0][1] = 0.0;
					 c2[0][2] = -0.5;

					 //z derivative 
					 if(k==Nz-1){

								for(ii=0;ii<5;ii++){
										  if(wall_inf) Qup[2][ii] = Qo_top[ii];
										  else Qup[2][ii] = Qsurftop[5*index+ii];
										  Qdown[2][ii] = Q0[5*km + ii];
								}

								c[2][0] = eight3rd;
								c[2][1] = -4;
								c[2][2] = four3rd;

								c2[2][0] = four3rd;
								c2[2][1] = -1.0;
								c2[2][2] = -third;
					 }
					 else if(k==0 || (cav_on && cav_flag[km]==-1) ){

								for(ii=0;ii<5;ii++){
										  Qup[2][ii] = Q0[5*kp + ii];

										  if(k==0){
													 if(wall_inf) Qdown[2][ii] = Qo_bottom[ii]; 
													 else Qdown[2][ii] = Qsurfbot[5*index+ii]; 
										  }
								}

								c[2][0] = four3rd;
								c[2][1] = -4;
								c[2][2] = eight3rd;

								c2[2][0] = third;
								c2[2][1] = +1.0;
								c2[2][2] = -four3rd;

					 }
					 else{ //Bulk points
								for(ii=0;ii<5;ii++){
										  Qup[2][ii] = Q0[5*kp + ii];
										  Qdown[2][ii] = Q0[5*km + ii];
								}

								c[2][0] =  1;
								c[2][1] = -2;
								c[2][2] =  1;

								c2[2][0] = 0.5;
								c2[2][1] = 0.0;
								c2[2][2] = -0.5;
					 }

					 //y-derivatives
					 if(cav_flag[index]==1){
								//forward
								if(cav_flag[jm]==-1){

										  for(ii=0;ii<5;ii++){

													 Qup[1][ii] = Q0[5*jp + ii];

													 if(cav_on && cav_inf) Qdown[1][ii] = Qsides[ii]; 
													 else Qdown[1][ii] = Qsurfbot[5*index+ii]; 
										  }


										  c[1][0] = four3rd;
										  c[1][1] = -4;
										  c[1][2] = eight3rd;

										  c2[1][0] = third;
										  c2[1][1] = +1;
										  c2[1][2] = -four3rd;

								}
								else{//backward

								}

					 }
					 else{

								for(ii=0;ii<5;ii++){
										  Qup[1][ii] = Q0[5*jp+ii];
										  Qdown[1][ii] = Q0[5*jm+ii];
								}

								// y first derivative coefficients

								c[1][0] =  1;
								c[1][1] = -2;
								c[1][2] =  1;

								//y second derivative coefficients

								c2[1][0] = 0.5;
								c2[1][1] = 0.0;
								c2[1][2] = -0.5;
					 }


					 for(ii=0;ii<3;ii++){
								for(jj=0;jj<5;jj++){

										  d2[5*ii+jj] = c[ii][0]*Qup[ii][jj] + c[ii][1]*Qijk[jj] + c[ii][2]*Qdown[ii][jj];
										  ddQ[5*ii+jj] = c2[ii][0]*Qup[ii][jj] + c2[ii][1]*Qijk[jj] + c2[ii][2]*Qdown[ii][jj];

								}
					 }
		  }
}
