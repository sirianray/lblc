#include "lb.h"
#include <math.h>

//Calculate S_eq
double getS()
{
		  double S=0;
		  S = 0.25+0.75*sqrt(1.0-eight3rd/U);
		  return S;
}

//Kronecker delta
double Delta(int i, int j)
{
		  if(i==j) return 1.0;
		  
		  return 0;
}

double OneThirdDelta(int i, int j)
{

		  if(i==j) return third;
		  else return 0.0;
}

void cal_fequ(double *f, double *u0)
{
		  int i,j,k,ii=0,jj=0,kk=0, index;

		  double E2[3][3], E1[3][3];
		  double trsigma=0, temp;
		  double A=0, A0=0, A1=0, A2=0;
		  double B=0, B0=0, B1=0, B2=0;
		  double C=0, C0=0, C1=0, C2=0;
		  double D=0, D0=0, D1=0, D2=0;
		  double temp0=0;

		  double sigma_mat[3][3] = {0};

		  B2 = one24th;
		  B1 = third;
		  B0 = 0.0;

		  C2 = -one24th;
		  C1 = -one12th;
		  C0 = -two3rd;

		  D2 = 0.0625;
		  D1 = 0.5;
		  D0 = 0.0;

		  for (i=0; i<Nx; i++) {
					 for (j=0; j<Ny; j++) {
								for (k=0; k<Nz; k++) {

										  index = i + j*Nx + k*Nx*Ny;

										  kk=0;
										  for(ii=0;ii<3;ii++){
													 for(jj=ii;jj<3;jj++){
																sigma_mat[ii][jj] = sigma[6*index+kk];
																++kk;
													 }
										  }

										  sigma_mat[1][0] = sigma_mat[0][1];
										  sigma_mat[2][0] = sigma_mat[0][2];
										  sigma_mat[2][1] = sigma_mat[1][2];

										  trsigma = sigma_mat[0][0] + sigma_mat[1][1] + sigma_mat[2][2];	  
										  A2 = -one30th*trsigma;
										  A1 =  A2;
										  A0 =  Rho[index] - 14.0 * A2;

										  for(ii=0;ii<3;ii++){
													 for(jj=0;jj<3;jj++){
																E2[ii][jj] = 0.0625*(-sigma_mat[ii][jj] + trsigma*OneThirdDelta(ii,jj));
																E1[ii][jj] = 8.0*E2[ii][jj];
													 }
										  }

										  //Calculate f_eq[i][j][k] 
										  for(ii=0;ii<15;ii++){
													 temp = 0;
													 if(ii==0){
																A = A0;
																B = B0;
																C = C0;
																D = D0;						
													 }
													 else if(ii<=6) {
																A = A1;
																B = B1;
																C = C1;
																D = D1;
																for(jj=0;jj<3;jj++){
																		  for(kk=0;kk<3;kk++){
																					 temp+=E1[jj][kk]*e[ii][jj]*e[ii][kk];
																		  }
																}
													 }
													 else {
																A = A2;
																B = B2;
																C = C2;
																D = D2;
																for(jj=0;jj<3;jj++){
																		  for(kk=0;kk<3;kk++){
																					 temp+=E2[jj][kk]*e[ii][jj]*e[ii][kk];
																		  }
																}
													 }

													 temp += A/rho;

													 for(jj=0;jj<3;jj++){
																temp += B*u0[3*index+jj]*e[ii][jj]+C*u0[3*index+jj]*u0[3*index+jj];

																for(kk=0;kk<3;kk++){
																		  temp += D*u0[3*index+jj]*u0[3*index+kk]*e[ii][jj]*e[ii][kk];
																}
													 }

													 f[15*index+ii] = temp*rho;
										  }

								}
					 }
		  }

}

void cal_feqf_new()
{
		  int i,j,k,ii=0,jj=0,kk=0, index;
		  int compressibility=0;		// 0: incompressible flow; nonzero: compressible flow

		  double E[3][3], E1, E2, rho0=rho;
		  double trsigma=0, trsigmar=0, temp;
		  double A=0, A0=0, A1=0, A2=0;
		  double B=0, B0=0, B1=0, B2=0;
		  double C=0, C0=0, C1=0, C2=0;
		  double D=0, D0=0, D1=0, D2=0;
		  double udote, usq, udote_sq;
		  double temp1, temp2, rhototal=0;

		  double sigma_mat[3][3] = {0};

		  B2 = one24th;
		  B1 = third;
		  //	B0 = 0.0;

		  C2 = -one24th;
		  C1 = -one12th;
		  C0 = -two3rd;

		  D2 = 0.0625;
		  D1 = 0.5;
		  //	D0 = 0.0;

		  for (k=0; k<Nz; k++) {
					 for (j=0; j<Ny; j++) {
								for (i=0; i<Nx; i++) {
										  
										  index = i + j*Nx + k*Nx*Ny;

										  kk=0;
										  for(ii=0;ii<3;ii++){
													 for(jj=ii;jj<3;jj++){
																sigma_mat[ii][jj] = sigma[6*index+kk];
																++kk;
													 }
										  }

										  sigma_mat[1][0] = sigma_mat[0][1];
										  sigma_mat[2][0] = sigma_mat[0][2];
										  sigma_mat[2][1] = sigma_mat[1][2];

										  trsigma = sigma_mat[0][0] + sigma_mat[1][1] + sigma_mat[2][2];

										  A2 = -one30th*trsigma;
										  A1 = A2;
										  A0 = Rho[index] - 14.0 * A2;

										  usq = 0;
										  for(ii=0;ii<3;ii++){

													 for(jj=0;jj<3;jj++){
																E[ii][jj] = 0.0625*(-sigma_mat[ii][jj] + trsigma*OneThirdDelta(ii,jj) );
													 }

													 usq += u[3*index+ii]*u[3*index+ii];
										  }

										  //Calculate f_eq[i][j][k] 
										  for(ii=0;ii<15;ii++){

													 udote = 0;
													 E2    = 0;

													 for(jj=0; jj<3; jj++){

																udote += u[3*index+jj]*e[ii][jj];

																for(kk=0; kk<3; kk++){
																		  E2 += E[jj][kk] * e[ii][jj] * e[ii][kk];
																}
													 }

													 udote_sq = udote * udote;
													 E1       =   8.0 * E2;

													 if(compressibility!=0){
																rho0 = Rho[index];
													 }
													 if(ii==0){
																f_eq[15*index+ii] = A0 + rho0 * (           C0*usq                   );
													 } else if(ii<7) {
																f_eq[15*index+ii] = A1 + rho0 * (B1*udote + C1*usq + D1*udote_sq + E1);
													 } else {
																f_eq[15*index+ii] = A2 + rho0 * (B2*udote + C2*usq + D2*udote_sq + E2);
													 }
													 rhototal += f_eq[15*index+ii];
										  }
								}
					 }
		  }
}

void cal_p(double *f)
{
		  int i,j,k,ii=0,jj=0,kk=0, index;
		  double temp, td, tu, dx, dy, dz, T1=0, T2=0;

		  T2 = 1.0/24.0;
		  T1 = 8*T2;

		  for (k=0; k<Nz; k++) {
					 for (j=0; j<Ny; j++) {
								for (i=0; i<Nx; i++) {

										  index = i + j*Nx + k*Nx*Ny;

										  for (ii=0; ii<15; ii++) {

													 if(yforce>0){
																temp = yforce * e[ii][1];

																if(ii>0 && ii<7){
																		  p[15*index+ii] = T1 * temp;
																}else if(ii>6){
																		  p[15*index+ii] = T2 * temp;
																}else {
																		  p[15*index+ii] = 0 * temp;
																}
													 }
													 else p[15*index+ii] = 0;


													 // This is the portion of stress that go into forcing term
													 if(ii>0 && ii<7) p[15*index+ii] += T1 * (dtau[3*index+0]*e[ii][0] +dtau[3*index+1]*e[ii][1] +dtau[3*index+2]*e[ii][2] ); 
													 else if (ii>6) p[15*index+ii] += T2 * (dtau[3*index+0]*e[ii][0] +dtau[3*index+1]*e[ii][1] +dtau[3*index+2]*e[ii][2] );
										  }
								}
					 }
		  }
}


void cal_W()
{
		  int i, j, k, index;
		  int im=0, jm=0, km=0;
		  int ip=0, jp=0, kp=0;

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

										  W[9*index+0] =  0.5*(u[3*ip+0] - u[3*im+0]);
										  W[9*index+3] =  0.5*(u[3*ip+1] - u[3*im+1]);
										  W[9*index+6] =  0.5*(u[3*ip+2] - u[3*im+2]);

										  W[9*index+1] =  0.5*(u[3*jp+0] - u[3*jm+0]);
										  W[9*index+4] =  0.5*(u[3*jp+1] - u[3*jm+1]);
										  W[9*index+7] =  0.5*(u[3*jp+2] - u[3*jm+2]);

										  if(k==0){

													 W[9*index+2] = third * u[3*kp+0] + u[3*index+0];
													 W[9*index+5] = third * u[3*kp+1] + u[3*index+1] - four3rd * uy_bottom;
													 W[9*index+8] = third * u[3*kp+2] + u[3*index+2];

										  } else if(k==Nz-1) {

													 W[9*index+2] =-third * u[3*km+0] - u[3*index+0];
													 W[9*index+5] =-third * u[3*km+1] - u[3*index+1]+ four3rd * uy_top;
													 W[9*index+8] =-third * u[3*km+2] - u[3*index+2];

										  } else {

													 W[9*index+2] =  0.5 * (u[3*kp+0] - u[3*km+0]);
													 W[9*index+5] =  0.5 * (u[3*kp+1] - u[3*km+1]);
													 W[9*index+8] =  0.5 * (u[3*kp+2] - u[3*km+2]);
										  }
								}
					 }
		  }

}


void cal_rho(double *f)
{
		  int i,j;
		  double points, temp=0;

		  points = Nx*Ny*Nz;

		  for(i=0;i<points;i++){

					 if(cav_flag[i]>-1){

								temp = 0;
								for(j=0;j<15;j++) temp += f[15*i+j];
								Rho[i] = temp;

					 } else {
								Rho[i] = 1;
					 }
					 
		  }

}

void cal_u(double *f)
{
		  int i,j,k;
		  int points = Nx*Ny*Nz;
		  double temp;

		  cal_rho(f);

		  for(i=0;i<points;i++){

					 if(cav_flag[i]>-1) {

								for(j=0;j<3;j++){ 

										  temp=0;
										  for(k=0;k<15;k++){
													 temp += f[15*i+k]*(double)e[k][j];
										  }

										  
										  u[3*i+j] = temp/Rho[i];
								}
					 }
					 else {
								for(j=0;j<3;j++) u[3*i+j] = 0;
					 }

		  }

}

void cal_sigma()
{
		  int i,j,points = Nx*Ny*Nz;

		  for (i=0; i<points; i++) {

					 for(j=0;j<6;j++){
								if(j==0 || j==3 || j==5) sigma[6*i+j] = -Rho[i] * T + sigma_q[6*i+j];
								else                     sigma[6*i+j] =  sigma_q[6*i+j];
					 }

		  }

		  
		  
}


void monitor(int flag, int print)
{
		  int i, j, k, ii, index;
		  double utotal[3]={0}, rhototal=0.0, rhotemp, utemp, vtemp, wtemp;
		  double q_diff_sq=0, k_eng_new=0, k_diff;

		  double temp1 = 0;

		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Nx + k*Nx*Ny;

										  if(!cav_on || (cav_on && cav_flag[index]>-1) ){

													 for(ii=0;ii<5;ii++){

																q_diff_sq += (Q1[5*index+ii]-Q2[5*index+ii])*(Q1[5*index+ii]-Q2[5*index+ii]);

																if( (q_diff_sq != q_diff_sq) || !isfinite(q_diff_sq) ){

																		  if(q_diff_sq != q_diff_sq) printf("NAN detected in loop\n");
																		  if(!isfinite(q_diff_sq)) printf("INF detected in loop\n");

																		  printf("ii i j k = %d %d %d %d\n",ii,i,j,k);
																		  printf("Q1 = %lf Q2 = %lf\n",Q1[5*index+ii], Q2[5*index+ii]);
																		  exit(1);
																}
													 }

													 for(ii=0;ii<3;ii++) {
																k_eng_new+=0.5*(u[3*index+ii]*u[3*index+ii]+u[3*index+ii]*u[3*index+ii]+u[3*index+ii]*u[3*index+ii])*Rho[index];
																utotal[ii]+=u[3*index+ii];
													 }

													 rhototal+=Rho[index];

										  }
								}
					 }
		  }

		  if(print!=0 && u_on) printf("\nvel=%20.15f %20.15f %20.15f\n",utotal[0], utotal[1], utotal[2]);
		  if(print!=0 && u_on) printf("rho=%20.15f\n",rhototal);

		  if(flag<=0){

					 if(print!=0) printf("dif_Q^2=%30.26f\n",q_diff_sq);

					 if(q_diff_sq != q_diff_sq){
								printf("Nan detected in Q-field\n");
							  	exit(1);
					 }
		  }

		  if(q_diff_sq < 1e-20 && flag<=0) Q_converge=1;
		  
		  if(flag>=0){
					 if(k_eng<1e-15 && k_eng>-1e-15){
								k_diff = k_eng_new;
					 }else{
								k_diff = (k_eng_new-k_eng)/k_eng;
					 }

					 if(((k_eng-k_eng_new<1e-18 && k_eng-k_eng_new>-1e-18)||(k_diff<1e-13 && k_diff>-1e-13)) && flag>=0) u_converge=1;
					 if(print!=0 && flag>=0)printf("diff_k =%30.25f\n",k_diff);
		  }

		  k_eng=k_eng_new;
}

int bounce(int i)
{	
		  switch( i ) 
		  {
					 case 0:
								return 0;
								break;
					 case 1 :
								return 3;
								break;
					 case 2 :
								return 4;
								break;
					 case 3 :
								return 1;
								break;
					 case 4 :
								return 2;
								break;
					 case 5 :
								return 6;
								break;
					 case 6 :
								return 5;
								break;
					 case 7 :
								return 13;
								break;
					 case 8 :
								return 14;
								break;
					 case 9 :
								return 11;
								break;
					 case 10 :
								return 12;
								break;
					 case 11 :
								return 9;
								break;
					 case 12 :
								return 10;
								break;
					 case 13 :
								return 7;
								break;
					 case 14 :
								return 8;
								break;
					 default:
								printf("error: index of bounce() is out of range");
								return -1;
		  }
}


double QQ(double M[3][3], int i, int j)
{
		  return M[i][0]*M[0][j]+M[i][1]*M[1][j]+M[i][2]*M[2][j];
}

void write_restart(double *Q, double *Rho, double *u)
{
		  FILE* frestart;
		  int i, j, k, ii, jj;
		  int points = Nx*Ny*Nz;

		  frestart=fopen("restart.out","wb");

		  for(i=0; i<points; i++){
					 for(ii=0;ii<5;ii++) fwrite(&Q[5*i+ii],sizeof(double),1,frestart);
		  }

		  for(i=0; i<points; i++){
					 fwrite(&Rho[i],sizeof(double),1,frestart);
		  }

		  for(i=0; i<points; i++){
					 for(ii=0; ii<3; ii++) fwrite(&u[5*i+ii],sizeof(double),1,frestart);
		  }

		  fclose(frestart);
}

void read_restart(double *Q, double *Rho, double *u)
{
		  FILE* frestart;
		  int points = Nx*Ny*Nz;

		  int i,ii;
		  frestart=fopen("restart.out","rb");


		  for(i=0; i<points; i++){
					 for(ii=0; ii<5; ii++) fread(&Q[5*i+ii], sizeof(double),1,frestart);
		  }

		  for(i=0; i<points; i++){
					 fread(&Rho[i], sizeof(double),1,frestart);
		  }

		  for(i=0; i<points; i++){
					 for(ii=0; ii<3; ii++) fread(&u[3*i+ii], sizeof(double),1,frestart);
		  }


		  fclose(frestart);
}


void print_visual(int input)
{
		  int index, i,j,k, points=Nx*Ny*Nz;
		  FILE *qfile, *ufile, *grid, *sfile, *dfile, *ffile;

		  if(input){
					 qfile=fopen("Q.out","w");
					 if(u_on) ufile=fopen("vel.out","w");
					 grid =fopen("grid.out","w");
					 sfile=fopen("stress.out","w");
					 dfile=fopen("den.out","w");
					 ffile=fopen("f.out","w");
		  }
		  else{
					 qfile=fopen("Q.out","a");
					 if(u_on) ufile=fopen("vel.out","a");
					 sfile=fopen("stress.out","a");
					 dfile=fopen("den.out","a");

					 grid =fopen("grid.out","a");
		  }

		  if(input) fprintf(grid,"Nx Ny Nz %d %d %d\n",Nx,Ny,Nz);

		  for(k=0;k<Nz;k++){
					 for(j=0;j<Ny;j++){
								for(i=0;i<Nx;i++){

										  index = i + j*Ny + k*Nx*Ny;

										  fprintf(qfile,"%e %e %e %e %e\n",Q1[5*index+0],Q1[5*index+1],Q1[5*index+2],Q1[5*index+3],Q1[5*index+4]);
										  fprintf(sfile,"%e %e %e\n",dtau[3*index+0],dtau[3*index+1],dtau[3*index+2]);
										  fprintf(dfile,"%e\n",Rho[index]);

										  int ii;
										  for(ii=0;ii<15;ii++) fprintf(ffile,"%lf\t",f2[15*index+ii]);
										  fprintf(ffile,"\n");

										  if(u_on) fprintf(ufile,"%e %e %e\n",u[3*index+0],u[3*index+1],u[3*index+2]);
										  if(input) fprintf(grid,"%d %d %d %d\n",i,j,k,cav_flag[index]);
								}
					 }
		  }

		  fclose(qfile);
		  fclose(sfile);
		  fclose(dfile);
		  fclose(ffile);

		  if(u_on) fclose(ufile);
		  if(input) fclose(grid);
}

void cav_surf(){

		  int i, j, k, h, w, start, stop, index;

		  h = cav_height;
		  w = cav_width;

		  //Beggning and end of cavity
		  start = floor(0.5*(Ny-w));
		  stop  = floor(0.5*(Ny+w));

		  //Surface is inbetween points with flag=-1 and flag=+1
		  if(cav_on){
					 for(k=0;k<Nz;k++){
								for(j=0;j<Ny;j++){
										  for(i=0;i<Nx;i++){
													 index = i + Ny*j + Nx*Ny*k;
													 if( (j==start && k<=h) || (j==stop && k<=h) || (j>=start && j<=stop && k==h) ) cav_flag[index] = 1;
													 if( j>start && j<stop && k<h) cav_flag[index] = -1;
										  }
								}
					 }
		  }
}
