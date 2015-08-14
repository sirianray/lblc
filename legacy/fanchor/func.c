#include "lb.h"
#include "particle.h"

//Calculate S_eq
real getS()
{
	real S=0;
	S = 0.25+0.75*sqrt(1.0-eight3rd/U);
	return S;
}

//Kronecker delta
real Delta(int i, int j)
{
	if(i==j){
		return 1.0;
	}
	return 0.0;
}
real OneThirdDelta(int i, int j)
{
	if(i==j) {
		return third;
	}
	return 0.0;
}

void clear2d(real A[3][3])
{
	int i, j;
	for (i=0; i<3; i++){
		for (j=0; j<3; j++) A[i][j]=0;
	}
}

void cal_W(real u0[Nx][Ny][Nz][3])
{
	int i, j, k;
	int i_minus=0, j_minus=0, k_minus=0;
	int i_plus=0, j_plus=0, k_plus=0;

	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {

				i_minus = i-1;
		  		j_minus = j-1;
		  		k_minus = k-1;

				i_plus = i+1;
		  		j_plus = j+1;
		  		k_plus = k+1;

		  		if(i==0) i_minus = Nx-1;
		  		if(j==0) j_minus = Ny-1;
		  		if(k==0) k_minus = Nz-1;

		  		if(i==Nx-1) i_plus = 0;
		  		if(j==Ny-1) j_plus = 0;
		  		if(k==Nz-1) k_plus = 0;

				W[i][j][k][0][0] =  0.5*(u0[i_plus][j][k][0] - u0[i_minus][j][k][0]);
				W[i][j][k][1][0] =  0.5*(u0[i_plus][j][k][1] - u0[i_minus][j][k][1]);
				W[i][j][k][2][0] =  0.5*(u0[i_plus][j][k][2] - u0[i_minus][j][k][2]);

				W[i][j][k][0][1] =  0.5*(u0[i][j_plus][k][0] - u0[i][j_minus][k][0]);
				W[i][j][k][1][1] =  0.5*(u0[i][j_plus][k][1] - u0[i][j_minus][k][1]);
				W[i][j][k][2][1] =  0.5*(u0[i][j_plus][k][2] - u0[i][j_minus][k][2]);

				if(k==0 && wall_on!=0){
					W[i][j][k][0][2] = third * u0[i][j][k_plus][0] + u0[i][j][k][0];
					W[i][j][k][1][2] = third * u0[i][j][k_plus][1] + u0[i][j][k][1] - four3rd * uy_bottom;
					W[i][j][k][2][2] = third * u0[i][j][k_plus][2] + u0[i][j][k][2];
				} else if(k==Nz-1 && wall_on!=0) {
					W[i][j][k][0][2] =-third * u0[i][j][k_minus][0] - u0[i][j][k][0];
					W[i][j][k][1][2] =-third * u0[i][j][k_minus][1] - u0[i][j][k][1]+ four3rd * uy_top;
					W[i][j][k][2][2] =-third * u0[i][j][k_minus][2] - u0[i][j][k][2];
				} else {
					W[i][j][k][0][2] =  0.5 * (u0[i][j][k_plus][0] - u0[i][j][k_minus][0]);
					W[i][j][k][1][2] =  0.5 * (u0[i][j][k_plus][1] - u0[i][j][k_minus][1]);
					W[i][j][k][2][2] =  0.5 * (u0[i][j][k_plus][2] - u0[i][j][k_minus][2]);
				}
			}
		}
	}
	
/*	printf("W:");
	for(i=0;i<3;i++){
		for(j=0;j<3;j++)printf(" %f",W[0][0][0][i][j]);
	}
	printf("\n");
 */
}


void cal_rho(real rho_local[Nx][Ny][Nz], real f[Nx][Ny][Nz][15])
{
	int i,j,k,ii=0;
	real *f0, temp;
	
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){
				f0=&f[i][j][k][0];
				temp=0.;
				for(ii=0;ii<15;ii++)temp+=f0[ii];
				rho_local[i][j][k]=temp;
			}
		}
	}
}

void cal_u(real u0[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15])
{
	int ii=0,i,j,k,l;
	real *f0, irho;
	
	cal_mom(u0,f);
	cal_rho(Rho,f);
	
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){
				irho = 1.0/Rho[i][j][k];
				u0[i][j][k][0]=u0[i][j][k][0] * irho;
				u0[i][j][k][1]=u0[i][j][k][1] * irho;
				u0[i][j][k][2]=u0[i][j][k][2] * irho;
			}
		}
	}

}


void cal_mom(real m[Nx][Ny][Nz][3], real f[Nx][Ny][Nz][15])
{
	int ii=0,i,j,k,l;
	real *f0, m0, m1, m2;
	
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){
				f0=&f[i][j][k][0];
				m0=0.0;
				m1=0.0;
				m2=0.0;
				for (ii=0; ii<15; ii++) {
					m0 += f0[ii] * e[ii][0];
					m1 += f0[ii] * e[ii][1];
					m2 += f0[ii] * e[ii][2];
				}
				m[i][j][k][0] = m0;
				m[i][j][k][1] = m1;
				m[i][j][k][2] = m2;
			}
		}
	}
	
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
			printf("error: index of bounce() is out of range!");
			return -1;
	}
}


real QQ(real T[3][3], int i, int j)
{
//	real qq=0.0;
//	qq = T[i][0]*T[0][j];

	return T[i][0]*T[0][j]+T[i][1]*T[1][j]+T[i][2]*T[2][j];
}


void write_restart(real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3])
{
	FILE* frestart;
	int i, j, k, ii, jj;

	frestart=fopen("restart.dat","wb");

	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
				for(ii=0; ii<3; ii++){
					for(jj=0; jj<3; jj++){
						fwrite(&Q[i][j][k][ii][jj],sizeof(real),1,frestart);
					}
				}
			}
		}
	}

	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
				fwrite(&Rho[i][j][k],sizeof(real),1,frestart);
			}
		}
	}

	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
//		    	u[i][j][k][0]=0;
//		    	u[i][j][k][2]=0;
				for(ii=0; ii<3; ii++)fwrite(&u[i][j][k][ii],sizeof(real),1,frestart);
			}
		}
	}

	fclose(frestart);
}

void read_restart(real Q[Nx][Ny][Nz][3][3], real Rho[Nx][Ny][Nz], real u[Nx][Ny][Nz][3])
{
	FILE* frestart;

	int i,j,k,ii,jj;
	frestart=fopen("restart.dat","rb");


	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
				for(ii=0; ii<3; ii++)for(jj=0; jj<3; jj++)fread(&Q[i][j][k][ii][jj],sizeof(double),1,frestart);
			}
		}
	}

	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
				fread(&Rho[i][j][k],sizeof(double),1,frestart);
//				Rho[i][j][k]=rho;
			}
		}
	}

	for(i=0; i<Nx; i++){
	    for(j=0; j<Ny; j++){
		    for(k=0; k<Nz; k++){
				for(ii=0; ii<3; ii++){
					fread(&u[i][j][k][ii],sizeof(double),1,frestart);
				}
			}
		}
	}


	fclose(frestart);
}


void print_visual(int new)
{
	int i, j, k,  ip, jp, kp, im, jm, km, ipar, inside;
	FILE *qfile, *ufile, *grid, *sfile, *dfile;
	real dx, dy, dz, rsq;
	
	if(new==1) {
		qfile=fopen("Q_3d","w");
		ufile=fopen("u_3d","w");
		grid =fopen("grid","w");
		sfile=fopen("stress_3d","w");
		dfile=fopen("rho_3d","w");
	} else {
		qfile=fopen("Q_3d","a");
        	ufile=fopen("u_3d","a");
        	sfile=fopen("stress_3d","a");
        	dfile=fopen("rho_3d","a");
	}	

	if(new==1)fprintf(grid,"Nx Ny Nz %d %d %d\n",Nx,Ny,Nz);
	
	for(k=0;k<Nz;k++){
	    for(j=0;j<Ny;j++){
	    	for(i=0;i<Nx;i++){
	    		inside=-1;
	    		for(ipar=0;ipar<npar;ipar++){
	    			rsq= p_rad[ipar]*p_rad[ipar];
	    			dx = (real)i - p_pos[ipar][0];
	    			dy = (real)j - p_pos[ipar][1];
	    			dz = (real)k - p_pos[ipar][2];
	    			if(dx < -0.5*Nx) dx += (real)Nx;
					if(dx >  0.5*Nx) dx -= (real)Nx;
					if(dy < -0.5*Ny) dy += (real)Ny;
					if(dy >  0.5*Ny) dy -= (real)Ny;
					if(dx*dx+dy*dy+dz*dz<rsq)inside=1;
	    		}
	    		if(new==1)fprintf(grid,"%d %d %d %d\n",i,j,k,inside);
				fprintf(qfile,"%e %e %e %e %e\n",Q1[i][j][k][0][0],Q1[i][j][k][0][1],Q1[i][j][k][0][2],Q1[i][j][k][1][1],Q1[i][j][k][1][2]);
	    		fprintf(ufile,"%e %e %e\n",u[i][j][k][0],u[i][j][k][1],u[i][j][k][2]);
	    		fprintf(sfile,"%e %e %e\n",sigma_p[i][j][k][0],sigma_p[i][j][k][1],sigma_p[i][j][k][2]);
	    		fprintf(dfile,"%e\n",Rho[i][j][k]);
	    	}
	    }
	}
	
	fclose(qfile);
	fclose(ufile);
	if(new==1)fclose(grid);
	fclose(sfile);
	fclose(dfile);
}


void period_image(int *i, int *j, int *k)
{
	if (*i<0) {
		*i += Nx;
	}else if (*i>=Nx){
		*i -= Nx;
	}
	
	if (*j<0) {
		*j += Ny;
	}else if (*j>=Ny){
		*j -= Ny;
	}
	
	if (*k<0) {
		if (wall_on==0) {
			*k += Nz;
		}else {
			*k = 0;
			printf("warning: in period_image(), lattic point out of box!\n");
		}

	}else if (*k>=Nz){
		if (wall_on==0) {
			*k -= Nz;
		}else {
			*k = Nz-1;
			printf("warning: in period_image(), lattic point out of box!\n");
		}
	}
}


void cal_fequ(real f_eq[Nx][Ny][Nz][15], real Rho0[Nx][Ny][Nz], real u0[Nx][Ny][Nz][3], real sigma0[Nx][Ny][Nz][3][3])
{
	int i,j,k,ii=0,jj=0,kk=0;
	
	real E2[3][3], E1[3][3];
	real trsigma=0, temp;
	real A=0, A0=0, A1=0, A2=0;
	real B=0, B0=0, B1=0, B2=0;
	real C=0, C0=0, C1=0, C2=0;
	real D=0, D0=0, D1=0, D2=0;
	real temp0=0;
	
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
				trsigma =sigma0[i][j][k][0][0] + sigma0[i][j][k][1][1] + sigma0[i][j][k][2][2];	  
				A2 = -one30th*trsigma;
				A1 =  A2;
				A0 =  Rho0[i][j][k] - 14.0 * A2;
				
				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						E2[ii][jj] = 0.0625*(-sigma0[i][j][k][ii][jj] + trsigma*OneThirdDelta(ii,jj));
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
						temp += B*u0[i][j][k][jj]*e[ii][jj]+C*u0[i][j][k][jj]*u0[i][j][k][jj];
						
						for(kk=0;kk<3;kk++){
							temp += D*u0[i][j][k][jj]*u0[i][j][k][kk]*e[ii][jj]*e[ii][kk];
						}
					}
					f_eq[i][j][k][ii] = temp*rho;
					temp0 += f_eq[i][j][k][ii];
				}
				
			}
		}
	}
//	printf("rrrrho=%20.15f\n",temp0);
	
}


void cal_feqf_new(real f_eq0[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real sigma[Nx][Ny][Nz][3][3])
{
	int i,j,k,ii=0,jj=0,kk=0;
	int compressibility=1;		// 0: incompressible flow; nonzero: compressible flow
	
	real E[3][3], E1, E2, rho0=rho;
	real trsigma=0, trsigmar=0, temp;
	real A=0, A0=0, A1=0, A2=0;
	real B=0, B0=0, B1=0, B2=0;
	real C=0, C0=0, C1=0, C2=0;
	real D=0, D0=0, D1=0, D2=0;
	real udote, usq, udote_sq;
	real temp1, temp2;
	
	B2 = one24th;
	B1 = third;
	//	B0 = 0.0;
	
	C2 = -one24th;
	C1 = -one12th;
	C0 = -two3rd;
	
	D2 = 0.0625;
	D1 = 0.5;
	//	D0 = 0.0;
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				
				trsigma = sigma[i][j][k][0][0] + sigma[i][j][k][1][1] + sigma[i][j][k][2][2];
				A2 = -one30th*trsigma;
				A1 = A2;
				A0 = Rho[i][j][k] - 14.0 * A2;
				
				usq = 0;
				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						E[ii][jj] = 0.0625*(-sigma[i][j][k][ii][jj] + trsigma*OneThirdDelta(ii,jj) );
					}
					usq += u[i][j][k][ii]*u[i][j][k][ii];
		  		}
				
		  		//Calculate f_eq[i][j][k] 
		  		for(ii=0;ii<15;ii++){
					udote = 0;
					E2    = 0;
					for(jj=0; jj<3; jj++){
						udote += u[i][j][k][jj]*e[ii][jj];
						for(kk=0; kk<3; kk++){
							E2 += E[jj][kk] * e[ii][jj] * e[ii][kk];
						}
					}
					udote_sq = udote * udote;
					E1       =   8.0 * E2;
					if(compressibility!=0){
						rho0 = Rho[i][j][k];
					}
					if(ii==0){
						f_eq0[i][j][k][ii] = A0 + rho0 * (           C0*usq                   );
					} else if(ii<7) {
						f_eq0[i][j][k][ii] = A1 + rho0 * (B1*udote + C1*usq + D1*udote_sq + E1);
					} else {
						f_eq0[i][j][k][ii] = A2 + rho0 * (B2*udote + C2*usq + D2*udote_sq + E2);
					}
					
				}
			}
		}
	}

}


void cal_sigma(real Q0[Nx][Ny][Nz][3][3])
{
	int i,j,k,ii=0,jj=0, mm=0, nn=0;
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				for(ii=0;ii<3;ii++) {
					for(jj=0;jj<3;jj++){
//						sigma[i][j][k][ii][jj] = -Rho[i][j][k] * T * Delta(ii,jj) + sigma_q[i][j][k][ii][jj];
                                                sigma[i][j][k][ii][jj] = -Rho[i][j][k] * T * Delta(ii,jj);
					}
				}
			}
		}
	}
}


void cal_p(real p[Nx][Ny][Nz][15], real f[Nx][Ny][Nz][15], real feq[Nx][Ny][Nz][15])
{
	int i,j,k,ii=0,jj=0,kk=0,ipar,in_particle;
	real temp, td, tu, dx, dy, dz, T1, T2;
	
	T2 = 1.0/24.0;
	T1 = 8*T2;
	//	T0 = 0.0;
	
	//	cal_detau();
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				in_particle=0;
				if(particle_on!=0){					
					for(ipar=0; ipar<npar; ipar++){
						dx = i - p_pos[ipar][0];
						dy = j - p_pos[ipar][1];
						dz = k - p_pos[ipar][2];
						if(dx<-0.5*Lx)dx+=Lx;
						if(dx> 0.5*Lx)dx-=Lx;
						if(dy<-0.5*Ly)dy+=Ly;
						if(dy> 0.5*Ly)dy-=Ly;
						if(wall_on==0){
							if(dz<-0.5*Lz)dz+=Lz;
							if(dz> 0.5*Lz)dz-=Lz;
						}
						if(dx*dx+dy*dy+dz*dz<=p_rad[ipar]*p_rad[ipar])in_particle=1;
					}
				}
				if(in_particle==0){
					for(ii=0;ii<15;ii++){
						td = 0;
						tu = 0;
						for (jj=0; jj<15; jj++) {
							td += f[i][j][k][jj];
							tu += f[i][j][k][jj]*e[jj][1];
						}
						for (ii=0; ii<15; ii++) {
//							temp = 3.0*(1-0.5*itau_f)*(e[ii][1]-tu*(1.0/td))*yforce*feq[i][j][k][ii];
//							p[i][j][k][ii] = temp;
							temp = yforce * e[ii][1];
							if(ii>0 && ii<7){
								p[i][j][k][ii] = T1 * temp;
							}else if(ii>6){
								p[i][j][k][ii] = T2 * temp;
							}else {
								p[i][j][k][ii] = 0 * temp;
							}
							
							if(t_current>-1){
//							if(t_current==-1){
								if(ii>0 && ii<7){
									p[i][j][k][ii] += T1 * (sigma_p[i][j][k][0]*e[ii][0] +sigma_p[i][j][k][1]*e[ii][1] +sigma_p[i][j][k][2]*e[ii][2] ); 
								}
								else if (ii>6){
									p[i][j][k][ii] += T2 * (sigma_p[i][j][k][0]*e[ii][0] +sigma_p[i][j][k][1]*e[ii][1] +sigma_p[i][j][k][2]*e[ii][2] );
								}
							}
						}
					}
		  		}
			}
		}
	}
	/*	
	 printf("p:");
	 for(i=0;i<3;i++){
	 printf(" %f",p[0][0][0][i]);
	 }
	 printf("\n");
	 */
	/*	
	 if(t_current-t_max==-1) {
	 printf("p:");
	 for(i=0;i<Nz;i++)printf(" %21.15le %21.15le %21.15le\n",p[Nx/2][Ny/2][i][0],p[Nx/2][Ny/2][i][1],p[Nx/2][Ny/2][i][2]);
	 printf("\n");
	 }
	 */
}


void monitor(int flag, int print)
{
	int i, j, k, ii;
	real utotal[3], rhototal=0.0, rhotemp, utemp, vtemp, wtemp;
	real q_diff_sq=0, k_eng_new=0, k_diff;
		
	if(print!=0)printf("duz/dz=%20.15f\n",W[Nx/2][Ny/2][0][2][2]);
		
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){
				if (flag<=0) {
					q_diff_sq+=(Q1[i][j][k][0][0]-Q2[i][j][k][0][0])*(Q1[i][j][k][0][0]-Q2[i][j][k][0][0]);
					q_diff_sq+=(Q1[i][j][k][0][1]-Q2[i][j][k][0][1])*(Q1[i][j][k][0][1]-Q2[i][j][k][0][1]);
					q_diff_sq+=(Q1[i][j][k][0][2]-Q2[i][j][k][0][2])*(Q1[i][j][k][0][2]-Q2[i][j][k][0][2]);
					q_diff_sq+=(Q1[i][j][k][1][1]-Q2[i][j][k][1][1])*(Q1[i][j][k][1][1]-Q2[i][j][k][1][1]);
					q_diff_sq+=(Q1[i][j][k][1][2]-Q2[i][j][k][1][2])*(Q1[i][j][k][1][2]-Q2[i][j][k][1][2]);
				}
				if (flag>=0) {
					k_eng_new+=(u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2])*Rho[i][j][k];
				}				
			}
		}
	}
	k_eng_new *= 0.5;
	
	if (print!=0) {
		utotal[0]=0.0;
		utotal[1]=0.0;
		utotal[2]=0.0;
		for(i=0;i<Nx;i++){
			for(j=0;j<Ny;j++){
				for(k=0;k<Nz;k++){
					utotal[0]+=u[i][j][k][0];
					utotal[1]+=u[i][j][k][1];
					utotal[2]+=u[i][j][k][2];
					rhototal+=Rho[i][j][k];
				}
			}
		}
		printf("vel=%20.15f %20.15f %20.15f\n",utotal[0], utotal[1], utotal[2]);
		printf("rho=%20.15f\n",rhototal);
	}
	
	if (flag<=0) {		// monitor Q
		if(print!=0)printf("dif_Q^2=%30.26f\n",q_diff_sq);
		if(q_diff_sq != q_diff_sq) exit(1); //check if q_diff_sq is NaN		
		if(q_diff_sq<qthreshold)Q_on=0;
	}
		
	if(flag>=0){		// monitor f's
		if(k_eng<1e-15 && k_eng>-1e-15){
			k_diff = k_eng_new;
		}else{
			k_diff = (k_eng_new-k_eng)/k_eng;
		}
		if((k_eng-k_eng_new<1e-18 && k_eng-k_eng_new>-1e-18)||(k_diff<1e-13 && k_diff>-1e-13)) flow_on=0;
		if(print!=0)printf("diff_k =%30.25f\n",k_diff);
		k_eng=k_eng_new;
	}
}


void print_u1d(FILE *FILE_rho, FILE *FILE_ux, FILE *FILE_uy, FILE *FILE_uz)
{
	int i, j, k, ii;
	real utotal[3], rhototal=0.0, rhotemp, utemp, vtemp, wtemp;
	real q_diff_sq=0, k_eng_new=0, k_diff;
	
	for(k=0;k<Nz;k++){
		rhotemp=0.0;
		utemp=0.0;
		vtemp=0.0;
		wtemp=0.0;
//		ptemp=0.0;
		for(j=Ny/2;j<Ny/2+1;j++){
			for(i=Nx/2;i<Nx/2+1;i++){
				rhotemp += Rho[i][j][k];
				utemp +=   u[i][j][k][0];
				vtemp +=   u[i][j][k][1];
				wtemp +=   u[i][j][k][2];
//						  ptemp +=  Po[i][j][k];
			}
		}
//				ptemp = Po[Nx/2][Ny/2][kk];
		fprintf(FILE_rho,"%20.15f ",rhotemp);
		fprintf(FILE_ux ,"%20.15f ",utemp);
		fprintf(FILE_uy ,"%20.15f ",vtemp);
		fprintf(FILE_uz ,"%20.15f ",wtemp);
//				fprintf(out5,"%20.15f ",ptemp);
	}
	fprintf(FILE_rho,"\n");
	fprintf(FILE_ux,"\n");
	fprintf(FILE_uy,"\n");
	fprintf(FILE_uz,"\n");
//			fprintf(out5,"\n");
//		  printf("\n");
//			print_Q(1);
//			print_Q2(1);
}


void print_Q1do(int flag)
{
	int i,j,k,ii=0;
	real q1, q2, q3, q4, q5;
	FILE *fout1, *fout2, *fout3, *fout4, *fout5;
	
	if(flag==0){
		fout1=fopen("q1","w");
		fout2=fopen("q2","w");
		fout3=fopen("q3","w");
		fout4=fopen("q4","w");
		fout5=fopen("q5","w");
	} else {
		fout1=fopen("q1","a");
		fout2=fopen("q2","a");
		fout3=fopen("q3","a");
		fout4=fopen("q4","a");
		fout5=fopen("q5","a");
	}
	
	for(k=0;k<Nz;k++){
		q1=0.;
		q2=0.;
		q3=0.;
		q4=0.;
		q5=0.;
		for(i=0;i<Nx;i++){
			for(j=0;j<Ny;j++){		
				q1 += Q1[i][j][k][0][0];
				q2 += Q1[i][j][k][0][1];
				q3 += Q1[i][j][k][0][2];
				q4 += Q1[i][j][k][1][1];
				q5 += Q1[i][j][k][1][2];				
			}
		}
		q1=q1/(real)Nx/(real)Ny;
		q2=q2/(real)Nx/(real)Ny;
		q3=q3/(real)Nx/(real)Ny;
		q4=q4/(real)Nx/(real)Ny;
		q5=q5/(real)Nx/(real)Ny;
		fprintf(fout1,"%21.15f ",q1);
		fprintf(fout2,"%21.15f ",q2);
		fprintf(fout3,"%21.15f ",q3);
		fprintf(fout4,"%21.15f ",q4);
		fprintf(fout5,"%21.15f ",q5);
	}
	fprintf(fout1,"\n");
	fprintf(fout2,"\n");
	fprintf(fout3,"\n");
	fprintf(fout4,"\n");
	fprintf(fout5,"\n");
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);
	fclose(fout5);
}



void print_Q1d(int flag)
{
	int i,j,k,ii=0;
	real q1, q2, q3, q4, q5;
	FILE *fout1, *fout2, *fout3, *fout4, *fout5;
	
	if(flag==0){
		fout1=fopen("q1","w");
		fout2=fopen("q2","w");
		fout3=fopen("q3","w");
		fout4=fopen("q4","w");
		fout5=fopen("q5","w");
	} else {
		fout1=fopen("q1","a");
		fout2=fopen("q2","a");
		fout3=fopen("q3","a");
		fout4=fopen("q4","a");
		fout5=fopen("q5","a");
	}
	
	for(k=0;k<Nz;k++){
		q1=0.;
		q2=0.;
		q3=0.;
		q4=0.;
		q5=0.;
		//		for(i=0;i<Nx;i++){
		//			for(j=0;j<Ny;j++){
		i = Nx/2;
		j = Ny/2;
		q1 += Q1[i][j][k][0][0];
		q2 += Q1[i][j][k][0][1];
		q3 += Q1[i][j][k][0][2];
		q4 += Q1[i][j][k][1][1];
		q5 += Q1[i][j][k][1][2];				
		//			}
		//		}
		/*		q1=q1/(real)Nx/(real)Ny;
		 q2=q2/(real)Nx/(real)Ny;
		 q3=q3/(real)Nx/(real)Ny;
		 q4=q4/(real)Nx/(real)Ny;
		 q5=q5/(real)Nx/(real)Ny;*/
		fprintf(fout1,"%21.15le ",q1);
		fprintf(fout2,"%21.15le ",q2);
		fprintf(fout3,"%21.15le ",q3);
		fprintf(fout4,"%21.15le ",q4);
		fprintf(fout5,"%21.15le ",q5);
	}
	fprintf(fout1,"\n");
	fprintf(fout2,"\n");
	fprintf(fout3,"\n");
	fprintf(fout4,"\n");
	fprintf(fout5,"\n");
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);
	fclose(fout5);
}

void print_Q2d(int dim, int x, int flag)
{
	int i, j, k, ii, jj, iihi, jjhi;
	real q1, q2, q3, q4, q5;
	FILE *fout1, *fout2, *fout3, *fout4, *fout5;
	
	if(flag==0){
		fout1=fopen("Q1","w");
		fout2=fopen("Q2","w");
		fout3=fopen("Q3","w");
		fout4=fopen("Q4","w");
		fout5=fopen("Q5","w");
	} else {
		fout1=fopen("Q1","a");
		fout2=fopen("Q2","a");
		fout3=fopen("Q3","a");
		fout4=fopen("Q4","a");
		fout5=fopen("Q5","a");
	}
	
	if (dim==0) {
		i = x;
		iihi = Ny;
		jjhi = Nz;
	} else if (dim==1) {
		j = x;
		iihi = Nx;
		jjhi = Nz;
	} else if (dim==2) {
		k = x;
		iihi = Nx;
		jjhi = Ny;
	} else {
		printf("error: using print_Q3!\n");
	}
	
	
	for(ii=0;ii<iihi;ii++){
		for(jj=0;jj<jjhi;jj++){
			if (dim==0) {
				j = ii;
				k = jj;
			} else if (dim==1) {
				i = ii;
				k = jj;
			} else if (dim==2) {
				i = ii;
				j = jj;
			} 
			q1 = Q1[i][j][k][0][0];
			q2 = Q1[i][j][k][0][1];
			q3 = Q1[i][j][k][0][2];
			q4 = Q1[i][j][k][1][1];
			q5 = Q1[i][j][k][1][2];				
			
			fprintf(fout1,"%21.15le ",q1);
			fprintf(fout2,"%21.15le ",q2);
			fprintf(fout3,"%21.15le ",q3);
			fprintf(fout4,"%21.15le ",q4);
			fprintf(fout5,"%21.15le ",q5);
		}
		fprintf(fout1,"\n");
		fprintf(fout2,"\n");
		fprintf(fout3,"\n");
		fprintf(fout4,"\n");
		fprintf(fout5,"\n");
	}
	
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);
	fclose(fout5);
	
}



void print_u2d(int dim, int x, int flag)
{
	int i, j, k, ii, jj, iihi, jjhi, ix, iy;
	real cuy, cuz, c_rho;
	FILE *fout1, *fout2, *fout3;
	
	if(flag==0){
		fout1=fopen("cuy","w");
		fout2=fopen("cuz","w");
		fout3=fopen("crho","w");
	} else {
		fout1=fopen("cuy","a");
		fout2=fopen("cuz","a");
		fout3=fopen("crho","a");
	}
	ix = 1;
	iy = 2;
	
	if (dim==0) {
		i = x;
		iihi = Ny;
		jjhi = Nz;
	} else if (dim==1) {
		j = x;
		iihi = Nx;
		jjhi = Nz;
	} else if (dim==2) {
		k = x;
		iihi = Ny;
		jjhi = Nz;
	} else {
		printf("error: using print_u2d!\n");
	}
	
	
	for(ii=0;ii<iihi;ii++){
		for(jj=0;jj<jjhi;jj++){
			if (dim==0) {
				j = ii;
				k = jj;
			} else if (dim==1) {
				i = ii;
				k = jj;
			} else if (dim==2) {
				i = ii;
				j = jj;
			} 
			cuy = u[i][j][k][ix];
			cuz = u[i][j][k][iy];
			c_rho = Rho[i][j][k];
			
			fprintf(fout1,"%21.15le ",cuy);
			fprintf(fout2,"%21.15le ",cuz);
			fprintf(fout3,"%21.15le ",c_rho);
		}
		fprintf(fout1,"\n");
		fprintf(fout2,"\n");
		fprintf(fout3,"\n");
	}
	
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
}
