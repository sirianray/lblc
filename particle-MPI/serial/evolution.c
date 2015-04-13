/*
 *  evolution.c
 *  
 *
 *  Created by Sirius on 3/2/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */
#include "main.h"
#include "particle.h"

void evol_f(real *fin, real *fout)
{
	int id, ip;
	double Cf1, Cf2;
	
//	cal_sigma();
	cal_p();
	streaming(fin,fout);
	
	cal_feqf_diff(fin,Cf);
	for (id=0; id<points*15; id++) {
		Cf1      = -itau_f * Cf[id] + p[id];
		ip       = nextf[id];
		fin[id] += 0.5 * dt * Cf1;
		fout[ip]+= dt * Cf1;
	}
	
//	cal_rho(fout);
	cal_u(fout);
//	cal_sigma();
	cal_p();
	cal_feqf_diff(fout,fout);
	for (id=0; id<points*15; id++) {
		ip       = nextf[id];
		Cf2      = -itau_f * fout[ip] + p[ip];
		fin[id] += dt * 0.5 * Cf2;
	}
	
	streaming(fin,fout);
	
//	cal_rho(fout);
	cal_u(fout);
}

void streaming(real *fin, real *fout)
{
	int ip, id, idn, k, kn, ii;
	double Bx2;
	
	for (ip=0; ip<points; ip++) {
		for (ii=0; ii<15; ii++) {
			id = ip*15+ii;
			idn= nextf[id];
			fout[idn]=fin[id];
			if (wall_on!=0) {
				k = ip/bulk0;
				kn = k + e[ii][2];
				if ( kn==-1 || kn==Nz ) {
					if(ii<7) {
						Bx2 = two3rd;
					} else {
						Bx2 = one12th;
					}
					if (k==0) {
						fout[idn] += -Bx2*Rho[ip]*(uy_bottom*e[ii][1]);
					} else {
						fout[idn] += -Bx2*Rho[ip]*(uy_top*e[ii][1]);
					}
				}
			}
			if (npar>0 && idn/15==ip) {
				if(ii<7) {
					Bx2 = two3rd  * rho;
				} else {
					Bx2 = one12th * rho;
				}
				fout[idn] += -Bx2 * ubounce[id];
			}
		}
	}
}

void cal_W()
{
	int id, ip, iw, idxm, idym, idzm, idxp, idyp, idzp, bulk3 = 3*bulk0;
	
	for (id=0; id<points*3; id+=3) {
		ip   = 5*id;
		idxp = (int)(nextf[ip+2]/15)*3;
		idxm = (int)(nextf[ip+1]/15)*3;
		idyp = (int)(nextf[ip+4]/15)*3;
		idym = (int)(nextf[ip+3]/15)*3;
		idzp = (int)(nextf[ip+6]/15)*3;
		idzm = (int)(nextf[ip+5]/15)*3;
		
		iw   = 3*id;
		
		if (npar==0 || idxm!=id && idxp!=id) {
			W[iw]  = 0.5 * (u[idxp]  - u[idxm]);
			W[iw+1]= 0.5 * (u[idxp+1]- u[idxm+1]);
			W[iw+2]= 0.5 * (u[idxp+2]- u[idxm+2]);
		} else if (idxm==id) {			
			W[iw]  =  third*u[idxp]  + u[id]  ;
			W[iw+1]=  third*u[idxp+1]+ u[id+1];
			W[iw+2]=  third*u[idxp+2]+ u[id+2];
		} else if (idxp==id) {
			W[iw]  =-(u[id]  + third*u[idxm]  );
			W[iw+1]=-(u[id+1]+ third*u[idxm+1]);
			W[iw+2]=-(u[id+2]+ third*u[idxm+2]);
		}
		
		if (npar==0 || idym!=id && idyp!=id) {
			W[iw+3]= 0.5 * (u[idyp]  - u[idym]);
			W[iw+4]= 0.5 * (u[idyp+1]- u[idym+1]);
			W[iw+5]= 0.5 * (u[idyp+2]- u[idym+2]);
		} else if (idyp==id) {
			W[iw+3]=-(u[id]  + third*u[idym]  );
			W[iw+4]=-(u[id+1]+ third*u[idym+1]);
			W[iw+5]=-(u[id+2]+ third*u[idym+2]);
		} else if (idym==id) {
			W[iw+3]=  third*u[idyp]  + u[id]  ;
			W[iw+4]=  third*u[idyp+1]+ u[id+1];
			W[iw+5]=  third*u[idyp+2]+ u[id+2];
		}
		
		if (idzp!=id && idzm!=id) {
			W[iw+6]= 0.5 * (u[idzp]  - u[idzm]);
			W[iw+7]= 0.5 * (u[idzp+1]- u[idzm+1]);
			W[iw+8]= 0.5 * (u[idzp+2]- u[idzm+2]);
		} else if (wall_on!=0 && id/bulk3==0) {
			W[iw+6]=  u[id]   + third*u[idzp];
			W[iw+7]=  u[id+1] + third*u[idzp+1] - four3rd * uy_bottom;
			W[iw+8]=  u[id+2] + third*u[idzp+2];
		} else if (wall_on!=0 && id/bulk3==Nz-1) {
			W[iw+6]= -u[id]   - third*u[idzm];
			W[iw+7]= -u[id+1] - third*u[idzm+1] + four3rd * uy_top;
			W[iw+8]= -u[id+2] - third*u[idzm+2];
		} else if (npar>0 && idzp==id) {
			W[iw+6]=-(u[id]  + third*u[idzm]  );
			W[iw+7]=-(u[id+1]+ third*u[idzm+1]);
			W[iw+8]=-(u[id+2]+ third*u[idzm+2]);
		} else if (npar>0 && idzm==id) {
			W[iw+6]=  third*u[idzp]  + u[id]  ;
			W[iw+7]=  third*u[idzp+1]+ u[id+1];
			W[iw+8]=  third*u[idzp+2]+ u[id+2];
		}
	}
}

void cal_rho(real *f0)
{
	int i, j, k;
	
	for (i=0; i<points; i++) {
		k      = i*15;
		Rho[i] = f0[k];
		for (j=1; j<15; j++) {
			Rho[i] += f0[k+j];
		}
	}
}

void cal_u(real *f0)
{
	int i, j, k, l, m;
	double irho;
	
	cal_rho(f0);
	
	for (i=0; i<points; i++) {
		irho  = 1.0/Rho[i];
		k     = i*3;
		l     = k*5;
		u[k]  = f0[l]*e[0][0];
		u[k+1]= f0[l]*e[0][1];
		u[k+2]= f0[l]*e[0][2];
		for (j=1; j<15; j++) {
			m      = l+j;
			u[k]  += f0[m]*e[j][0];
			u[k+1]+= f0[m]*e[j][1];
			u[k+2]+= f0[m]*e[j][2];
		}

		u[k]  *= irho;
		u[k+1]*= irho;
		u[k+2]*= irho;
	}
}


void cal_fequ(real *f0)
{
	int ii=0,jj=0,kk=0,id;
	
	real E2[3][3], E1[3][3];
	real trsigma=0, temp;
	real A=0, A0=0, A1=0, A2=0;
	real B=0, B0=0, B1=0, B2=0;
	real C=0, C0=0, C1=0, C2=0;
	real D=0, D0=0, D1=0, D2=0;
	real temp0=0;
	
	double sigma[3][3]={0.};
	
	B2 = one24th;
	B1 = third;
	B0 = 0.0;
	
	C2 = -one24th;
	C1 = -one12th;
	C0 = -two3rd;
	
	D2 = 0.0625;
	D1 = 0.5;
	D0 = 0.0;
	

	for (id=0; id<points; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
		A2 = -one30th*trsigma;
		A1 =  A2;
		A0 =  Rho[id] - 14.0 * A2;
		
		for(ii=0;ii<3;ii++){
			for(jj=0;jj<3;jj++){
				E2[ii][jj] = 0.0625*(-sigma[ii][jj] + trsigma*One3rdDelta(ii,jj));
				E1[ii][jj] = 8.0*E2[ii][jj];
			}
		}
		
		//Calculate f_eq 
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
				temp += B*u[id*3+jj]*e[ii][jj]+C*u[id*3+jj]*u[id*3+jj];
				
				for(kk=0;kk<3;kk++){
					temp += D*u[id*3+jj]*u[id*3+kk]*e[ii][jj]*e[ii][kk];
				}
			}
			f0[id*15+ii] = temp*rho;
		}
	}
}

void cal_feqf_new(real *f0)
{
	int i,j,k,ii=0,jj=0,kk=0, id;
	int compressibility=1;		// 0: incompressible flow; nonzero: compressible flow
	
	real E[3][3], E1, E2, rho0=rho;
	real trsigma=0, trsigmar=0, temp;
	real A=0, A0=0, A1=0, A2=0;
	real B=0, B0=0, B1=0, B2=0;
	real C=0, C0=0, C1=0, C2=0;
	real D=0, D0=0, D1=0, D2=0;
	real udote, usq, udote_sq;
	real temp1, temp2;
	
	double sigma[3][3]={0.};
	
	B2 = one24th;
	B1 = third;
//	B0 = 0.0;
	
	C2 = -one24th;
	C1 = -one12th;
	C0 = -two3rd;
	
	D2 = 0.0625;
	D1 = 0.5;
//	D0 = 0.0;
	
	for (id=0; id<points; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
		A2 = -one30th*trsigma;
		A1 = A2;
		A0 = Rho[id] - 14.0 * A2;
		
		usq = 0;
		for(ii=0;ii<3;ii++){
			for(jj=0;jj<3;jj++){
				E[ii][jj] = 0.0625*(-sigma[ii][jj] + trsigma*One3rdDelta(ii,jj) );
			}
			usq += u[id*3+ii]*u[id*3+ii];
		}
		
		//Calculate f_eq
		for(ii=0;ii<15;ii++){
			udote = 0;
			E2    = 0;
			for(jj=0; jj<3; jj++){
				udote += u[id*3+jj]*e[ii][jj];
				for(kk=0; kk<3; kk++){
					E2 += E[jj][kk] * e[ii][jj] * e[ii][kk];
				}
			}
			udote_sq = udote * udote;
			E1       =   8.0 * E2;
			if(compressibility!=0){
				rho0 = Rho[id];
			}
			if(ii==0){
				f0[id*15+ii] = A0 + rho0 * (           C0*usq                   );
			} else if(ii<7) {
				f0[id*15+ii] = A1 + rho0 * (B1*udote + C1*usq + D1*udote_sq + E1);
			} else {
				f0[id*15+ii] = A2 + rho0 * (B2*udote + C2*usq + D2*udote_sq + E2);
			}
			
		}
	}
	
}


void cal_feqf_diff(real *fin, real *fout)
{
	int i,j,k,ii=0,jj=0,kk=0, id;
	int compressibility=1;		// 0: incompressible flow; nonzero: compressible flow
	
	real E[3][3]={0.}, E1, E2, rho0=rho;
	real trsigma=0, trsigmar=0, temp;
	real A=0, A0=0, A1=0, A2=0;
	real B=0, B0=0, B1=0, B2=0;
	real C=0, C0=0, C1=0, C2=0;
	real D=0, D0=0, D1=0, D2=0;
	real udote, usq, udote_sq;
	real temp1, temp2;
	
	double sigma[3][3]={0.};
	
	B2 = one24th;
	B1 = third;
//	B0 = 0.0;
	
	C2 = -one24th;
	C1 = -one12th;
	C0 = -two3rd;
	
	D2 = 0.0625;
	D1 = 0.5;
//	D0 = 0.0;
	
	for (id=0; id<points; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
		A2 = -one30th*trsigma;
		A1 = A2;
		A0 = Rho[id] - 14.0 * A2;
		
		usq = 0;
		for(ii=0;ii<3;ii++){
//			for(jj=0;jj<3;jj++){
//				E[ii][jj] = 0.0625*(-sigma[id*9+ii*3+jj] + trsigma*One3rdDelta(ii,jj) );
//			}
			usq += u[id*3+ii]*u[id*3+ii];
		}
		
//		Calculate f_eq
		for(ii=0;ii<15;ii++){
			udote = 0;
			E2    = 0;
			for(jj=0; jj<3; jj++){
				udote += u[id*3+jj]*e[ii][jj];
//				for(kk=0; kk<3; kk++){
//					E2 += E[jj][kk] * e[ii][jj] * e[ii][kk];
//				}
			}
			udote_sq = udote * udote;
			E1       =   8.0 * E2;
			if(compressibility!=0){
				rho0 = Rho[id];
			}
			if(ii==0){
				fout[id*15+ii] = fin[id*15+ii] - (A0 + rho0 * (           C0*usq                   ));
			} else if(ii<7) {
				fout[id*15+ii] = fin[id*15+ii] - (A1 + rho0 * (B1*udote + C1*usq + D1*udote_sq + E1));
			} else {
				fout[id*15+ii] = fin[id*15+ii] - (A2 + rho0 * (B2*udote + C2*usq + D2*udote_sq + E2));
			}
			
		}
	}

}

void cal_p()
{
	int i, ii=0;
	double temp, T1, T2;
	
	T2 = 1.0/24.0;
	T1 = 8*T2;
	
	for (i=0; i<points; i++) {
		for (ii=0; ii<15; ii++) {
			temp = yforce * e[ii][1];
			if(ii>0 && ii<7){
				p[i*15+ii] = T1 * temp;
			}else if(ii>6){
				p[i*15+ii] = T2 * temp;
			}else {
				p[i*15+ii] = 0 * temp;
			}
			
			if (Q_on!=0) {
				if(ii>0 && ii<7){
					p[i*15+ii] += T1 * (sigma_p[i*3]*e[ii][0] +sigma_p[i*3+1]*e[ii][1] +sigma_p[i*3+2]*e[ii][2] ); 
				}
				else if (ii>6){
					p[i*15+ii] += T2 * (sigma_p[i*3]*e[ii][0] +sigma_p[i*3+1]*e[ii][1] +sigma_p[i*3+2]*e[ii][2] );
				}
			}
			
		}
	}
	
}