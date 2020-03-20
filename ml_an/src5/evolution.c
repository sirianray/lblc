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
	
	cal_p();
	streaming(fin,fout);

	cal_feqf_diff(fin,Cf);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winf);
        MPI_Win_fence(0, winf2);
	MPI_Win_fence(0, winp);
	for (id=0; id<lpoint*15; id++) {
		Cf1      = -itau_f * Cf[id] + p[id];
		ip       = nextf[id];
		fin[id] += 0.5 * dt * Cf1;
		fout[ip]+= dt * Cf1;
	}
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winf);
        MPI_Win_fence(0, winf2);

//	cal_rho(fout);
	cal_u(fout);
	cal_p();
	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winu);
        MPI_Win_fence(0, winr);
	cal_feqf_diff(fout,fout);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winf);
        MPI_Win_fence(0, winf2);
	MPI_Win_fence(0, winp);

	for (id=0; id<lpoint*15; id++) {
		ip       = nextf[id];
		Cf2      = -itau_f * fout[ip] + p[ip];
		fin[id] += dt * 0.5 * Cf2;
	}
	
	streaming(fin,fout);

	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winf);
        MPI_Win_fence(0, winf2);	
//	cal_rho(fout);
	cal_u(fout);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winu);
        MPI_Win_fence(0, winr);
}

void streaming(real *fin, real *fout)
{
	int ip, id, idn, i, in, j, jn, k, kn, ii;
	double Bx2;
	
	for (ip=0; ip<lpoint; ip++) {
		for (ii=0; ii<15; ii++) {
			id = ip*15+ii;
			idn= nextf[id];
			fout[idn]=fin[id];
			if(ii<7) {
				Bx2 = two3rd;
			} else {
				Bx2 = one12th;
			}
                        if (wall_x!=0) {
                                i = (ip+myid*point)%Nx;
                                in = i + e[ii][0];
                                if ( in==-1 || in==Nx ) {
                                        if (i==0) {
                                                fout[idn] += -Bx2*Rho[ip]*(uy_lo*e[ii][1]+uz_lo*e[ii][2]);
                                        } else {
                                                fout[idn] += -Bx2*Rho[ip]*(uy_hi*e[ii][1]+uz_hi*e[ii][2]);
                                        }
                                }
                        }
                        if (wall_y!=0) {
                                j = ((ip+myid*point)/Nx)%Ny;
                                jn = j + e[ii][1];
                                if ( jn==-1 || jn==Ny ) {
                                        if (j==0) {
                                                fout[idn] += -Bx2*Rho[ip]*(ux_lo*e[ii][0]+uz_lo*e[ii][2]);
                                        } else {
                                                fout[idn] += -Bx2*Rho[ip]*(ux_hi*e[ii][0]+uz_hi*e[ii][2]);
                                        }
                                }
                        }
			if (wall_z!=0) {
				k = (ip+myid*point)/bulk0;
				kn = k + e[ii][2];
				if ( kn==-1 || kn==Nz ) {
					if (k==0) {
						fout[idn] += -Bx2*Rho[ip]*(ux_lo*e[ii][0]+uy_lo*e[ii][1]);
					} else {
						fout[idn] += -Bx2*Rho[ip]*(ux_hi*e[ii][0]+uy_hi*e[ii][1]);
					}
				}
			}
			if (npar>0 && abs(id-idn)<15) {
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
	
	for (id=0; id<lpoint*3; id+=3) {
		ip   = 5*id;
		idxp = nextf[ip+2]>=0?(int)(nextf[ip+2]/15)*3:(int)(nextf[ip+2]/15-1)*3;
		idxm = nextf[ip+1]>=0?(int)(nextf[ip+1]/15)*3:(int)(nextf[ip+1]/15-1)*3;
		idyp = nextf[ip+4]>=0?(int)(nextf[ip+4]/15)*3:(int)(nextf[ip+4]/15-1)*3;
		idym = nextf[ip+3]>=0?(int)(nextf[ip+3]/15)*3:(int)(nextf[ip+3]/15-1)*3;
		idzp = nextf[ip+6]>=0?(int)(nextf[ip+6]/15)*3:(int)(nextf[ip+6]/15-1)*3;
		idzm = nextf[ip+5]>=0?(int)(nextf[ip+5]/15)*3:(int)(nextf[ip+5]/15-1)*3;
		
		iw   = 3*id;
		
		if (idxm!=id && idxp!=id) {
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
		
		if (idym!=id && idyp!=id) {
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
		} else if (wall_z!=0 && (id+myid*point*3)/bulk3==0) {
			W[iw+6]=  u[id]   + third*u[idzp]   - four3rd * ux_lo;
			W[iw+7]=  u[id+1] + third*u[idzp+1] - four3rd * uy_lo;
			W[iw+8]=  u[id+2] + third*u[idzp+2];
		} else if (wall_z!=0 && (id+myid*point*3)/bulk3==Nz-1) {
			W[iw+6]= -u[id]   - third*u[idzm]   + four3rd * ux_hi;
			W[iw+7]= -u[id+1] - third*u[idzm+1] + four3rd * uy_hi;
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
	
	for (i=0; i<lpoint; i++) {
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
	
	for (i=0; i<lpoint; i++) {
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
	
	C2 = -one24th*0.5;
	C1 = -one12th*2.0;
	C0 = -third;
	
	D2 = 0.0625;
	D1 = 0.5;
	D0 = 0.0;
	

	for (id=0; id<lpoint; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
                A2 = -one24th*trsigma*third;
                A1 = -third*trsigma*third;
                A0 =  Rho[id] - 8.0 * A2 - 6.0 * A1;
		
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

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winf);
        MPI_Win_fence(0, winf2);
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
	
	C2 = -one24th*0.5;
	C1 = -one12th*2.0;
	C0 = -third;
	
	D2 = 0.0625;
	D1 = 0.5;
//	D0 = 0.0;
	
	for (id=0; id<lpoint; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
                A2 = -one24th*trsigma*third;
                A1 = -third*trsigma*third;
                A0 =  Rho[id] - 8.0 * A2 - 6.0 * A1;
		
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
	
	C2 = -one24th*0.5;
	C1 = -one12th*2.0;
	C0 = -third;
	
	D2 = 0.0625;
	D1 = 0.5;
//	D0 = 0.0;
	
	for (id=0; id<lpoint; id++) {
		for (ii=0; ii<3; ii++) sigma[ii][ii] = -Rho[id]*Temp;
		trsigma = -Rho[id];
                A2 = -one24th*trsigma*third;
                A1 = -third*trsigma*third;
                A0 =  Rho[id] - 8.0 * A2 - 6.0 * A1;
		
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
	double temp, T1, T2, T0, edotu;
	
        T2 = one24th;
        T1 = third;
        T0 = two3rd;
	
	for (i=0; i<lpoint; i++) {
		for (ii=0; ii<15; ii++) {
                        edotu = e[ii][0]*u[i*3]+e[ii][1]*u[i*3+1]+e[ii][2]*u[i*3+2];
//                        p[i*15+ii] = ( (e[ii][1] - u[i*3+1]) + 3.0 * edotu * e[ii][1] )*yforce;
			p[i*15+ii] = (e[ii][0] - u[i*3])*(xforce-fr*u[i*3]) + (e[ii][1] - u[i*3+1])*(yforce-fr*u[i*3+1]) + (e[ii][2] - u[i*3+2])*0. + 3.0 * edotu * (e[ii][0]*(xforce-fr*u[i*3]) + e[ii][1]*(yforce-fr*u[i*3+1]) + e[ii][2]*0.);

                        if (Q_on!=0) {
                                p[i*15+ii]+= ( (e[ii][0]-u[i*3])*sigma_p[i*3] + (e[ii][1]-u[i*3+1])*(sigma_p[i*3+1]) +  (e[ii][2]-u[i*3+2])*0. );
                                p[i*15+ii]+= 3.0 * edotu *(e[ii][0]*sigma_p[i*3] + e[ii][1]*sigma_p[i*3+1] + e[ii][2]*0.);
                        }

                        if (ii>6) {
                                p[i*15+ii] = T2*p[i*15+ii];
                        } else if(ii>0) {
                                p[i*15+ii] = T1*p[i*15+ii];
                        } else {
                                p[i*15+ii] = T0*p[i*15+ii];
                        }	
		}
	}
	
}
