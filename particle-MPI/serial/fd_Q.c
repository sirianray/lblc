/*
 *  fd_Q.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */
#include "main.h"

void evol_Q()
{
	int id, id0=0, iq, ii, jj, kk, ip, iw=-9;
	double W1[3][3]={0.}, M[3][3]={0.}, S[3][3]={0.}, trQW=0.0, temp;
	
	if (t_current%t_print==0) Q_diff=0;
	
	for (id=0, iq=0; id<qpoints; id++, iq+=5) {
		if (info[id]==-1) {
			if (flow_on!=0) {
				iw       = 9*(id-bulk0);
				W1[0][0] = xi * W[iw];
				W1[0][1] = xi1* W[iw+3] + xi2* W[iw+1];
				W1[1][0] = xi1* W[iw+1] + xi2* W[iw+3];
				W1[0][2] = xi1* W[iw+6] + xi2* W[iw+2];
				W1[2][0] = xi1* W[iw+2] + xi2* W[iw+6];
				W1[1][1] = xi * W[iw+4];
				W1[1][2] = xi1* W[iw+7] + xi2* W[iw+5];
				W1[2][1] = xi1* W[iw+5] + xi2* W[iw+7];
				W1[2][2] = xi * W[iw+8];
				
				M[0][0]  = Q[iq]+third;
				M[0][1]  = Q[iq+1];
				M[1][0]  = M[0][1];
				M[0][2]  = Q[iq+2];
				M[2][0]  = M[0][2];
				M[1][1]  = Q[iq+3]+third;
				M[1][2]  = Q[iq+4];
				M[2][1]  = M[1][2];
				M[2][2]  = 1.0-M[0][0]-M[1][1];
				
				trQW = Q[iq]*(W[iw]-W[iw+8])+Q[iq+3]*(W[iw+4]-W[iw+8])+Q[iq+1]*(W[iw+1]+W[iw+3])+Q[iq+2]*(W[iw+2]+W[iw+6])+Q[iq+4]*(W[iw+5]+W[iw+7]);
				
				for(ii=0;ii<2;ii++){
					for(jj=ii;jj<3;jj++){
						temp = 0;
						for(kk=0;kk<3;kk++){
							temp += W1[ii][kk]*M[kk][jj] + M[ii][kk]*W1[jj][kk];
						}
						S[ii][jj] = temp - 2.0*xi*M[ii][jj]*trQW;
					}
				}
			}
			
			Q[iq]  += qdt * ( H[iq]   + S[0][0] );
			Q[iq+1]+= qdt * ( H[iq+1] + S[0][1] );
			Q[iq+2]+= qdt * ( H[iq+2] + S[0][2] );
			Q[iq+3]+= qdt * ( H[iq+3] + S[1][1] );
			Q[iq+4]+= qdt * ( H[iq+4] + S[1][2] );
//			Q[iq+5] = -Q[iq] - Q[iq+3];
//			printf("%d\n",iq);
			if (t_current%t_print==0) {
				Q_diff += (H[iq]+S[0][0])*(H[iq]+S[0][0])+(H[iq+1]+S[0][1])*(H[iq+1]+S[0][1])+(H[iq+2]+S[0][2])*(H[iq+2]+S[0][2])+(H[iq+3]+S[1][1])*(H[iq+3]+S[1][1])+(H[iq+4]+S[1][2])*(H[iq+4]+S[1][2])+(H[iq]+S[0][0])*(H[iq+3]+S[1][1]);
			}
		} else {
			ip = info[id];
			if (ip>=0 && surf[ip]>0) {
				for (ii=0; ii<5; ii++) {
					Q[iq+ii] += -qdt*Gamma_rot*(-kappa*H[iq+ii]+surf[ip+1]*(Q[iq+ii]-surf[ip+5+ii]));
				}
//				Q[iq+5] = -Q[iq] - Q[iq+3];
			}
		}
	}
}

void cal_dQ()
{
	int id, iq, ib, iu, ii, nid[6], flag=0, ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, eld, eel;
	real *p=NULL;
	
	if (t_current%t_print==0) {
		flag  = 1;
		e_tot = 0;
		e_ld  = 0;
		e_el  = 0;
		e_sf  = 0;
	}
	
	for (iq=0, ib=0, id=0; id<qpoints; iq+=5, ib+=6, id++){
		if (info[id]==-1) {
			p  = &Q[iq];
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighbor[ib+ii];
			}
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2xQ[ii]= Q[ip2+ii] - 2.0 * Q[iq+ii] + Q[ip1+ii];
				}
			} else if (ip1<0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2xQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Q[-ip1-1+ii];
				}
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = four3rd*Q[-ip2-1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2xQ[ii]= eight3rd*Q[-ip2-1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2yQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1<0 && ip2>=0) {
//				printf("%d",-ip1-1+ii);
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2yQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Q[-ip1-1+ii];
//					printf(" %e",Q[-ip1-1+ii]);
				}
//				printf("\n");
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = four3rd*Q[-ip2-1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2yQ[ii]= eight3rd*Q[-ip2-1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2zQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1<0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2zQ[ii]= four3rd*(Q[ip2+ii]-*(p+ii))+eight3rd*(Q[-ip1-1+ii]-*(p+ii));					
				}
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] =-(third*Q[ip1+ii]+*(p+ii)-four3rd*Q[-ip2-1+ii]);
					d2zQ[ii]= four3rd*(Q[ip1+ii]-*(p+ii))+eight3rd*(Q[-ip2-1+ii]-*(p+ii));
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			trqq = trQQ(p);
			
			for (ii=0; ii<5; ii++) {
				H[iq+ii] = -1*(A_ldg*(1.0-third*U_lc)*(*(p+ii)) - A_ldg*U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - kappa * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
//				printf("%e ",d2xQ[ii]+d2yQ[ii]+d2zQ[ii]);
			}
//			printf("\n");
			
			if (flag==1) {
				qqq = QQQ(p);
				eld = 0.5 * (1.0-third*U_lc)*trqq - third*U_lc*qqq + 0.25*U_lc*trqq*trqq;
				eel = dxQ[0]*dxQ[0] + dyQ[0]*dyQ[0] + dzQ[0]*dzQ[0];				
				for (ii=1; ii<5; ii++) {
					eel += dxQ[ii]*dxQ[ii] + dyQ[ii]*dyQ[ii] + dzQ[ii]*dzQ[ii];
				}
				eel += dxQ[0]*dxQ[3] + dyQ[0]*dyQ[3] + dzQ[0]*dzQ[3];
				eel *= kappa;
				e_ld+= eld;
				e_el+= eel;
			}
//			printf("%e %e %e %e %e\n",H[iq],H[iq+1],H[iq+2],H[iq+3],H[iq+4]);
//			combine H and convective of Q			
			iu = (id - bulk0) * 3;
			if (flow_on!=0) {
				for (ii=0; ii<5; ii++) {
					H[iq+ii] = Gamma_rot * H[iq+ii] - u[iu]*dxQ[ii] - u[iu+1]*dyQ[ii] - u[iu+2]*dzQ[ii];
				}
			} else {
				for (ii=0; ii<5; ii++) {
					H[iq+ii] = Gamma_rot * H[iq+ii] ;
				}
			}
//			printf("%e %e %e %e\n",Gamma_rot,u[iu],u[iu+1],u[iu+2]);
		} else {
			ip = info[id];
			if (ip>=0 && surf[ip]>0) {
				if (surf[ip+2]!=0) {
					ip1 = neighbor[ib];
					ip2 = neighbor[ib+1];
					if (ip1>=0) {
						for (ii=0; ii<5; ii++) {
							dxQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Q[iq+ii];
						}
					} else if (ip2>0) {
						for (ii=0; ii<5; ii++) {
							dxQ[ii] = 0.5*(Q[ip2+ii] - Q[-ip1-1+ii]);
						}
					} else {
						for (ii=0; ii<5; ii++) {
							dxQ[ii] = -0.5*Q[-ip2-1+ii] + 2.0*Q[-ip1-1+ii] - 1.5*Q[iq+ii];
						}
					}
				}
				if (surf[ip+3]!=0) {
					ip1 = neighbor[ib+2];
					ip2 = neighbor[ib+3];
					if (ip1>=0) {
						for (ii=0; ii<5; ii++) {
							dyQ[ii] =-third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Q[iq+ii];
						}
					} else if (ip2>0) {
						for (ii=0; ii<5; ii++) {
							dyQ[ii] = 0.5*(Q[ip2+ii] - Q[-ip1-1+ii]);
						}
					} else {
						for (ii=0; ii<5; ii++) {
							dyQ[ii] =-0.5*Q[-ip2-1+ii] + 2.0*Q[-ip1-1+ii] - 1.5*Q[iq+ii];
						}
					}					
				}
				if (surf[ip+4]!=0) {
					ip1 = neighbor[ib+4];
					ip2 = neighbor[ib+5];
					if (ip1>=0) {
						for (ii=0; ii<5; ii++) {
							dzQ[ii] =-third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Q[iq+ii];
						}
					} else if (ip2>0) {
						for (ii=0; ii<5; ii++) {
							dzQ[ii] = 0.5*(Q[ip2+ii] - Q[-ip1-1+ii]);
						}
					} else {
						for (ii=0; ii<5; ii++) {
							dzQ[ii] =-0.5*Q[-ip2-1+ii] + 2.0*Q[-ip1-1+ii] - 1.5*Q[iq+ii];
						}
					}
				}
				for (ii=0; ii<5; ii++) {
					H[iq+ii] = dxQ[ii]*fabs(surf[ip+2])+dyQ[ii]*fabs(surf[ip+3])+dzQ[ii]*fabs(surf[ip+4]);
				}
				if (flag==1) {
					for (ii=0; ii<4; ii++) {
						e_sf += surf[ip+1]*(Q[iq+ii]-surf[ip+5+ii])*(Q[iq+ii]-surf[ip+5+ii]);
					}
					e_sf += surf[ip+1]*(Q[iq]-surf[ip+5])*(Q[iq+3]-surf[ip+5+3]);
				}
			}
		}
	}
	e_tot = e_ld + e_el + e_sf;
}


void cal_stress()
{
	int id, iq, ib, ii, jj, mm, nid[6], ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, h[3][3], M[3][3], dQ[3][3], dQ2, temp, sum;
	real *p=NULL;
	
	for (id=0, iq=bulk0*5, ib=bulk0*6; id<points; id++, iq+=5, ib+=6){
		if (info[id+bulk0]==-1) {
			p  = &Q[iq];
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighbor[ib+ii];
			}
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2xQ[ii]= Q[ip2+ii] - 2.0 * Q[iq+ii] + Q[ip1+ii];
				}				
			} else if (ip1<0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2xQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Q[-ip1-1+ii];
				}
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = four3rd*Q[-ip2-1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2xQ[ii]= eight3rd*Q[-ip2-1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2yQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1<0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2yQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Q[-ip1-1+ii];
				}
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = four3rd*Q[-ip2-1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2yQ[ii]= eight3rd*Q[-ip2-1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1>=0 && ip2>=0) {
				for (ii=0; ii<6; ii++) {
					dzQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2zQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1<0 && ip2>=0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Q[-ip1-1+ii];
					d2zQ[ii]= four3rd*(Q[ip2+ii]-*(p+ii))+eight3rd*(Q[-ip1-1+ii]-*(p+ii));
				}
			} else if (ip1>=0 && ip2<0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] =-(*(p+ii)+third*Q[ip1+ii]-four3rd*Q[-ip2-1+ii]);
					d2zQ[ii]= eight3rd*(Q[-ip2-1+ii]-*(p+ii))+four3rd*(Q[ip1+ii]-*(p+ii));
				}
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			trqq = trQQ(p);
			
			for (ii=0; ii<5; ii++) {
				H[iq+ii] = -1*(A_ldg*(1.0-third*U_lc)*(*(p+ii)) - A_ldg*U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - kappa * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
			}
			
			dQ[0][0]=2.*(dxQ[0]*dxQ[0]+dxQ[1]*dxQ[1]+dxQ[2]*dxQ[2]+dxQ[3]*dxQ[3]+dxQ[4]*dxQ[4]+dxQ[0]*dxQ[3]);
			dQ[1][1]=2.*(dyQ[0]*dyQ[0]+dyQ[1]*dyQ[1]+dyQ[2]*dyQ[2]+dyQ[3]*dyQ[3]+dyQ[4]*dyQ[4]+dyQ[0]*dyQ[3]);
			dQ[2][2]=2.*(dzQ[0]*dzQ[0]+dzQ[1]*dzQ[1]+dzQ[2]*dzQ[2]+dzQ[3]*dzQ[3]+dzQ[4]*dzQ[4]+dzQ[0]*dzQ[3]);
			dQ[0][1]=2.*(dxQ[0]*dyQ[0]+dxQ[1]*dyQ[1]+dxQ[2]*dyQ[2]+dxQ[3]*dyQ[3]+dxQ[4]*dyQ[4]+dxQ[0]*dyQ[3]);
			dQ[0][2]=2.*(dxQ[0]*dzQ[0]+dxQ[1]*dzQ[1]+dxQ[2]*dzQ[2]+dxQ[3]*dzQ[3]+dxQ[4]*dzQ[4]+dxQ[0]*dzQ[3]);
			dQ[1][2]=2.*(dyQ[0]*dzQ[0]+dyQ[1]*dzQ[1]+dyQ[2]*dzQ[2]+dyQ[3]*dzQ[3]+dyQ[4]*dzQ[4]+dyQ[0]*dzQ[3]);
			dQ[1][0]=dQ[0][1];
			dQ[2][0]=dQ[0][2];
			dQ[2][1]=dQ[1][2];
			
			dQ2 = dQ[0][0]+dQ[1][1]+dQ[2][2];
			
			M[0][0] = Q[iq]+third;
			M[0][1] = Q[iq+1];
			M[1][0] = M[0][1];
			M[0][2] = Q[iq+2];
			M[2][0] = M[0][2];
			M[1][1] = Q[iq+3]+third;
			M[1][2] = Q[iq+4];
			M[2][1] = M[1][2];
			M[2][2] = 1.0-M[0][0]-M[1][1];
			
			h[0][0] = H[iq];
			h[0][1] = H[iq+1];
			h[1][0] = h[0][1];
			h[0][2] = H[iq+2];
			h[2][0] = h[0][2];
			h[1][1] = H[iq+3];
			h[1][2] = H[iq+4];
			h[2][1] = h[1][2];
			h[2][2] =-H[iq]-H[iq+3];
			
			sum = 2.0*(Q[iq]*H[iq]+Q[iq+3]*H[iq+3]+Q[iq+1]*H[iq+1]+Q[iq+2]*H[iq+2]+Q[iq+4]*H[iq+4])+Q[iq]*H[iq+3]+Q[iq+3]*H[iq];
			
			for(ii=0;ii<3;ii++){
				for(jj=0;jj<3;jj++){
					temp = 0.5*kappa*dQ2*Delta(ii,jj) + 2*xi*M[ii][jj]*sum;
					
					for(mm=0;mm<3;mm++){
						temp += (1.0-xi)*M[ii][mm]*h[mm][jj]-(xi+1.0)*h[ii][mm]*M[mm][jj];
					}
					
					temp += -kappa * dQ[ii][jj];
					sigma_q[id*9+ii*3+jj] = temp;
					
				}
			}
		}
	}
	
}


void cal_sigma_p()
{
	int id0, id, ip, is, idxm, idxp, idym, idyp, idzm, idzp, idzmm, idzpp;
	
	for (id0=0; id0<points; id0++) {
		is = id0 * 3;
		id = is  * 3;
		ip = is  * 5;
		
		if (info[bulk0+id0]==-1) {
			idxp = (int)(nextf[ip+2]/15)*9;
			idxm = (int)(nextf[ip+1]/15)*9;
			idyp = (int)(nextf[ip+4]/15)*9;
			idym = (int)(nextf[ip+3]/15)*9;
			idzp = (int)(nextf[ip+6]/15)*9;
			idzm = (int)(nextf[ip+5]/15)*9;
			idzpp= (int)(nextf[nextf[ip+6]]/15)*9;
			idzmm= (int)(nextf[nextf[ip+5]]/15)*9;
			
			if (npar==0 || idxm!=id && idxp!=id) {
				sigma_p[is]  = 0.5 * (sigma_q[idxp]  - sigma_q[idxm]);
				sigma_p[is+1]= 0.5 * (sigma_q[idxp+3]- sigma_q[idxm+3]);
				sigma_p[is+2]= 0.5 * (sigma_q[idxp+6]- sigma_q[idxm+6]);
			} else if (idxm==id) {
				sigma_p[is]  = sigma_q[idxp]  - sigma_q[id];
				sigma_p[is+1]= sigma_q[idxp+3]- sigma_q[id+3];
				sigma_p[is+2]= sigma_q[idxp+6]- sigma_q[id+6];
			} else if (idxp==id) {
				sigma_p[is]  = sigma_q[id]  - sigma_q[idxm];
				sigma_p[is+1]= sigma_q[id+3]- sigma_q[idxm+3];
				sigma_p[is+2]= sigma_q[id+6]- sigma_q[idxm+6];
			}
			
			if (npar==0 || idym!=id && idyp!=id) {
				sigma_p[is]  += 0.5 * (sigma_q[idyp+1]- sigma_q[idym+1]);
				sigma_p[is+1]+= 0.5 * (sigma_q[idyp+4]- sigma_q[idym+4]);
				sigma_p[is+2]+= 0.5 * (sigma_q[idyp+7]- sigma_q[idym+7]);
			} else if (idym==id) {
				sigma_p[is]  += sigma_q[idyp+1]- sigma_q[id+1];
				sigma_p[is+1]+= sigma_q[idyp+4]- sigma_q[id+4];
				sigma_p[is+2]+= sigma_q[idyp+7]- sigma_q[id+7];
			} else if (idyp==id) {
				sigma_p[is]  += sigma_q[id+1]- sigma_q[idym+1];
				sigma_p[is+1]+= sigma_q[id+4]- sigma_q[idym+4];
				sigma_p[is+2]+= sigma_q[id+7]- sigma_q[idym+7];
			}
			
			if (idzp!=id && idzm!=id) {
				sigma_p[is]  += 0.5 * (sigma_q[idzp+2] - sigma_q[idzm+2]);
				sigma_p[is+1]+= 0.5 * (sigma_q[idzp+5] - sigma_q[idzm+5]);
				sigma_p[is+2]+= 0.5 * (sigma_q[idzp+8] - sigma_q[idzm+8]);
			} else if (wall_on!=0 && id0/bulk0==Nz-1) {
				sigma_p[is]  +=  0.5*sigma_q[idzmm+2] + 1.5*sigma_q[id+2] - 2.0*sigma_q[idzm+2];
				sigma_p[is+1]+=  0.5*sigma_q[idzmm+5] + 1.5*sigma_q[id+5] - 2.0*sigma_q[idzm+5];
				sigma_p[is+2]+=  0.5*sigma_q[idzmm+8] + 1.5*sigma_q[id+8] - 2.0*sigma_q[idzm+8];
			} else if (wall_on!=0 && id0/bulk0==0) {
				sigma_p[is]  +=-(0.5*sigma_q[idzpp+2] + 1.5*sigma_q[id+2] - 2.0*sigma_q[idzp+2]);
				sigma_p[is+1]+=-(0.5*sigma_q[idzpp+5] + 1.5*sigma_q[id+5] - 2.0*sigma_q[idzp+5]);
				sigma_p[is+2]+=-(0.5*sigma_q[idzpp+8] + 1.5*sigma_q[id+8] - 2.0*sigma_q[idzp+8]);
			} else if (npar>0 && idzp==id) {
				sigma_p[is]  += sigma_q[id+2]- sigma_q[idzm+2];
				sigma_p[is+1]+= sigma_q[id+5]- sigma_q[idzm+5];
				sigma_p[is+2]+= sigma_q[id+8]- sigma_q[idzm+8];
			} else if (npar>0 && idzm==id) {
				sigma_p[is]  += sigma_q[idzp+2]- sigma_q[id+2];
				sigma_p[is+1]+= sigma_q[idzp+5]- sigma_q[id+5];
				sigma_p[is+2]+= sigma_q[idzp+8]- sigma_q[id+8];
			}
		} else {
			sigma_p[is]  = 0;
			sigma_p[is+1]= 0;
			sigma_p[is+2]= 0;
		}
	}
}
