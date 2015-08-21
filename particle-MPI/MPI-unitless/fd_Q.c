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
	double W1[3][3]={0.}, M[3][3]={0.}, S[3][3]={0.}, trQW=0.0, temp, Q_diffi;
	int nn, flag=0;
	double nu[3], Qs[3][3], Qtilde[3][3], Pmat[3][3], nuQnu, Qperp[3][3], Qtildevec[5], Qperpvec[5], e_sfi=0, e_sfn;
	
	if (t_current%t_print==0) {
		Q_diffi= 0;
		flag   = 1;
	}
	
	for (id=0, iq=0; id<lpoint; id++, iq+=5) {
		if (info[id]==-1) {
			if (flow_on!=0) {
				iw       = 9*id;
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
			if (flag==1) {
				Q_diffi += (H[iq]+S[0][0])*(H[iq]+S[0][0])+(H[iq+1]+S[0][1])*(H[iq+1]+S[0][1])+(H[iq+2]+S[0][2])*(H[iq+2]+S[0][2])+(H[iq+3]+S[1][1])*(H[iq+3]+S[1][1])+(H[iq+4]+S[1][2])*(H[iq+4]+S[1][2])+(H[iq]+S[0][0])*(H[iq+3]+S[1][1]);
			}
		}
	}
	MPI_Win_fence(0,winq);
	MPI_Reduce(&Q_diffi, &Q_diff, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

	if (wall_on!=0 || npar>0) {
		for (id=0; id<node; id++) {
			if (id+myid*node<nsurf){
				iq = id * 5;
				ip = iq * 2;
				if (surf[ip]==1) {
					for (ii=0; ii<5; ii++) {
						Qsurf[iq+ii] += -qdt*Gamma_rot*(-Hsurf[iq+ii]/l0+surf[ip+1]*(Qsurf[iq+ii]-surf[ip+5+ii]));
					}
				} else if (surf[ip]==2) {
					for(ii=0; ii<3; ii++) nu[ii] = surf[ip+2+ii];
					for(kk=0, ii=0; ii<2; ii++) {
						for(jj=ii; jj<3; jj++){
							Qs[ii][jj]=Qsurf[iq+kk];
							kk++;
						}
					}
					Qs[1][0] = Qs[0][1];
					Qs[2][0] = Qs[0][2];
					Qs[2][1] = Qs[1][2];
					Qs[2][2] =-Qs[0][0]-Qs[1][1];

					for(ii=0;ii<3;ii++) {
						for(jj=0;jj<3;jj++) {
							Qtilde[ii][jj] = Qs[ii][jj] + One3rdDelta(ii,jj)*S_lc;
							Pmat[ii][jj] = Delta(ii,jj) - nu[ii]*nu[jj];
						}
					}	

					nuQnu = 0;
					for(ii=0;ii<3;ii++) {
						for(jj=0;jj<3;jj++) {
							nuQnu += nu[ii]*Qtilde[ii][jj]*nu[jj];
						}
					}

					for(ii=0;ii<3;ii++) {
						for(jj=0;jj<3;jj++) {
							Qperp[ii][jj]=0;
						}
					}

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								for(nn=0;nn<3;nn++){
									Qperp[ii][jj] += Pmat[ii][kk]*Qtilde[kk][nn]*Pmat[nn][jj];
								}
							}
						}
					}

					kk=0;
					for(ii=0;ii<2;ii++){
						for(jj=ii;jj<3;jj++){
							 Qtildevec[kk] = Qtilde[ii][jj];
							 Qperpvec[kk] = Qperp[ii][jj];
							 ++kk;
						}
					}

					for(ii=0;ii<5;ii++){
						if(ii==0 || ii==3) temp = third*nuQnu;
						else temp = 0;
						Qsurf[iq+ii] += - qdt*Gamma_rot * (-Hsurf[iq+ii]/l0 + 2.0*surf[ip+1]*( (Qtildevec[ii]-Qperpvec[ii]) - temp));
						if(flag==1)e_sfi+= surf[ip+1]*(Qtildevec[ii]-Qperpvec[ii])*(Qtildevec[ii]-Qperpvec[ii]);
					}
					if(flag==1)e_sfi += surf[ip+1]*(Qtildevec[0]-Qperpvec[0])*(Qtildevec[3]-Qperpvec[3]);
				}
			}
		}
		MPI_Win_fence(0,winQsurf);
		if (flag==1) {
			MPI_Reduce(&e_sfi, &e_sfn, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
			if (myid==root) {
				e_sf += e_sfn;
				e_tot+= e_sfn;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void cal_dQ()
{
	int id, iq, ib, iu, ii, nid[6], flag=0, ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, eld, eel;
	double e_ldi, e_eli, e_sfi;
	real *p=NULL;
	
	if (t_current%t_print==0) {
		flag  = 1;
		e_tot = 0;
		e_ldi = 0;
		e_eli = 0;
		e_sfi = 0;
	}
	
	for (iq=0, ib=0, id=0; id<lpoint; iq+=5, ib+=6, id++){
		if (info[id]==-1) {
			p  = &Q[iq];
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[ib+ii];
			}
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2xQ[ii]= Q[ip2+ii] - 2.0 * Q[iq+ii] + Q[ip1+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2xQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Qsurf[ip1+1+ii];
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2xQ[ii]= eight3rd*Qsurf[ip2+1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dxQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2xQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2yQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2yQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Qsurf[ip1+1+ii];
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2yQ[ii]= eight3rd*Qsurf[ip2+1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dyQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2yQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2zQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2zQ[ii]= four3rd*(Q[ip2+ii]-*(p+ii))+eight3rd*(Qsurf[ip1+1+ii]-*(p+ii));					
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] =-(third*Q[ip1+ii]+*(p+ii)-four3rd*Qsurf[ip2+1+ii]);
					d2zQ[ii]= four3rd*(Q[ip1+ii]-*(p+ii))+eight3rd*(Qsurf[ip2+1+ii]-*(p+ii));
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dzQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2zQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			trqq = trQQ(p);
			
			for (ii=0; ii<5; ii++) {
				H[iq+ii] = -1*((1.0-third*U_lc)*(*(p+ii)) - U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - kappa * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
			}
			
			if (flag==1) {
				qqq = QQQ(p);
				eld = 0.5 * (1.0-third*U_lc)*trqq - third*U_lc*qqq + 0.25*U_lc*trqq*trqq;
				eel = dxQ[0]*dxQ[0] + dyQ[0]*dyQ[0] + dzQ[0]*dzQ[0];				
				for (ii=1; ii<5; ii++) {
					eel += dxQ[ii]*dxQ[ii] + dyQ[ii]*dyQ[ii] + dzQ[ii]*dzQ[ii];
				}
				eel  += dxQ[0]*dxQ[3] + dyQ[0]*dyQ[3] + dzQ[0]*dzQ[3];
				eel  *= kappa;
				e_ldi+= eld;
				e_eli+= eel;
			}
//			combine H and convective of Q
			iu = id * 3;
			if (flow_on!=0) {
				for (ii=0; ii<5; ii++) {
					H[iq+ii] = Gamma_rot * H[iq+ii] - u[iu]*dxQ[ii] - u[iu+1]*dyQ[ii] - u[iu+2]*dzQ[ii];
				}
			} else {
				for (ii=0; ii<5; ii++) {
					H[iq+ii] = Gamma_rot * H[iq+ii] ;
				}
			}
		}
	}

	if (wall_on!=0 || npar>0) {
		for (id=0; id<node; id++) {
			if (id+myid*node<nsurf) {
				iq = id * 5;
				ip = iq * 2;
				if (surf[ip]>0) {
					ib = id * 6;
					if (surf[ip+2]!=0) {
						ip1 = neighbsurf[ib];
						ip2 = neighbsurf[ib+1];
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dxQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
							}
						} else if (ip1%5==4 || ip1%5==-1) {
							for (ii=0; ii<5; ii++) {
								dxQ[ii] = 0.5*(Q[ip2+ii] - Q[ip1+1+ii]);
							}
						} else if (ip1%5==3 || ip1%5==-2) {
							for (ii=0; ii<5; ii++) {
								dxQ[ii] = -0.5*Q[ip2+ii] + 2.0*Q[ip1+2+ii] - 1.5*Qsurf[iq+ii];
							}
						}
					}

					if (surf[ip+3]!=0) {
						ip1 = neighbsurf[ib+2];
						ip2 = neighbsurf[ib+3];
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dyQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
							}
						} else if (ip1%5==4 || ip1%5==-1) {
							for (ii=0; ii<5; ii++) {
								dyQ[ii] = 0.5*(Q[ip2+ii] - Q[ip1+1+ii]);
							}
						} else if (ip1%5==3 || ip1%5==-2) {
							for (ii=0; ii<5; ii++) {
								dyQ[ii] = -0.5*Q[ip2+ii] + 2.0*Q[ip1+2+ii] - 1.5*Qsurf[iq+ii];
							}
						}
					}

					if (surf[ip+4]!=0) {
						ip1 = neighbsurf[ib+4];
						ip2 = neighbsurf[ib+5];
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dzQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
							}
						} else if (ip1%5==4 || ip1%5==-1) {
							for (ii=0; ii<5; ii++) {
								dzQ[ii] = 0.5*(Q[ip2+ii] - Q[ip1+1+ii]);
							}
						} else if (ip1%5==3 || ip1%5==-2) {
							for (ii=0; ii<5; ii++) {
								dzQ[ii] = -0.5*Q[ip2+ii] + 2.0*Q[ip1+2+ii] - 1.5*Qsurf[iq+ii];
							}
						}
					}
				}

				for (ii=0; ii<5; ii++) {
					Hsurf[iq+ii] = dxQ[ii]*fabs(surf[ip+2])+dyQ[ii]*fabs(surf[ip+3])+dzQ[ii]*fabs(surf[ip+4]);
				}
				if (flag==1 && surf[ip]==1) {
					for (ii=0; ii<5; ii++) {
						e_sfi += surf[ip+1]*(Qsurf[iq+ii]-surf[ip+5+ii])*(Qsurf[iq+ii]-surf[ip+5+ii]);
					}
					e_sfi += surf[ip+1]*(Qsurf[iq]-surf[ip+5])*(Qsurf[iq+3]-surf[ip+5+3]);
				}
			}
		}
	}
	MPI_Reduce(&e_ldi, &e_ld, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_eli, &e_el, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_sfi, &e_sf, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if (myid==root) e_tot = e_ld + e_el + e_sf;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Win_fence(0, winHsurf);	
}


void cal_stress()
{
	int id, iq, ib, ii, jj, mm, nid[6], ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, h[3][3], M[3][3], dQ[3][3], dQ2, temp, sum, eld;
	real *p=NULL;
	
	for (id=0, iq=0, ib=0; id<lpoint; id++, iq+=5, ib+=6){
		if (info[id]==-1) {
			p  = &Q[iq];
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[ib+ii];
			}
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2xQ[ii]= Q[ip2+ii] - 2.0 * Q[iq+ii] + Q[ip1+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2xQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Qsurf[ip1+1+ii];
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dxQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2xQ[ii]= eight3rd*Qsurf[ip2+1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dxQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2xQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2yQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2yQ[ii]= four3rd*Q[ip2+ii]-4.0*(*(p+ii))+eight3rd*Qsurf[ip1+1+ii];
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dyQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
					d2yQ[ii]= eight3rd*Qsurf[ip2+1+ii]-4.0*(*(p+ii))+four3rd*Q[ip1+ii];
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dyQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2yQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
					d2zQ[ii]= Q[ip1+ii] - 2.0 * Q[iq+ii] + Q[ip2+ii];
				}
			} else if (ip1%5!=0 && ip2%5==0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
					d2zQ[ii]= four3rd*(Q[ip2+ii]-*(p+ii))+eight3rd*(Qsurf[ip1+1+ii]-*(p+ii));
				}
			} else if (ip1%5==0 && ip2%5!=0) {
				for (ii=0; ii<5; ii++) {
					dzQ[ii] =-(*(p+ii)+third*Q[ip1+ii]-four3rd*Qsurf[ip2+1+ii]);
					d2zQ[ii]= eight3rd*(Qsurf[ip2+1+ii]-*(p+ii))+four3rd*(Q[ip1+ii]-*(p+ii));
				}
                        } else if (ip1%5!=0 && ip2%5!=0) {
                                for (ii=0; ii<5; ii++) {
                                        dzQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                                        d2zQ[ii]= 4.0*(Qsurf[ip2+1+ii] - 2.0 * Q[iq+ii] + Qsurf[ip1+1+ii]);
                                }
			} else {
				printf("surface too close\n");
				exit(-1);
			}
			
			trqq = trQQ(p);
			
			for (ii=0; ii<5; ii++) {
				H[iq+ii] = -1*((1.0-third*U_lc)*(*(p+ii)) - U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - kappa * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
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

			sigma_p[id*3]  =( -2.0 * (H[iq]*dxQ[0] + H[iq+1]*dxQ[1] + H[iq+2]*dxQ[2] + H[iq+3]*dxQ[3] + H[iq+4]*dxQ[4]) - H[iq]*dxQ[3] - H[iq+3]*dxQ[0] )*strcf;
			sigma_p[id*3+1]=( -2.0 * (H[iq]*dyQ[0] + H[iq+1]*dyQ[1] + H[iq+2]*dyQ[2] + H[iq+3]*dyQ[3] + H[iq+4]*dyQ[4]) - H[iq]*dyQ[3] - H[iq+3]*dyQ[0] )*strcf;
			sigma_p[id*3+2]=( -2.0 * (H[iq]*dzQ[0] + H[iq+1]*dzQ[1] + H[iq+2]*dzQ[2] + H[iq+3]*dzQ[3] + H[iq+4]*dzQ[4]) - H[iq]*dzQ[3] - H[iq+3]*dzQ[0] )*strcf; 
			
			sum = 2.0*(Q[iq]*H[iq]+Q[iq+3]*H[iq+3]+Q[iq+1]*H[iq+1]+Q[iq+2]*H[iq+2]+Q[iq+4]*H[iq+4])+Q[iq]*H[iq+3]+Q[iq+3]*H[iq];

                        trqq = trQQ(p);
                        qqq = QQQ(p);
                        eld = 0.5 * (1.0-third*U_lc)*trqq - third*U_lc*qqq + 0.25*U_lc*trqq*trqq - Fld0;
			
			for(ii=0;ii<3;ii++){
				for(jj=0;jj<3;jj++){
					temp = (0.*kappa*dQ2)*Delta(ii,jj) + 2*xi*M[ii][jj]*sum;

					
					for(mm=0;mm<3;mm++){
						temp += (1.0-xi)*M[ii][mm]*h[mm][jj]-(xi+1.0)*h[ii][mm]*M[mm][jj];
					}
					
//					temp += -kappa * dQ[ii][jj];
					sigma_q[id*9+ii*3+jj] = temp * strcf;
				}
			}
		}
	}
	MPI_Win_fence(0, wins);
	MPI_Barrier(MPI_COMM_WORLD);
}


void cal_sigma_p()
{
	int id0, id, ip, is, idxm, idxp, idym, idyp, idzm, idzp, idzmm, idzpp;
	
	for (id0=0; id0<lpoint; id0++) {
		is = id0 * 3;
		id = is  * 3;
		ip = is  * 5;
		
		if (info[id0]==-1) {
                        idxp = nextf[ip+2]>=0?(int)(nextf[ip+2]/15)*9:(int)(nextf[ip+2]/15-1)*9;
                        idxm = nextf[ip+1]>=0?(int)(nextf[ip+1]/15)*9:(int)(nextf[ip+1]/15-1)*9;
                        idyp = nextf[ip+4]>=0?(int)(nextf[ip+4]/15)*9:(int)(nextf[ip+4]/15-1)*9;
                        idym = nextf[ip+3]>=0?(int)(nextf[ip+3]/15)*9:(int)(nextf[ip+3]/15-1)*9;
                        idzp = nextf[ip+6]>=0?(int)(nextf[ip+6]/15)*9:(int)(nextf[ip+6]/15-1)*9;
                        idzm = nextf[ip+5]>=0?(int)(nextf[ip+5]/15)*9:(int)(nextf[ip+5]/15-1)*9;
			idzpp= nextf[nextf[ip+6]];
                        idzpp= idzpp>=0?(int)(idzpp/15)*9:(int)(idzpp/15-1)*9;
			idzmm= nextf[nextf[ip+5]];
                        idzmm= idzmm>=0?(int)(idzmm/15)*9:(int)(idzmm/15-1)*9;
			
			if (npar==0 || idxm!=id && idxp!=id) {
				sigma_p[is]  += 0.5 * (sigma_q[idxp]  - sigma_q[idxm]);
				sigma_p[is+1]+= 0.5 * (sigma_q[idxp+3]- sigma_q[idxm+3]);
				sigma_p[is+2]+= 0.5 * (sigma_q[idxp+6]- sigma_q[idxm+6]);
			} else if (idxm==id) {
				sigma_p[is]  += sigma_q[idxp]  - sigma_q[id];
				sigma_p[is+1]+= sigma_q[idxp+3]- sigma_q[id+3];
				sigma_p[is+2]+= sigma_q[idxp+6]- sigma_q[id+6];
			} else if (idxp==id) {
				sigma_p[is]  += sigma_q[id]  - sigma_q[idxm];
				sigma_p[is+1]+= sigma_q[id+3]- sigma_q[idxm+3];
				sigma_p[is+2]+= sigma_q[id+6]- sigma_q[idxm+6];
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
			} else if (wall_on!=0 && (id0+myid*point)/bulk0==Nz-1) {
				sigma_p[is]  +=  0.5*sigma_q[idzmm+2] + 1.5*sigma_q[id+2] - 2.0*sigma_q[idzm+2];
				sigma_p[is+1]+=  0.5*sigma_q[idzmm+5] + 1.5*sigma_q[id+5] - 2.0*sigma_q[idzm+5];
				sigma_p[is+2]+=  0.5*sigma_q[idzmm+8] + 1.5*sigma_q[id+8] - 2.0*sigma_q[idzm+8];
			} else if (wall_on!=0 && (id0+myid*point)/bulk0==0) {
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
