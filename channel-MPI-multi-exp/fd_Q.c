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

	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		for (id=0; id<node; id++) {
			if (id+myid*node<nsurf){
				iq = id * 5;
				ip = iq * 2;
				if (surf[ip]==1) {
					for (ii=0; ii<5; ii++) {
						Qsurf[iq+ii] += -qdt*Gamma_rot*(-Hsurf[iq+ii]+surf[ip+1]*(Qsurf[iq+ii]-surf[ip+5+ii]));
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
						Qsurf[iq+ii] += - qdt*Gamma_rot * (-Hsurf[iq+ii] + 2.0*surf[ip+1]*( (Qtildevec[ii]-Qperpvec[ii]) - temp));
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
	int id, iq, ib, iu, ii, nid[6], flag=0, ip, ip1, ip2, ip3, bid, sii, sjj, herid, mp1, mp2, inext;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, eld, eel, ech, t1[6], t2[6];
	double e_ldi, e_eli, e_chi, e_sfi;
	real *p=NULL;

	int ip0, dxy, dyx, dxz, dzx, dyz, dzy;
	double dxyQ[5], dxzQ[5], dyzQ[5], dQm[5], dQp[5], Qikckj[3][3], dQ[3][3], QdQ[5], Qijcj[3], Qc2[5], dQ2;
	double e_L1i, e_L2i, e_L3i, e_L4i;

	double ttemp, tt, nu_dot_Q[3], nu_dot_dQ, nu[3];
	int i, j, k, iim, iip, jjm, jjp, kkm, kkp, id1, id2, id3, id4;

	if (t_current%t_print==0) {
		flag  = 1;
		e_tot = 0;
		e_ldi = 0;
		e_eli = 0;
		e_chi = 0;
		e_sfi = 0;
        e_L1i = 0;
        e_L2i = 0;
        e_L3i = 0;
		e_L4i = 0;
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
				H[iq+ii] = -1*(A_ldg*(1.0-third*U_lc)*(*(p+ii)) - A_ldg*U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - L1 * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
			}

			if (q_ch>0) {
				H[iq]  -= 2.*( dyQ[2] - dzQ[1] ) * twqL_ch;
				H[iq+3]-= 2.*( dzQ[1] - dxQ[4] ) * twqL_ch;
				H[iq+1]-= (  dyQ[4] - dzQ[3] + dzQ[0] - dxQ[2] ) * twqL_ch;
				H[iq+2]-= (-(dyQ[0] + dyQ[3])- dzQ[4] + dxQ[1] - dyQ[0] ) * twqL_ch;
				H[iq+4]-= (  dzQ[2] +(dxQ[0] + dxQ[3])+ dxQ[3] - dyQ[1] ) * twqL_ch;
			}

			dQ[0][0]=2.*(dxQ[0]*dxQ[0]+dxQ[1]*dxQ[1]+dxQ[2]*dxQ[2]+dxQ[3]*dxQ[3]+dxQ[4]*dxQ[4]+dxQ[0]*dxQ[3]);
			dQ[1][1]=2.*(dyQ[0]*dyQ[0]+dyQ[1]*dyQ[1]+dyQ[2]*dyQ[2]+dyQ[3]*dyQ[3]+dyQ[4]*dyQ[4]+dyQ[0]*dyQ[3]);
			dQ[2][2]=2.*(dzQ[0]*dzQ[0]+dzQ[1]*dzQ[1]+dzQ[2]*dzQ[2]+dzQ[3]*dzQ[3]+dzQ[4]*dzQ[4]+dzQ[0]*dzQ[3]);
			dQ[0][1]=2.*(dxQ[0]*dyQ[0]+dxQ[1]*dyQ[1]+dxQ[2]*dyQ[2]+dxQ[3]*dyQ[3]+dxQ[4]*dyQ[4]) + dxQ[0]*dyQ[3] + dxQ[3]*dyQ[0];
			dQ[0][2]=2.*(dxQ[0]*dzQ[0]+dxQ[1]*dzQ[1]+dxQ[2]*dzQ[2]+dxQ[3]*dzQ[3]+dxQ[4]*dzQ[4]) + dxQ[0]*dzQ[3] + dxQ[3]*dzQ[0];
			dQ[1][2]=2.*(dyQ[0]*dzQ[0]+dyQ[1]*dzQ[1]+dyQ[2]*dzQ[2]+dyQ[3]*dzQ[3]+dyQ[4]*dzQ[4]) + dyQ[0]*dzQ[3] + dyQ[3]*dzQ[0];
			if (L2+L4!=0 || L3!=0) {

				if (nid[0]%5==0 && nid[1]%5==0) {
					ip0 = nid[0];
					ip1 = next_neighbor(id,0,2);
					ip2 = next_neighbor(id,0,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[1];
					ip1 = next_neighbor(id,1,2);
					ip2 = next_neighbor(id,1,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxyQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[2]%5==0 && nid[3]%5==0) {
					ip0 = nid[2];
					ip1 = next_neighbor(id,2,0);
					ip2 = next_neighbor(id,2,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[3];
					ip1 = next_neighbor(id,3,0);
					ip2 = next_neighbor(id,3,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxyQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dxyQ[ii] = 0;

					dxy = -1;
					if (nid[0]%5==0) {
						dxy = 0;
						ip0 = nid[0];
					} else if (nid[1]%5==0) {
						dxy = 1;
						ip0 = nid[1];
					} 
					if (dxy>=0) {
						ip1 = next_neighbor(id,dxy,2);
						ip2 = next_neighbor(id,dxy,3);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxyQ[ii] += (2*dxy-1) * ( dQm[ii] - dyQ[ii] );
					}

					dyx = -1;
					if (nid[2]%5==0) {
						dyx = 2;
						ip0 = nid[2];
					} else if (nid[3]%5==0) {
						dyx = 3;
						ip0 = nid[3];
					} 
					if (dyx>=0) {
						ip1 = next_neighbor(id,dyx,0);
						ip2 = next_neighbor(id,dyx,1);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxyQ[ii] += (2*dyx-5) * ( dQp[ii] - dxQ[ii] );
					}

					if (dxy>=0 && dyx>=0) {
						 for (ii=0; ii<5; ii++) dxyQ[ii] *= 0.5;
					} else if (dxy<0 && dyx<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}	

				if (nid[2]%5==0 && nid[3]%5==0) {
					ip0 = nid[2];
					ip1 = next_neighbor(id,2,4);
					ip2 = next_neighbor(id,2,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[3];
					ip1 = next_neighbor(id,3,4);
					ip2 = next_neighbor(id,3,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dyzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[4]%5==0 && nid[5]%5==0) {
					ip0 = nid[4];
					ip1 = next_neighbor(id,4,2);
					ip2 = next_neighbor(id,4,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[5];
					ip1 = next_neighbor(id,5,2);
					ip2 = next_neighbor(id,5,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dyzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dyzQ[ii] = 0;

					dyz = -1;
					if (nid[2]%5==0) {
						dyz = 2;
						ip0 = nid[2];
					} else if (nid[3]%5==0) {
						dyz = 3;
						ip0 = nid[3];
					} 
					if (dyz>=0) {
						ip1 = next_neighbor(id,dyz,4);
						ip2 = next_neighbor(id,dyz,5);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dyzQ[ii] += (2*dyz-5) * ( dQm[ii] - dzQ[ii] );
					}

					dzy = -1;
					if (nid[4]%5==0) {
						dzy = 4;
						ip0 = nid[4];
					} else if (nid[5]%5==0) {
						dzy = 5;
						ip0 = nid[5];
					} 
					if (dzy>=0) {
						ip1 = next_neighbor(id,dzy,2);
						ip2 = next_neighbor(id,dzy,3);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dyzQ[ii] += (2*dzy-9) * ( dQp[ii] - dyQ[ii] );
					}

					if (dyz>=0 && dzy>=0) {
						 for (ii=0; ii<5; ii++) dyzQ[ii] *= 0.5;
					} else if (dyz<0 && dzy<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}	

				if (nid[0]%5==0 && nid[1]%5==0) {
					ip0 = nid[0];
					ip1 = next_neighbor(id,0,4);
					ip2 = next_neighbor(id,0,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[1];
					ip1 = next_neighbor(id,1,4);
					ip2 = next_neighbor(id,1,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[4]%5==0 && nid[5]%5==0) {
					ip0 = nid[4];
					ip1 = next_neighbor(id,4,0);
					ip2 = next_neighbor(id,4,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[5];
					ip1 = next_neighbor(id,5,0);
					ip2 = next_neighbor(id,5,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dxzQ[ii] = 0;

					dxz = -1;
					if (nid[0]%5==0) {
						dxz = 0;
						ip0 = nid[0];
					} else if (nid[1]%5==0) {
						dxz = 1;
						ip0 = nid[1];
					} 
					if (dxz>=0) {
						ip1 = next_neighbor(id,dxz,4);
						ip2 = next_neighbor(id,dxz,5);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxzQ[ii] += (2*dxz-1) * ( dQm[ii] - dzQ[ii] );
					}

					dzx = -1;
					if (nid[4]%5==0) {
						dzx = 4;
						ip0 = nid[4];
					} else if (nid[5]%5==0) {
						dzx = 5;
						ip0 = nid[5];
					} 
					if (dzx>=0) {
						ip1 = next_neighbor(id,dzx,0);
						ip2 = next_neighbor(id,dzx,1);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxzQ[ii] += (2*dzx-9) * ( dQp[ii] - dxQ[ii] );
					}

					if (dxz>0 && dzx>0) {
						 for (ii=0; ii<5; ii++) dxzQ[ii] *= 0.5;
					} else if (dxz<0 && dzx<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}

				if (1==0) {
					i = (id+myid*point)%Nx;
					j =((id+myid*point)/Nx)%Ny;
					k =((id+myid*point)/Nx)/Ny;
					iim = i-1;
					iip = i+1;
					if (iim<0)   iim+=Nx;
					if (iip>=Nx) iip-=Nx;
					jjm = j-1;
					jjp = j+1;
					if (jjm<0)   jjm+=Ny;
					if (jjp>=Ny) jjp-=Ny;
					kkm = k-1;
					kkp = k+1;
					if (kkm<0)   kkm+=Nz;
					if (kkp>=Nz) kkp-=Nz;
					
					id1 = iim + (jjm + k*Ny)*Nx - myid*point;
					id2 = iim + (jjp + k*Ny)*Nx - myid*point;
					id3 = iip + (jjm + k*Ny)*Nx - myid*point;
					id4 = iip + (jjp + k*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dxyQ[ii]<-1e-15 || tt-dxyQ[ii]>1e-15) printf("wrong xy! diff=%e\n",tt-dxyQ[ii]);
					}
					id1 = iim + (j + kkm*Ny)*Nx - myid*point;
					id2 = iim + (j + kkp*Ny)*Nx - myid*point;
					id3 = iip + (j + kkm*Ny)*Nx - myid*point;
					id4 = iip + (j + kkp*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dxzQ[ii]<-1e-15 || tt-dxzQ[ii]>1e-15) printf("wrong xz! diff=%e\n",tt-dxzQ[ii]);
					}
					id1 = i + (jjm + kkm*Ny)*Nx - myid*point;
					id2 = i + (jjm + kkp*Ny)*Nx - myid*point;
					id3 = i + (jjp + kkm*Ny)*Nx - myid*point;
					id4 = i + (jjp + kkp*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dyzQ[ii]<-1e-15 || tt-dyzQ[ii]>1e-15) printf("wrong yz! diff=%e\n",tt-dyzQ[ii]);
					}
				}

				if (L2+L4!=0) {
					Qikckj[0][0] = d2xQ[0] + dxyQ[1] + dxzQ[2];
					Qikckj[0][1] = dxyQ[0] + d2yQ[1] + dyzQ[2];
					Qikckj[0][2] = dxzQ[0] + dyzQ[1] + d2zQ[2];
					Qikckj[1][0] = d2xQ[1] + dxyQ[3] + dxzQ[4];
					Qikckj[1][1] = dxyQ[1] + d2yQ[3] + dyzQ[4];
					Qikckj[1][2] = dxzQ[1] + dyzQ[3] + d2zQ[4];
					Qikckj[2][0] = d2xQ[2] + dxyQ[4] - dxzQ[0] - dxzQ[3];
					Qikckj[2][1] = dxyQ[2] + d2yQ[4] - dyzQ[0] - dyzQ[3];
					Qikckj[2][2] = dxzQ[2] + dyzQ[4] - d2zQ[0] - d2zQ[3];
					dQ2     = Qikckj[0][0] + Qikckj[1][1] + Qikckj[2][2];

					H[iq]  += (L2+L4)*(      Qikckj[0][0] - third * dQ2   );
					H[iq+1]+= (L2+L4)*( 0.5*(Qikckj[0][1] + Qikckj[1][0]) );
                                        H[iq+2]+= (L2+L4)*( 0.5*(Qikckj[0][2] + Qikckj[2][0]) );
					H[iq+3]+= (L2+L4)*(      Qikckj[1][1] - third * dQ2   );
					H[iq+4]+= (L2+L4)*( 0.5*(Qikckj[1][2] + Qikckj[2][1]) );
				}

				if (L3!=0) {
					dQ2     = dQ[0][0] + dQ[1][1] + dQ[2][2];
					QdQ[0]  = Q[iq]*(d2xQ[0]-d2zQ[0]) + Q[iq+3]*(d2yQ[0]-d2zQ[0]) + 2.0 * ( Q[iq+1]*dxyQ[0] + Q[iq+2]*dxzQ[0] + Q[iq+4]*dyzQ[0] );
					QdQ[1]  = Q[iq]*(d2xQ[1]-d2zQ[1]) + Q[iq+3]*(d2yQ[1]-d2zQ[1]) + 2.0 * ( Q[iq+1]*dxyQ[1] + Q[iq+2]*dxzQ[1] + Q[iq+4]*dyzQ[1] );
					QdQ[2]  = Q[iq]*(d2xQ[2]-d2zQ[2]) + Q[iq+3]*(d2yQ[2]-d2zQ[2]) + 2.0 * ( Q[iq+1]*dxyQ[2] + Q[iq+2]*dxzQ[2] + Q[iq+4]*dyzQ[2] );
					QdQ[3]  = Q[iq]*(d2xQ[3]-d2zQ[3]) + Q[iq+3]*(d2yQ[3]-d2zQ[3]) + 2.0 * ( Q[iq+1]*dxyQ[3] + Q[iq+2]*dxzQ[3] + Q[iq+4]*dyzQ[3] );
					QdQ[4]  = Q[iq]*(d2xQ[4]-d2zQ[4]) + Q[iq+3]*(d2yQ[4]-d2zQ[4]) + 2.0 * ( Q[iq+1]*dxyQ[4] + Q[iq+2]*dxzQ[4] + Q[iq+4]*dyzQ[4] );

					Qijcj[0]= dxQ[0] + dyQ[1] + dzQ[2];
					Qijcj[1]= dxQ[1] + dyQ[3] + dzQ[4];
					Qijcj[2]= dxQ[2] + dyQ[4] - dzQ[0] - dzQ[3];
					Qc2[0]  = dxQ[0]*Qijcj[0] + dyQ[0]*Qijcj[1] + dzQ[0]*Qijcj[2];
                    Qc2[1]  = dxQ[1]*Qijcj[0] + dyQ[1]*Qijcj[1] + dzQ[1]*Qijcj[2];
                    Qc2[2]  = dxQ[2]*Qijcj[0] + dyQ[2]*Qijcj[1] + dzQ[2]*Qijcj[2];
                    Qc2[3]  = dxQ[3]*Qijcj[0] + dyQ[3]*Qijcj[1] + dzQ[3]*Qijcj[2];
                    Qc2[4]  = dxQ[4]*Qijcj[0] + dyQ[4]*Qijcj[1] + dzQ[4]*Qijcj[2];

					H[iq]  +=-L3*( 0.5*(dQ[0][0] - third*dQ2) - QdQ[0] - Qc2[0] );
					H[iq+1]+=-L3*( 0.5* dQ[0][1]              - QdQ[1] - Qc2[1] );
					H[iq+2]+=-L3*( 0.5* dQ[0][2]              - QdQ[2] - Qc2[2] );
					H[iq+3]+=-L3*( 0.5*(dQ[1][1] - third*dQ2) - QdQ[3] - Qc2[3] );
					H[iq+4]+=-L3*( 0.5* dQ[1][2]              - QdQ[4] - Qc2[4] );
				}
			}
			
			if (flag==1) {
				qqq = QQQ(p);
				eld = A_ldg * ( 0.5 * (1.0-third*U_lc)*trqq - third*U_lc*qqq + 0.25*U_lc*trqq*trqq );
				eel = dxQ[0]*dxQ[0] + dyQ[0]*dyQ[0] + dzQ[0]*dzQ[0];				
				for (ii=1; ii<5; ii++) {
					eel += dxQ[ii]*dxQ[ii] + dyQ[ii]*dyQ[ii] + dzQ[ii]*dzQ[ii];
				}
				eel  += dxQ[0]*dxQ[3] + dyQ[0]*dyQ[3] + dzQ[0]*dzQ[3];
				eel  *= L1;

                                e_L1i+= dQ[0][0] + dQ[1][1] + dQ[2][2];
				e_L2i+= (dxQ[0]+dyQ[1]+dzQ[2])*(dxQ[0]+dyQ[1]+dzQ[2])+(dxQ[1]+dyQ[3]+dzQ[4])*(dxQ[1]+dyQ[3]+dzQ[4])+(dxQ[2]+dyQ[4]-dzQ[0]-dzQ[3])*(dxQ[2]+dyQ[4]-dzQ[0]-dzQ[3]);
				e_L3i+= Q[iq]*(dQ[0][0]-dQ[2][2]) + Q[iq+3]*(dQ[1][1]-dQ[2][2]) + 2.* ( Q[iq+1]*dQ[0][1] + Q[iq+2]*dQ[0][2] + Q[iq+4]*dQ[1][2] );
				e_L4i+= dxQ[0]*dxQ[0]+dxQ[1]*dxQ[1]+dxQ[2]*dxQ[2]+dyQ[1]*dyQ[1]+dyQ[3]*dyQ[3]+dyQ[4]*dyQ[4]+dzQ[2]*dzQ[2]+dzQ[4]*dzQ[4]+(dzQ[0]+dzQ[3])*(dzQ[0]+dzQ[3]);
				e_L4i+= 2.*( dyQ[0]*dxQ[1] + dyQ[1]*dxQ[3] + dyQ[2]*dxQ[4] );
				e_L4i+= 2.*( dzQ[0]*dxQ[2] + dzQ[1]*dxQ[4] - dzQ[2]*(dxQ[0]+dxQ[3]) ); 
				e_L4i+= 2.*( dzQ[1]*dyQ[2] + dzQ[3]*dyQ[4] - dzQ[4]*(dyQ[0]+dyQ[3]) );
				e_ldi+= eld;
				e_eli+= eel;
				if (q_ch>0) {
					ech   = Q[iq+0] * dyQ[2] + Q[iq+1] * dyQ[4] - Q[iq+2] *(dyQ[0]+dyQ[3]) \
                                               +Q[iq+1] * dzQ[0] + Q[iq+3] * dzQ[1] + Q[iq+4] * dzQ[2] \
                                               +Q[iq+2] * dxQ[1] + Q[iq+4] * dxQ[3]-(Q[iq]+Q[iq+3])*dxQ[4] \
                                               -Q[iq+0] * dzQ[1] - Q[iq+1] * dzQ[3] - Q[iq+2] * dzQ[4] \
                                               -Q[iq+2] * dyQ[0] - Q[iq+4] * dyQ[1]+(Q[iq]+Q[iq+3])*dyQ[2] \
                                               -Q[iq+1] * dxQ[2] - Q[iq+3] * dxQ[4] + Q[iq+4] *(dxQ[0]+dxQ[3]);
					e_chi+= ech*twqL_ch;
				}
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

	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		for (id=0; id<node; id++) {                     // id: surf point id
			if (id+myid*node<nsurf) {
				iq = id * 5;
				ip = iq * 2;
				if (surf[ip]>0) {                       // finite anchoring
					ib = id * 6;
                    for (ii=0; ii<3; ii++) {
                        nu[ii] = surf[ip+2+ii];         // nu[3]: surface normal
                        if (neighbsurf[ib+ii*2]%5==0) {
                            bid = neighbsurf[ib+ii*2]/5;// bid: bulk id nearest to surf point id
                            sii = ii*2;                 // orientation towards surface [0-5]
                            sjj = ii*2+1;               // orientation towards bulk [0-5]
                            if (nu[ii]<0) {
                                sii = ii*2 + 1;
                                sjj = ii*2;
                            }
                        }
                    }
                    herid = (bid+myid*point)/point;     // herid: processor id storing bid
                    for (ii=0; ii<6; ii++) {            // nid[6]: neighbors of bid (x5)
                        nid[ii] = neighb[bid*6+ii];
                        if (nid[ii]%5==0) nid[ii] += (herid-myid)*point*5;
                        else nid[ii] += (herid-myid)*node*5;
                    }
                    if (debug_on!=0 && nid[sii]!=iq-1) {
                        printf("error: bug in surface evolution!: %d %d\n",nid[sii],iq-1);
                        exit(-1);
                    }
                    
                    // x derivative for surf point id
                    ip1 = neighbsurf[ib];
                    ip2 = neighbsurf[ib+1];
                    mp1 = mod(ip1,5);                   // remainder after division of ip1 by 5 [-2~+2]
                    mp2 = mod(ip2,5);                   // remainder after division of ip2 by 5 [-2~+2]
                    if (mp1==0 && mp2==0) {             // arm (+1/2, +3/2), main direction
						for (ii=0; ii<5; ii++) dxQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
						if (nu[0]<0) for (ii=0; ii<5; ii++) dxQ[ii]=-dxQ[ii];
                        
                    } else if (mp1==0 && mp2==-1) {     // arm (+1/2, +1), main direction
                        for (ii=0; ii<5; ii++) dxQ[ii] = -Qsurf[ip2-mp2+ii] + 4*Q[ip1+ii] - 3*Qsurf[iq+ii];
						if (nu[0]<0) for (ii=0; ii<5; ii++) dxQ[ii]=-dxQ[ii];
                        
                    } else if (mp1==1 && mp2==-1) {     // arm (-1/2, +1/2)
						for (ii=0; ii<5; ii++) dxQ[ii] = Qsurf[ip2-mp2+ii] - Qsurf[ip1-mp1+ii];
                        
                    } else if (mp1<0 && mp2>0) {        // arm (-1, +1)
                        if (mp1==-1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii] + Q[nid[0]+ii]);
                        if (mp2== 1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else for (ii=0; ii<5; ii++) t2[ii] = 0.5*(Q[ip2-mp2+ii] + Q[nid[1]+ii]);
                        for (ii=0; ii<5; ii++) dxQ[ii] = 0.5*(t2[ii]-t1[ii]);
                        
                    } else if (mp1*mp2>0) {             // arm (+1, +2)
                        if (mp1>0) inext = 1;           // +x direction
                        else inext = 0;                 // -x direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        if (abs(mp2)==1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else {
                            ip3 = next_neighbor2(bid,inext,inext);
                            for (ii=0; ii<5; ii++) t2[ii]  = 0.5*(Q[ip2-mp2+ii]+Q[ip3+ii]);
                        }
                        for (ii=0; ii<5; ii++) dxQ[ii] = -0.5*t2[ii] + 2.*t1[ii] - 1.5*Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dxQ[ii]=-dxQ[ii];
                        
                    } else if (mp1!=0 && ip2==0) {      // arm (+1)
                        if (mp1>0) inext = 1;           // +x direction
                        else inext = 0;                 // -x direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii]  = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        for (ii=0; ii<5; ii++) dxQ[ii] = t1[ii] - Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dxQ[ii]=-dxQ[ii];
                        
                    } else {
                        printf("error in surface Q's x derivative calculation!");
                        exit(-1);
                    }

                    // y derivative for surf point id
                    ip1 = neighbsurf[ib+2];
                    ip2 = neighbsurf[ib+3];
                    mp1 = mod(ip1,5);                   // remainder after division of ip1 by 5 [-2~+2]
                    mp2 = mod(ip2,5);                   // remainder after division of ip2 by 5 [-2~+2]
                    if (mp1==0 && mp2==0) {             // arm (+1/2, +3/2), main direction
						for (ii=0; ii<5; ii++) dyQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
						if (nu[1]<0) for (ii=0; ii<5; ii++) dyQ[ii]=-dyQ[ii];
                        
                    } else if (mp1==0 && mp2==-1) {     // arm (+1/2, +1), main direction
                        for (ii=0; ii<5; ii++) dyQ[ii] = -Qsurf[ip2-mp2+ii] + 4.*Q[ip1+ii] - 3.*Qsurf[iq+ii];
						if (nu[1]<0) for (ii=0; ii<5; ii++) dyQ[ii]=-dyQ[ii];
                        
                    } else if (mp1==1 && mp2==-1) {     // arm (-1/2, +1/2)
						for (ii=0; ii<5; ii++) dyQ[ii] = Qsurf[ip2-mp2+ii] - Qsurf[ip1-mp1+ii];
                        
                    } else if (mp1<0 && mp2>0) {        // arm (-1, +1)
                        if (mp1==-1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii] + Q[nid[2]+ii]);
                        if (mp2== 1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else for (ii=0; ii<5; ii++) t2[ii] = 0.5*(Q[ip2-mp2+ii] + Q[nid[3]+ii]);
                        for (ii=0; ii<5; ii++) dyQ[ii] = 0.5*(t2[ii]-t1[ii]);
                        
                    } else if (mp1*mp2>0) {             // arm (+1, +2)
                        if (mp1>0) inext = 3;           // +y direction
                        else inext = 2;                 // -y direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        if (abs(mp2)==1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else {
                            ip3 = next_neighbor2(bid,inext,inext);
                            for (ii=0; ii<5; ii++) t2[ii]  = 0.5*(Q[ip2-mp2+ii]+Q[ip3+ii]);
                        }
                        for (ii=0; ii<5; ii++) dyQ[ii] = -0.5*t2[ii] + 2.*t1[ii] - 1.5*Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dyQ[ii]=-dyQ[ii];
                        
                    } else if (mp1!=0 && ip2==0) {      // arm (+1)
                        if (mp1>0) inext = 3;           // +y direction
                        else inext = 2;                 // -y direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii]  = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        for (ii=0; ii<5; ii++) dyQ[ii] = t1[ii] - Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dyQ[ii]=-dyQ[ii];
                        
                    } else {
                        printf("error in surface Q's y derivative calculation!");
                        exit(-1);
                    }

                    // z derivative for surf point id
                    ip1 = neighbsurf[ib+4];
                    ip2 = neighbsurf[ib+5];
                    mp1 = mod(ip1,5);                   // remainder after division of ip1 by 5 [-2~+2]
                    mp2 = mod(ip2,5);                   // remainder after division of ip2 by 5 [-2~+2]
                    if (mp1==0 && mp2==0) {             // arm (+1/2, +3/2), main direction
						for (ii=0; ii<5; ii++) dzQ[ii] = -third*Q[ip2+ii] + 3.0*Q[ip1+ii] -eight3rd*Qsurf[iq+ii];
						if (nu[2]<0) for (ii=0; ii<5; ii++) dzQ[ii]=-dzQ[ii];
                        
                    } else if (mp1==0 && mp2==-1) {     // arm (+1/2, +1), main direction
                        for (ii=0; ii<5; ii++) dzQ[ii] = -Qsurf[ip2-mp2+ii] + 4.*Q[ip1+ii] - 3.*Qsurf[iq+ii];
						if (nu[2]<0) for (ii=0; ii<5; ii++) dzQ[ii]=-dzQ[ii];
                        
                    } else if (mp1==1 && mp2==-1) {     // arm (-1/2, +1/2)
						for (ii=0; ii<5; ii++) dzQ[ii] = Qsurf[ip2-mp2+ii] - Qsurf[ip1-mp1+ii];
                        
                    } else if (mp1<0 && mp2>0) {        // arm (-1, +1)
                        if (mp1==-1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii] + Q[nid[4]+ii]);
                        if (mp2== 1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else for (ii=0; ii<5; ii++) t2[ii] = 0.5*(Q[ip2-mp2+ii] + Q[nid[5]+ii]);
                        for (ii=0; ii<5; ii++) dzQ[ii] = 0.5*(t2[ii]-t1[ii]);
                        
                    } else if (mp1*mp2>0) {             // arm (+1, +2)
                        if (mp1>0) inext = 5;           // +z direction
                        else inext = 4;                 // -z direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii] = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        if (abs(mp2)==1) for (ii=0; ii<5; ii++) t2[ii] = Qsurf[ip2-mp2+ii];
                        else {
                            ip3 = next_neighbor2(bid,inext,inext);
                            for (ii=0; ii<5; ii++) t2[ii]  = 0.5*(Q[ip2-mp2+ii]+Q[ip3+ii]);
                        }
                        for (ii=0; ii<5; ii++) dzQ[ii] = -0.5*t2[ii] + 2.*t1[ii] - 1.5*Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dzQ[ii]=-dzQ[ii];
                        
                    } else if (mp1!=0 && ip2==0) {      // arm (+1)
                        if (mp1>0) inext = 5;           // +z direction
                        else inext = 4;                 // -z direction
                        if (abs(mp1)==1) for (ii=0; ii<5; ii++) t1[ii] = Qsurf[ip1-mp1+ii];
                        else for (ii=0; ii<5; ii++) t1[ii]  = 0.5*(Q[ip1-mp1+ii]+Q[nid[inext]+ii]);
                        for (ii=0; ii<5; ii++) dzQ[ii] = t1[ii] - Qsurf[iq+ii]; 
                        if (mp1<0) for (ii=0; ii<5; ii++) dzQ[ii]=-dzQ[ii];
                        
                    } else {
                        printf("error in surface Q's z derivative calculation!");
                        exit(-1);
                    }

                    for (ii=0; ii<5; ii++) {
                    	Hsurf[iq+ii] = L1 * ( dxQ[ii]*surf[ip+2]+dyQ[ii]*surf[ip+3]+dzQ[ii]*surf[ip+4] );
                    }
                    if (L2!=0 || L4!=0) {
                    	Qijcj[0]    =  dxQ[0] + dyQ[1] + dzQ[2];
                    	Qijcj[1]    =  dxQ[1] + dyQ[3] + dzQ[4];
                    	Qijcj[2]    =  dxQ[2] + dyQ[4] - dzQ[0] - dzQ[3];
                    }
                    if (L2!=0) {
                    	Hsurf[iq]  += L2 * Qijcj[0] * surf[ip+2];
                    	Hsurf[iq+1]+= L2 *(Qijcj[0] * surf[ip+3] + Qijcj[1] * surf[ip+2]) * 0.5;
                    	Hsurf[iq+2]+= L2 *(Qijcj[0] * surf[ip+4] + Qijcj[2] * surf[ip+2]) * 0.5;
                    	Hsurf[iq+3]+= L2 * Qijcj[1] * surf[ip+3];
                    	Hsurf[iq+4]+= L2 *(Qijcj[1] * surf[ip+4] + Qijcj[2] * surf[ip+3]) * 0.5;
                    }
                    if (L3!=0) {
                    	nu_dot_Q[0] = surf[ip+2]*Qsurf[iq]  +surf[ip+3]*Qsurf[iq+1]+surf[ip+4]*Qsurf[iq+2];
                    	nu_dot_Q[1] = surf[ip+2]*Qsurf[iq+1]+surf[ip+3]*Qsurf[iq+3]+surf[ip+4]*Qsurf[iq+4];
                    	nu_dot_Q[2] = surf[ip+2]*Qsurf[iq+2]+surf[ip+3]*Qsurf[iq+4]-surf[ip+4]*(Qsurf[iq]+Qsurf[iq+3]);
                    	Hsurf[iq]  += L3 * (nu_dot_Q[0]*dxQ[0]+nu_dot_Q[1]*dyQ[0]+nu_dot_Q[2]*dzQ[0]);
                    	Hsurf[iq+1]+= L3 * (nu_dot_Q[0]*dxQ[1]+nu_dot_Q[1]*dyQ[1]+nu_dot_Q[2]*dzQ[1]);
                    	Hsurf[iq+2]+= L3 * (nu_dot_Q[0]*dxQ[2]+nu_dot_Q[1]*dyQ[2]+nu_dot_Q[2]*dzQ[2]);
                    	Hsurf[iq+3]+= L3 * (nu_dot_Q[0]*dxQ[3]+nu_dot_Q[1]*dyQ[3]+nu_dot_Q[2]*dzQ[3]);
                    	Hsurf[iq+4]+= L3 * (nu_dot_Q[0]*dxQ[4]+nu_dot_Q[1]*dyQ[4]+nu_dot_Q[2]*dzQ[4]);
                    }
                    if (L4!=0) {
                    	Hsurf[iq]  +=     L4* (surf[ip+2]* dxQ[0]        +surf[ip+3]* dxQ[1]        +surf[ip+4]* dxQ[2]);
                    	Hsurf[iq+1]+= 0.5*L4* (surf[ip+2]*(dxQ[1]+dyQ[0])+surf[ip+3]*(dxQ[3]+dyQ[1])+surf[ip+4]*(dxQ[4]+dyQ[2]));
                    	Hsurf[iq+2]+= 0.5*L4* (surf[ip+2]*(dxQ[2]+dzQ[0])+surf[ip+3]*(dxQ[4]+dzQ[1])+surf[ip+4]*(-dxQ[0]-dxQ[3]+dzQ[2]));
                    	Hsurf[iq+3]+=     L4* (surf[ip+2]* dyQ[1]        +surf[ip+3]* dyQ[3]        +surf[ip+4]* dyQ[4]);
                    	Hsurf[iq+4]+= 0.5*L4* (surf[ip+2]*(dyQ[2]+dzQ[1])+surf[ip+3]*(dyQ[4]+dzQ[3])+surf[ip+4]*(-dyQ[0]-dyQ[3]+dzQ[4]));
                    }
                    if (L2+L4!=0) {
                    	nu_dot_dQ   = surf[ip+2]*Qijcj[0]+surf[ip+3]*Qijcj[1]+surf[ip+4]*Qijcj[2];
                    	Hsurf[iq]  +=-(L2+L4)*third*nu_dot_dQ;
                    	Hsurf[iq+3]+=-(L2+L4)*third*nu_dot_dQ;
                    }
                    if (q_ch>0) {
                    	Hsurf[iq]  +=q_ch*L1*( 2.*(surf[ip+4]*Qsurf[iq+1]-surf[ip+3]*Qsurf[iq+2]) );
                    	Hsurf[iq+3]+=q_ch*L1*( 2.*(surf[ip+2]*Qsurf[iq+4]-surf[ip+4]*Qsurf[iq+1]) );
                    	Hsurf[iq+1]+=q_ch*L1*( surf[ip+4]* Qsurf[iq+3]           -surf[ip+3]* Qsurf[iq+4]           +surf[ip+2]*Qsurf[iq+2]-surf[ip+4]*Qsurf[iq]);
                    	Hsurf[iq+2]+=q_ch*L1*( surf[ip+4]* Qsurf[iq+4]           +surf[ip+3]*(Qsurf[iq]+Qsurf[iq+3])+surf[ip+3]*Qsurf[iq]  -surf[ip+2]*Qsurf[iq+1]);
                    	Hsurf[iq+4]+=q_ch*L1*(-surf[ip+2]*(Qsurf[iq]+Qsurf[iq+3])-surf[ip+4]* Qsurf[iq+2]           +surf[ip+3]*Qsurf[iq+1]-surf[ip+2]*Qsurf[iq+3]);
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
	}
	MPI_Reduce(&e_ldi, &e_ld, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_eli, &e_el, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_sfi, &e_sf, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	MPI_Reduce(&e_chi, &e_ch, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_L1i, &e_L1, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_L2i, &e_L2, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_L3i, &e_L3, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&e_L4i, &e_L4, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	e_L1 *= (0.5 * L1);
        e_L2 *= (0.5 * L2);
        e_L3 *= (0.5 * L3);
        e_L4 *= (0.5 * L4);
	e_el  = e_L1 + e_L2 + e_L3 + e_L4;
	if (myid==root) e_tot = e_ld + e_L1 + e_L2 + e_L3 + e_L4 + e_sf + e_ch;
	MPI_Barrier(MPI_COMM_WORLD);
	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) MPI_Win_fence(0, winHsurf);	
}


void cal_stress2()
{
	int id, iq, ib, iu, ii, nid[6], ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq, trH3rd, eld, eel, ech;
	double e_ldi, e_eli, e_chi, e_sfi;
	real *p=NULL;

	int ip0, dxy, dyx, dxz, dzx, dyz, dzy;
	double dxyQ[5], dxzQ[5], dyzQ[5], dQm[5], dQp[5], Qikckj[3][3], dQ[3][3], QdQ[5], Qijcj[3], Qc2[5], dQ2;
	double e_L1i, e_L2i, e_L3i, e_L4i;

	double ttemp, tt, nu_dot_Q[3], nu_dot_dQ;
	int i, j, k, iim, iip, jjm, jjp, kkm, kkp, id1, id2, id3, id4;
	
    int jj, mm;
    double temp, sum, h[3][3], M[3][3];

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
				H[iq+ii] = -1*(A_ldg*(1.0-third*U_lc)*(*(p+ii)) - A_ldg*U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - L1 * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
			}

			if (q_ch>0) {
				H[iq]  -= 2.*( dyQ[2] - dzQ[1] ) * twqL_ch;
				H[iq+3]-= 2.*( dzQ[1] - dxQ[4] ) * twqL_ch;
				H[iq+1]-= (  dyQ[4] - dzQ[3] + dzQ[0] - dxQ[2] ) * twqL_ch;
				H[iq+2]-= (-(dyQ[0] + dyQ[3])- dzQ[4] + dxQ[1] - dyQ[0] ) * twqL_ch;
				H[iq+4]-= (  dzQ[2] +(dxQ[0] + dxQ[3])+ dxQ[3] - dyQ[1] ) * twqL_ch;
			}

			dQ[0][0]=2.*(dxQ[0]*dxQ[0]+dxQ[1]*dxQ[1]+dxQ[2]*dxQ[2]+dxQ[3]*dxQ[3]+dxQ[4]*dxQ[4]+dxQ[0]*dxQ[3]);
			dQ[1][1]=2.*(dyQ[0]*dyQ[0]+dyQ[1]*dyQ[1]+dyQ[2]*dyQ[2]+dyQ[3]*dyQ[3]+dyQ[4]*dyQ[4]+dyQ[0]*dyQ[3]);
			dQ[2][2]=2.*(dzQ[0]*dzQ[0]+dzQ[1]*dzQ[1]+dzQ[2]*dzQ[2]+dzQ[3]*dzQ[3]+dzQ[4]*dzQ[4]+dzQ[0]*dzQ[3]);
			dQ[0][1]=2.*(dxQ[0]*dyQ[0]+dxQ[1]*dyQ[1]+dxQ[2]*dyQ[2]+dxQ[3]*dyQ[3]+dxQ[4]*dyQ[4]) + dxQ[0]*dyQ[3] + dxQ[3]*dyQ[0];
			dQ[0][2]=2.*(dxQ[0]*dzQ[0]+dxQ[1]*dzQ[1]+dxQ[2]*dzQ[2]+dxQ[3]*dzQ[3]+dxQ[4]*dzQ[4]) + dxQ[0]*dzQ[3] + dxQ[3]*dzQ[0];
			dQ[1][2]=2.*(dyQ[0]*dzQ[0]+dyQ[1]*dzQ[1]+dyQ[2]*dzQ[2]+dyQ[3]*dzQ[3]+dyQ[4]*dzQ[4]) + dyQ[0]*dzQ[3] + dyQ[3]*dzQ[0];
			if (L2+L4!=0 || L3!=0) {

				if (nid[0]%5==0 && nid[1]%5==0) {
					ip0 = nid[0];
					ip1 = next_neighbor(id,0,2);
					ip2 = next_neighbor(id,0,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[1];
					ip1 = next_neighbor(id,1,2);
					ip2 = next_neighbor(id,1,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxyQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[2]%5==0 && nid[3]%5==0) {
					ip0 = nid[2];
					ip1 = next_neighbor(id,2,0);
					ip2 = next_neighbor(id,2,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[3];
					ip1 = next_neighbor(id,3,0);
					ip2 = next_neighbor(id,3,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxyQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dxyQ[ii] = 0;

					dxy = -1;
					if (nid[0]%5==0) {
						dxy = 0;
						ip0 = nid[0];
					} else if (nid[1]%5==0) {
						dxy = 1;
						ip0 = nid[1];
					} 
					if (dxy>=0) {
						ip1 = next_neighbor(id,dxy,2);
						ip2 = next_neighbor(id,dxy,3);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxyQ[ii] += (2*dxy-1) * ( dQm[ii] - dyQ[ii] );
					}

					dyx = -1;
					if (nid[2]%5==0) {
						dyx = 2;
						ip0 = nid[2];
					} else if (nid[3]%5==0) {
						dyx = 3;
						ip0 = nid[3];
					} 
					if (dyx>=0) {
						ip1 = next_neighbor(id,dyx,0);
						ip2 = next_neighbor(id,dyx,1);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxyQ[ii] += (2*dyx-5) * ( dQp[ii] - dxQ[ii] );
					}

					if (dxy>=0 && dyx>=0) {
						 for (ii=0; ii<5; ii++) dxyQ[ii] *= 0.5;
					} else if (dxy<0 && dyx<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}	

				if (nid[2]%5==0 && nid[3]%5==0) {
					ip0 = nid[2];
					ip1 = next_neighbor(id,2,4);
					ip2 = next_neighbor(id,2,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[3];
					ip1 = next_neighbor(id,3,4);
					ip2 = next_neighbor(id,3,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dyzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[4]%5==0 && nid[5]%5==0) {
					ip0 = nid[4];
					ip1 = next_neighbor(id,4,2);
					ip2 = next_neighbor(id,4,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[5];
					ip1 = next_neighbor(id,5,2);
					ip2 = next_neighbor(id,5,3);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dyzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dyzQ[ii] = 0;

					dyz = -1;
					if (nid[2]%5==0) {
						dyz = 2;
						ip0 = nid[2];
					} else if (nid[3]%5==0) {
						dyz = 3;
						ip0 = nid[3];
					} 
					if (dyz>=0) {
						ip1 = next_neighbor(id,dyz,4);
						ip2 = next_neighbor(id,dyz,5);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dyzQ[ii] += (2*dyz-5) * ( dQm[ii] - dzQ[ii] );
					}

					dzy = -1;
					if (nid[4]%5==0) {
						dzy = 4;
						ip0 = nid[4];
					} else if (nid[5]%5==0) {
						dzy = 5;
						ip0 = nid[5];
					} 
					if (dzy>=0) {
						ip1 = next_neighbor(id,dzy,2);
						ip2 = next_neighbor(id,dzy,3);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dyzQ[ii] += (2*dzy-9) * ( dQp[ii] - dyQ[ii] );
					}

					if (dyz>=0 && dzy>=0) {
						 for (ii=0; ii<5; ii++) dyzQ[ii] *= 0.5;
					} else if (dyz<0 && dzy<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}	

				if (nid[0]%5==0 && nid[1]%5==0) {
					ip0 = nid[0];
					ip1 = next_neighbor(id,0,4);
					ip2 = next_neighbor(id,0,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[1];
					ip1 = next_neighbor(id,1,4);
					ip2 = next_neighbor(id,1,5);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else if (nid[4]%5==0 && nid[5]%5==0) {
					ip0 = nid[4];
					ip1 = next_neighbor(id,4,0);
					ip2 = next_neighbor(id,4,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					} 

					ip0 = nid[5];
					ip1 = next_neighbor(id,5,0);
					ip2 = next_neighbor(id,5,1);
					if (ip1%5==0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5==0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
						}
					} else if (ip1%5==0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
						}
					} else if (ip1%5!=0 && ip2%5!=0) {
						for (ii=0; ii<5; ii++) {
							dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
						}
					}

					for (ii=0; ii<5; ii++) dxzQ[ii] = 0.5 * ( dQp[ii] - dQm[ii] );
				} else {
					for (ii=0; ii<5; ii++) dxzQ[ii] = 0;

					dxz = -1;
					if (nid[0]%5==0) {
						dxz = 0;
						ip0 = nid[0];
					} else if (nid[1]%5==0) {
						dxz = 1;
						ip0 = nid[1];
					} 
					if (dxz>=0) {
						ip1 = next_neighbor(id,dxz,4);
						ip2 = next_neighbor(id,dxz,5);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQm[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxzQ[ii] += (2*dxz-1) * ( dQm[ii] - dzQ[ii] );
					}

					dzx = -1;
					if (nid[4]%5==0) {
						dzx = 4;
						ip0 = nid[4];
					} else if (nid[5]%5==0) {
						dzx = 5;
						ip0 = nid[5];
					} 
					if (dzx>=0) {
						ip1 = next_neighbor(id,dzx,0);
						ip2 = next_neighbor(id,dzx,1);
						if (ip1%5==0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5==0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = third*Q[ip2+ii]+Q[ip0+ii]-four3rd*Qsurf[ip1+1+ii];
							}
						} else if (ip1%5==0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] =-(third*Q[ip1+ii]+Q[ip0+ii]-four3rd*Qsurf[ip2+1+ii]);
							}
						} else if (ip1%5!=0 && ip2%5!=0) {
							for (ii=0; ii<5; ii++) {
								dQp[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
							}
						}
						for (ii=0; ii<5; ii++) dxzQ[ii] += (2*dzx-9) * ( dQp[ii] - dxQ[ii] );
					}

					if (dxz>0 && dzx>0) {
						 for (ii=0; ii<5; ii++) dxzQ[ii] *= 0.5;
					} else if (dxz<0 && dzx<0 && debug_on!=0) {
						printf("point without bulk neighbors detected\n");
					}
				}

				if (1==0) {
					i = (id+myid*point)%Nx;
					j =((id+myid*point)/Nx)%Ny;
					k =((id+myid*point)/Nx)/Ny;
					iim = i-1;
					iip = i+1;
					if (iim<0)   iim+=Nx;
					if (iip>=Nx) iip-=Nx;
					jjm = j-1;
					jjp = j+1;
					if (jjm<0)   jjm+=Ny;
					if (jjp>=Ny) jjp-=Ny;
					kkm = k-1;
					kkp = k+1;
					if (kkm<0)   kkm+=Nz;
					if (kkp>=Nz) kkp-=Nz;
					
					id1 = iim + (jjm + k*Ny)*Nx - myid*point;
					id2 = iim + (jjp + k*Ny)*Nx - myid*point;
					id3 = iip + (jjm + k*Ny)*Nx - myid*point;
					id4 = iip + (jjp + k*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dxyQ[ii]<-1e-15 || tt-dxyQ[ii]>1e-15) printf("wrong xy! diff=%e\n",tt-dxyQ[ii]);
					}
					id1 = iim + (j + kkm*Ny)*Nx - myid*point;
					id2 = iim + (j + kkp*Ny)*Nx - myid*point;
					id3 = iip + (j + kkm*Ny)*Nx - myid*point;
					id4 = iip + (j + kkp*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dxzQ[ii]<-1e-15 || tt-dxzQ[ii]>1e-15) printf("wrong xz! diff=%e\n",tt-dxzQ[ii]);
					}
					id1 = i + (jjm + kkm*Ny)*Nx - myid*point;
					id2 = i + (jjm + kkp*Ny)*Nx - myid*point;
					id3 = i + (jjp + kkm*Ny)*Nx - myid*point;
					id4 = i + (jjp + kkp*Ny)*Nx - myid*point;
					for (ii=0; ii<5; ii++) {
						tt = 0.25 * ( Q[id1*5+ii] - Q[id2*5+ii] - Q[id3*5+ii] + Q[id4*5+ii] );
						if (tt-dyzQ[ii]<-1e-15 || tt-dyzQ[ii]>1e-15) printf("wrong yz! diff=%e\n",tt-dyzQ[ii]);
					}
				}

				if (L2+L4!=0) {
					Qikckj[0][0] = d2xQ[0] + dxyQ[1] + dxzQ[2];
					Qikckj[0][1] = dxyQ[0] + d2yQ[1] + dyzQ[2];
					Qikckj[0][2] = dxzQ[0] + dyzQ[1] + d2zQ[2];
					Qikckj[1][0] = d2xQ[1] + dxyQ[3] + dxzQ[4];
					Qikckj[1][1] = dxyQ[1] + d2yQ[3] + dyzQ[4];
					Qikckj[1][2] = dxzQ[1] + dyzQ[3] + d2zQ[4];
					Qikckj[2][0] = d2xQ[2] + dxyQ[4] - dxzQ[0] - dxzQ[3];
					Qikckj[2][1] = dxyQ[2] + d2yQ[4] - dyzQ[0] - dyzQ[3];
					Qikckj[2][2] = dxzQ[2] + dyzQ[4] - d2zQ[0] - d2zQ[3];
					dQ2     = Qikckj[0][0] + Qikckj[1][1] + Qikckj[2][2];

					H[iq]  += (L2+L4)*(      Qikckj[0][0] - third * dQ2   );
					H[iq+1]+= (L2+L4)*( 0.5*(Qikckj[0][1] + Qikckj[1][0]) );
                                        H[iq+2]+= (L2+L4)*( 0.5*(Qikckj[0][2] + Qikckj[2][0]) );
					H[iq+3]+= (L2+L4)*(      Qikckj[1][1] - third * dQ2   );
					H[iq+4]+= (L2+L4)*( 0.5*(Qikckj[1][2] + Qikckj[2][1]) );
				}

				if (L3!=0) {
					dQ2     = dQ[0][0] + dQ[1][1] + dQ[2][2];
					QdQ[0]  = Q[iq]*(d2xQ[0]-d2zQ[0]) + Q[iq+3]*(d2yQ[0]-d2zQ[0]) + 2.0 * ( Q[iq+1]*dxyQ[0] + Q[iq+2]*dxzQ[0] + Q[iq+4]*dyzQ[0] );
					QdQ[1]  = Q[iq]*(d2xQ[1]-d2zQ[1]) + Q[iq+3]*(d2yQ[1]-d2zQ[1]) + 2.0 * ( Q[iq+1]*dxyQ[1] + Q[iq+2]*dxzQ[1] + Q[iq+4]*dyzQ[1] );
					QdQ[2]  = Q[iq]*(d2xQ[2]-d2zQ[2]) + Q[iq+3]*(d2yQ[2]-d2zQ[2]) + 2.0 * ( Q[iq+1]*dxyQ[2] + Q[iq+2]*dxzQ[2] + Q[iq+4]*dyzQ[2] );
					QdQ[3]  = Q[iq]*(d2xQ[3]-d2zQ[3]) + Q[iq+3]*(d2yQ[3]-d2zQ[3]) + 2.0 * ( Q[iq+1]*dxyQ[3] + Q[iq+2]*dxzQ[3] + Q[iq+4]*dyzQ[3] );
					QdQ[4]  = Q[iq]*(d2xQ[4]-d2zQ[4]) + Q[iq+3]*(d2yQ[4]-d2zQ[4]) + 2.0 * ( Q[iq+1]*dxyQ[4] + Q[iq+2]*dxzQ[4] + Q[iq+4]*dyzQ[4] );

					Qijcj[0]= dxQ[0] + dyQ[1] + dzQ[2];
					Qijcj[1]= dxQ[1] + dyQ[3] + dzQ[4];
					Qijcj[2]= dxQ[2] + dyQ[4] - dzQ[0] - dzQ[3];
					Qc2[0]  = dxQ[0]*Qijcj[0] + dyQ[0]*Qijcj[1] + dzQ[0]*Qijcj[2];
                    Qc2[1]  = dxQ[1]*Qijcj[0] + dyQ[1]*Qijcj[1] + dzQ[1]*Qijcj[2];
                    Qc2[2]  = dxQ[2]*Qijcj[0] + dyQ[2]*Qijcj[1] + dzQ[2]*Qijcj[2];
                    Qc2[3]  = dxQ[3]*Qijcj[0] + dyQ[3]*Qijcj[1] + dzQ[3]*Qijcj[2];
                    Qc2[4]  = dxQ[4]*Qijcj[0] + dyQ[4]*Qijcj[1] + dzQ[4]*Qijcj[2];

					H[iq]  +=-L3*( 0.5*(dQ[0][0] - third*dQ2) - QdQ[0] - Qc2[0] );
					H[iq+1]+=-L3*( 0.5* dQ[0][1]              - QdQ[1] - Qc2[1] );
					H[iq+2]+=-L3*( 0.5* dQ[0][2]              - QdQ[2] - Qc2[2] );
					H[iq+3]+=-L3*( 0.5*(dQ[1][1] - third*dQ2) - QdQ[3] - Qc2[3] );
					H[iq+4]+=-L3*( 0.5* dQ[1][2]              - QdQ[4] - Qc2[4] );
				}
			}
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

			sigma_p[id*3]  =-2.0 * (H[iq]*dxQ[0] + H[iq+1]*dxQ[1] + H[iq+2]*dxQ[2] + H[iq+3]*dxQ[3] + H[iq+4]*dxQ[4]) - H[iq]*dxQ[3] - H[iq+3]*dxQ[0];
			sigma_p[id*3+1]=-2.0 * (H[iq]*dyQ[0] + H[iq+1]*dyQ[1] + H[iq+2]*dyQ[2] + H[iq+3]*dyQ[3] + H[iq+4]*dyQ[4]) - H[iq]*dyQ[3] - H[iq+3]*dyQ[0];
			sigma_p[id*3+2]=-2.0 * (H[iq]*dzQ[0] + H[iq+1]*dzQ[1] + H[iq+2]*dzQ[2] + H[iq+3]*dzQ[3] + H[iq+4]*dzQ[4]) - H[iq]*dzQ[3] - H[iq+3]*dzQ[0]; 
			
			sum = 2.0*(Q[iq]*H[iq]+Q[iq+3]*H[iq+3]+Q[iq+1]*H[iq+1]+Q[iq+2]*H[iq+2]+Q[iq+4]*H[iq+4])+Q[iq]*H[iq+3]+Q[iq+3]*H[iq];
			
			for(ii=0;ii<3;ii++){
				for(jj=0;jj<3;jj++){
					temp = 2*xi*M[ii][jj]*sum;

					for(mm=0;mm<3;mm++){
						temp += (1.0-xi)*M[ii][mm]*h[mm][jj]-(xi+1.0)*h[ii][mm]*M[mm][jj];
					}
					
					sigma_q[id*9+ii*3+jj] = temp;
				}
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
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
				H[iq+ii] = -1*(A_ldg*(1.0-third*U_lc)*(*(p+ii)) - A_ldg*U_lc*( QQ(p,ii) - trqq*(*(p+ii) + OneThirdDelta(ii))) - L1 * (d2xQ[ii]+d2yQ[ii]+d2zQ[ii]));
			}

			if (q_ch>0) {
				H[iq]  -= 2.*( dyQ[2] - dzQ[1] ) * twqL_ch;
				H[iq+3]-= 2.*( dzQ[1] - dxQ[4] ) * twqL_ch;
				H[iq+1]-= (  dyQ[4] - dzQ[3] + dzQ[0] - dxQ[2] ) * twqL_ch;
				H[iq+2]-= (-(dyQ[0] + dyQ[3])- dzQ[4] + dxQ[1] - dyQ[0] ) * twqL_ch;
				H[iq+4]-= (  dzQ[2] +(dxQ[0] + dxQ[3])+ dxQ[3] - dyQ[1] ) * twqL_ch;
			}
			
			dQ[0][0]=2.*(dxQ[0]*dxQ[0]+dxQ[1]*dxQ[1]+dxQ[2]*dxQ[2]+dxQ[3]*dxQ[3]+dxQ[4]*dxQ[4]+dxQ[0]*dxQ[3]);
			dQ[1][1]=2.*(dyQ[0]*dyQ[0]+dyQ[1]*dyQ[1]+dyQ[2]*dyQ[2]+dyQ[3]*dyQ[3]+dyQ[4]*dyQ[4]+dyQ[0]*dyQ[3]);
			dQ[2][2]=2.*(dzQ[0]*dzQ[0]+dzQ[1]*dzQ[1]+dzQ[2]*dzQ[2]+dzQ[3]*dzQ[3]+dzQ[4]*dzQ[4]+dzQ[0]*dzQ[3]);
                        dQ[0][1]=2.*(dxQ[0]*dyQ[0]+dxQ[1]*dyQ[1]+dxQ[2]*dyQ[2]+dxQ[3]*dyQ[3]+dxQ[4]*dyQ[4]) + dxQ[0]*dyQ[3] + dxQ[3]*dyQ[0];
                        dQ[0][2]=2.*(dxQ[0]*dzQ[0]+dxQ[1]*dzQ[1]+dxQ[2]*dzQ[2]+dxQ[3]*dzQ[3]+dxQ[4]*dzQ[4]) + dxQ[0]*dzQ[3] + dxQ[3]*dzQ[0];
                        dQ[1][2]=2.*(dyQ[0]*dzQ[0]+dyQ[1]*dzQ[1]+dyQ[2]*dzQ[2]+dyQ[3]*dzQ[3]+dyQ[4]*dzQ[4]) + dyQ[0]*dzQ[3] + dyQ[3]*dzQ[0];
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

			sigma_p[id*3]  =-2.0 * (H[iq]*dxQ[0] + H[iq+1]*dxQ[1] + H[iq+2]*dxQ[2] + H[iq+3]*dxQ[3] + H[iq+4]*dxQ[4]) - H[iq]*dxQ[3] - H[iq+3]*dxQ[0];
			sigma_p[id*3+1]=-2.0 * (H[iq]*dyQ[0] + H[iq+1]*dyQ[1] + H[iq+2]*dyQ[2] + H[iq+3]*dyQ[3] + H[iq+4]*dyQ[4]) - H[iq]*dyQ[3] - H[iq+3]*dyQ[0];
			sigma_p[id*3+2]=-2.0 * (H[iq]*dzQ[0] + H[iq+1]*dzQ[1] + H[iq+2]*dzQ[2] + H[iq+3]*dzQ[3] + H[iq+4]*dzQ[4]) - H[iq]*dzQ[3] - H[iq+3]*dzQ[0]; 
			
			sum = 2.0*(Q[iq]*H[iq]+Q[iq+3]*H[iq+3]+Q[iq+1]*H[iq+1]+Q[iq+2]*H[iq+2]+Q[iq+4]*H[iq+4])+Q[iq]*H[iq+3]+Q[iq+3]*H[iq];

                        trqq = trQQ(p);
                        qqq = QQQ(p);
                        eld = A_ldg * ( 0.5 * (1.0-third*U_lc)*trqq - third*U_lc*qqq + 0.25*U_lc*trqq*trqq ) - Fld0;
			
			for(ii=0;ii<3;ii++){
				for(jj=0;jj<3;jj++){
					temp = (0.*L1*dQ2)*Delta(ii,jj) + 2*xi*M[ii][jj]*sum;

					
					for(mm=0;mm<3;mm++){
						temp += (1.0-xi)*M[ii][mm]*h[mm][jj]-(xi+1.0)*h[ii][mm]*M[mm][jj];
					}
					
//					temp += -L1 * dQ[ii][jj];
					sigma_q[id*9+ii*3+jj] = temp;
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
			
			if (idxm!=id && idxp!=id) {
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
		
			if (idym!=id && idyp!=id) {
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
			} else if (wall_z!=0 && (id0+myid*point)/bulk0==Nz-1) {
				sigma_p[is]  +=  0.5*sigma_q[idzmm+2] + 1.5*sigma_q[id+2] - 2.0*sigma_q[idzm+2];
				sigma_p[is+1]+=  0.5*sigma_q[idzmm+5] + 1.5*sigma_q[id+5] - 2.0*sigma_q[idzm+5];
				sigma_p[is+2]+=  0.5*sigma_q[idzmm+8] + 1.5*sigma_q[id+8] - 2.0*sigma_q[idzm+8];
			} else if (wall_z!=0 && (id0+myid*point)/bulk0==0) {
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
