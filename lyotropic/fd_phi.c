/*
 *  fd_phi.c
 *  
 *
 *  Created by Sirius on Nov/13/17.
 *  Copyright 2017 home. All rights reserved.
 *
 */
#include "main.h"

void cal_mu()
{
	int id, iq, ib, iu, ii, nid[6], ip, ip1, ip2;
	double dxQ[6], dyQ[6], dzQ[6], d2xQ[6], d2yQ[6], d2zQ[6], trqq, qqq;
	real *p=NULL;

	int ip0, dxy, dyx, dxz, dzx, dyz, dzy;
	double dxyQ[5], dxzQ[5], dyzQ[5], dQm[5], dQp[5], Qikckj[3][3], dQ[3][3], QdQ[5], Qijcj[3], Qc2[5], dQ2;
	double e_L1i, e_L2i, e_L3i, e_L4i;

	double ttemp, tt, nu_dot_Q[3], nu_dot_dQ;
	int i, j, k, iim, iip, jjm, jjp, kkm, kkp, id1, id2, id3, id4;

    double dx_phi, dy_phi, dz_phi, d2x_phi, d2y_phi, d2z_phi, phi, phi2, dxy_phi, dyz_phi, dxz_phi, fpp, fpm, fmp, fmm;
	
	for (iq=0, ib=0, id=0; id<lpoint; iq+=5, ib+=6, id++){
		if (info[id]==-1) {
			p  = &Q[iq];
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[ib+ii];
			}

            dx_phi  = 0.;
            dy_phi  = 0.;
            dz_phi  = 0.;
            d2x_phi = 0.;
            d2y_phi = 0.;
            d2z_phi = 0.;

            if (Q_on!=0) {
                ip1 = nid[0];
                ip2 = nid[1];
                if (ip1%5==0 && ip2%5==0) {
                	for (ii=0; ii<5; ii++) {
                		dxQ[ii] = 0.5 * (Q[ip2+ii] - Q[ip1+ii]);
                	}
                } else if (ip1%5!=0 && ip2%5==0) {
                	for (ii=0; ii<5; ii++) {
                		dxQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
                	}
                } else if (ip1%5==0 && ip2%5!=0) {
                	for (ii=0; ii<5; ii++) {
                		dxQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
                	}
                } else if (ip1%5!=0 && ip2%5!=0) {
                    for (ii=0; ii<5; ii++) {
                        dxQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
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
                	}
                } else if (ip1%5!=0 && ip2%5==0) {
                	for (ii=0; ii<5; ii++) {
                		dyQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
                	}
                } else if (ip1%5==0 && ip2%5!=0) {
                	for (ii=0; ii<5; ii++) {
                		dyQ[ii] = four3rd*Qsurf[ip2+1+ii]-*(p+ii)-third*Q[ip1+ii];
                	}
                } else if (ip1%5!=0 && ip2%5!=0) {
                    for (ii=0; ii<5; ii++) {
                        dyQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
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
                	}
                } else if (ip1%5!=0 && ip2%5==0) {
                	for (ii=0; ii<5; ii++) {
                		dzQ[ii] = third*Q[ip2+ii]+*(p+ii)-four3rd*Qsurf[ip1+1+ii];
                	}
                } else if (ip1%5==0 && ip2%5!=0) {
                	for (ii=0; ii<5; ii++) {
                		dzQ[ii] =-(third*Q[ip1+ii]+*(p+ii)-four3rd*Qsurf[ip2+1+ii]);
                	}
                } else if (ip1%5!=0 && ip2%5!=0) {
                    for (ii=0; ii<5; ii++) {
                        dzQ[ii] = (Qsurf[ip2+1+ii] - Qsurf[ip1+1+ii]);
                    }
                } else {
                	printf("surface too close\n");
                	exit(-1);
                }
            }
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0) {
                dx_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2x_phi = ly_phi[ip2/5] - 2.*ly_phi[id] + ly_phi[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dx_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2x_phi = ly_phi[ip2/5] - ly_phi[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dx_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2x_phi = ly_phi[ip1/5] - ly_phi[id];
			}

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0) {
                dy_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2y_phi = ly_phi[ip2/5] - 2.*ly_phi[id] + ly_phi[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dy_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2y_phi = ly_phi[ip2/5] - ly_phi[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dy_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2y_phi = ly_phi[ip1/5] - ly_phi[id];
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0) {
                dz_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2z_phi = ly_phi[ip2/5] - 2.*ly_phi[id] + ly_phi[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dz_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2z_phi = ly_phi[ip2/5] - ly_phi[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dz_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2z_phi = ly_phi[ip1/5] - ly_phi[id];
			}

            // cross derivatives
            if (ly_an!=0) {
                phi = ly_phi[id];
                if (nid[0]%5==0) {
                    ip1 = next_neighbor(id,0,2);
                    ip2 = next_neighbor(id,0,3);
                    if (ip1%5==0) fmm = ly_phi[ip1/5];
                    else fmm = ly_phi[nid[0]/5];
                    if (ip2%5==0) fmp = ly_phi[ip2/5];
                    else fmp = ly_phi[nid[0]/5];
                } else {
                    if (nid[2]%5==0) fmm = ly_phi[nid[2]/5];
                    else fmm = phi;
                    if (nid[3]%5==0) fmp = ly_phi[nid[3]/5];
                    else fmp = phi;
                }
                if (nid[1]%5==0) {
                    ip1 = next_neighbor(id,1,2);
                    ip2 = next_neighbor(id,1,3);
                    if (ip1%5==0) fpm = ly_phi[ip1/5];
                    else fpm = ly_phi[nid[1]/5];
                    if (ip2%5==0) fpp = ly_phi[ip2/5];
                    else fpp = ly_phi[nid[1]/5];
                } else {
                    if (nid[2]%5==0) fpm = ly_phi[nid[2]/5];
                    else fpm = phi;
                    if (nid[3]%5==0) fpp = ly_phi[nid[3]/5];
                    else fpp = phi;
                }
                dxy_phi = .25*(fmm+fpp-fmp-fpm);

                if (nid[2]%5==0) {
                    ip1 = next_neighbor(id,2,4);
                    ip2 = next_neighbor(id,2,5);
                    if (ip1%5==0) fmm = ly_phi[ip1/5];
                    else fmm = ly_phi[nid[2]/5];
                    if (ip2%5==0) fmp = ly_phi[ip2/5];
                    else fmp = ly_phi[nid[2]/5];
                } else {
                    if (nid[4]%5==0) fmm = ly_phi[nid[4]/5];
                    else fmm = phi;
                    if (nid[5]%5==0) fmp = ly_phi[nid[5]/5];
                    else fmp = phi;
                }
                if (nid[3]%5==0) {
                    ip1 = next_neighbor(id,3,4);
                    ip2 = next_neighbor(id,3,5);
                    if (ip1%5==0) fpm = ly_phi[ip1/5];
                    else fpm = ly_phi[nid[3]/5];
                    if (ip2%5==0) fpp = ly_phi[ip2/5];
                    else fpp = ly_phi[nid[3]/5];
                } else {
                    if (nid[4]%5==0) fpm = ly_phi[nid[4]/5];
                    else fpm = phi;
                    if (nid[5]%5==0) fpp = ly_phi[nid[5]/5];
                    else fpp = phi;
                }
                dyz_phi = .25*(fmm+fpp-fmp-fpm);

                if (nid[0]%5==0) {
                    ip1 = next_neighbor(id,0,4);
                    ip2 = next_neighbor(id,0,5);
                    if (ip1%5==0) fmm = ly_phi[ip1/5];
                    else fmm = ly_phi[nid[0]/5];
                    if (ip2%5==0) fmp = ly_phi[ip2/5];
                    else fmp = ly_phi[nid[0]/5];
                } else {
                    if (nid[4]%5==0) fmm = ly_phi[nid[4]/5];
                    else fmm = phi;
                    if (nid[5]%5==0) fmp = ly_phi[nid[5]/5];
                    else fmp = phi;
                }
                if (nid[1]%5==0) {
                    ip1 = next_neighbor(id,1,4);
                    ip2 = next_neighbor(id,1,5);
                    if (ip1%5==0) fpm = ly_phi[ip1/5];
                    else fpm = ly_phi[nid[1]/5];
                    if (ip2%5==0) fpp = ly_phi[ip2/5];
                    else fpp = ly_phi[nid[1]/5];
                } else {
                    if (nid[4]%5==0) fpm = ly_phi[nid[4]/5];
                    else fpm = phi;
                    if (nid[5]%5==0) fpp = ly_phi[nid[5]/5];
                    else fpp = phi;
                }
                dxz_phi = .25*(fmm+fpp-fmp-fpm);
            }

			trqq = trQQ(p);
			
            phi = ly_phi[id];
            phi2= phi*phi;
            ly_mu[id] = A_phi*(0.5*phi - 1.5*phi2 + phi*phi2) - ly_k*(d2x_phi+d2y_phi+d2z_phi);
            if (Q_on!=0) {
                ly_mu[id]-= A_ldg*ldg_b*S_lc*S_lc*trqq;
                if (ly_an!=0) {
                    ly_mu[id]-= 2.*ly_an*( d2x_phi*(Q[iq]+third*S_lc)+d2y_phi*(Q[iq+3]+third*S_lc)+d2z_phi*(-Q[iq]-Q[iq+3]+third*S_lc)+2.*(dxy_phi*Q[iq+1]+dxz_phi*Q[iq+2]+dyz_phi*Q[iq+4]) );
                    ly_mu[id]-= 2.*ly_an*( dx_phi*dxQ[0]+dy_phi*dyQ[3]-dz_phi*(dzQ[0]+dzQ[3])+dx_phi*(dyQ[1]+dzQ[2])+dy_phi*(dxQ[1]+dzQ[4])+dz_phi*(dxQ[2]+dyQ[4]) );
                }
            }
		}
	}

	MPI_Win_fence(0,winly_mu);
}



void evol_phi()
{
    int id, ib, ii, ip1, ip2, nid[6], flag;
    double dx_phi, dy_phi, dz_phi, d2x_mu, d2y_mu, d2z_mu, phi, phi2, phi_diffi, e_lyi;

	if (t_current%t_print==0) {
        e_lyi    = 0;
		phi_diffi= 0;
		flag     = 1;
	}

	for (id=0; id<lpoint; id++){
		if (info[id]==-1) {
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[id*6+ii];
			}
            dx_phi = 0.;
            dy_phi = 0.;
            dz_phi = 0.;
            d2x_mu = 0.;
            d2y_mu = 0.;
            d2z_mu = 0.;
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0) {
                dx_phi = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2x_mu = ly_mu[ip2/5] - 2.*ly_mu[id] + ly_mu[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dx_phi = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2x_mu = ly_mu[ip2/5] - ly_mu[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dx_phi = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2x_mu = ly_mu[ip1/5] - ly_mu[id];
            } else {
                printf("wrong!\n");
			}

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0) {
                dy_phi = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2y_mu = ly_mu[ip2/5] - 2.*ly_mu[id] + ly_mu[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dy_phi = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2y_mu = ly_mu[ip2/5] - ly_mu[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dy_phi = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2y_mu = ly_mu[ip1/5] - ly_mu[id];
            } else {
                printf("wrong!\n");
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0) {
                dz_phi = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
                d2z_mu = ly_mu[ip2/5] - 2.*ly_mu[id] + ly_mu[ip1/5];
			} else if (ip1%5!=0 && ip2%5==0) {
                dz_phi = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
                d2z_mu = ly_mu[ip2/5] - ly_mu[id];
			} else if (ip1%5==0 && ip2%5!=0) {
                dz_phi = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
                d2z_mu = ly_mu[ip1/5] - ly_mu[id];
            } else {
                printf("wrong!\n");
			}

            ly_dphi[id] = ly_Gamma*(d2x_mu+d2y_mu+d2z_mu);
//            temp += d2x_mu+d2y_mu+d2z_mu;
            if (flow_on!=0) ly_dphi[id]+= -u[id*3]*dx_phi-u[id*3+1]*dy_phi-u[id*3+2]*dz_phi; 

            if (flag==1) {
                phi    = ly_phi[id];
                phi2   = phi*phi;
                e_lyi += A_phi*(0.25*phi2*(1.-2.*phi+phi2)) + 0.5*ly_k*(dx_phi*dx_phi+dy_phi*dy_phi+dz_phi*dz_phi);
            }
        }
    }

	for (id=0; id<lpoint; id++){
		if (info[id]==-1) {
            ly_phi[id] += pdt * ly_dphi[id];
            if (flag==1) phi_diffi += ly_dphi[id]*ly_dphi[id];
        }
    }
	MPI_Win_fence(0,winly_phi);

	if(flag==1) MPI_Reduce(&phi_diffi, &phi_diff, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
	if(flag==1) MPI_Reduce(&e_lyi, &e_ly, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
}




void cal_stress_phi(int update)
{
    int id, ii, ip1, ip2, nid[6];
    double mu, dx_phi, dy_phi, dz_phi;

	for (id=0; id<lpoint; id++){
		if (info[id]==-1) {
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[id*6+ii];
			}

            dx_phi  = 0.;
            dy_phi  = 0.;
            dz_phi  = 0.;
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0) {
                dx_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
			} else if (ip1%5!=0 && ip2%5==0) {
                dx_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
			} else if (ip1%5==0 && ip2%5!=0) {
                dx_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
			}

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0) {
                dy_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
			} else if (ip1%5!=0 && ip2%5==0) {
                dy_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
			} else if (ip1%5==0 && ip2%5!=0) {
                dy_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
			}
			
			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0) {
                dz_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[ip1/5]);
			} else if (ip1%5!=0 && ip2%5==0) {
                dz_phi  = 0.5*(ly_phi[ip2/5] - ly_phi[id]);
			} else if (ip1%5==0 && ip2%5!=0) {
                dz_phi  = 0.5*(ly_phi[id] - ly_phi[ip1/5]);
			}

            mu = ly_mu[id];
			
            if (update!=1) {
                sigma_p[id*3]   = dx_phi*mu;
                sigma_p[id*3+1] = dy_phi*mu;
                sigma_p[id*3+2] = dz_phi*mu;
            } else {
                sigma_p[id*3]  += dx_phi*mu;
                sigma_p[id*3+1]+= dy_phi*mu;
                sigma_p[id*3+2]+= dz_phi*mu;
            }
		}
	}

	MPI_Win_fence(0, winsigma_p);
	MPI_Barrier(MPI_COMM_WORLD);
}
