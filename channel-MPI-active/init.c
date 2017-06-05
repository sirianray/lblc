#include "main.h"
#include "particle.h"
#include <time.h>

void init()
{
	int i, j, k, ii, ip, id, ipar, jj;
	time_t t;
	double nx, ny, nz, q[6], sita0, sita, lambda1, lambda2;
	double dx, dy, dz, ir, ir2, ir3, x, y, z;

	if (myid==root) printf("Initializing..\n");
	
	if (newrun_on!=0) {
		if (myid==0) {
			sita0 = n_top[0]*n_bot[0]+n_top[1]*n_bot[1]+n_top[2]*n_bot[2];
			sita0 = sita0/sqrt(n_top[0]*n_top[0]+n_top[1]*n_top[1]+n_top[2]*n_top[2]);        
			sita0 = sita0/sqrt(n_bot[0]*n_bot[0]+n_bot[1]*n_bot[1]+n_bot[2]*n_bot[2]);
			sita0 = acos(sita0);
			
			if (rand_init==0) {			                                        // random random seed
				srand((unsigned)time(&t));
			} else if (rand_init==1) {		                                    // random seed from input
				srand((unsigned)rand_seed);
			}

			if ( (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) && Q_on!=0) {
				for (i=0; i<nsurf; i++) {
					for (j=0; j<5; j++) Qsurf[i*5+j] = surf[i*10+5+j];
				}
			}
			
			for (k=0; k<Nz; k++) {
				if (sita0<1e-2 && sita0>-1e-2) {
					lambda1 = 0.5;
					lambda2 = 0.5;
				} else {
					sita = sita0 * ((double)k+0.5)/((double)Nz);
					lambda1=(cos(sita)-cos(sita0)*cos(sita0-sita));
					lambda2=(cos(sita0-sita)-cos(sita)*cos(sita0));
				}
				nx = lambda1 * n_bot[0] + lambda2 * n_top[0];
				ny = lambda1 * n_bot[1] + lambda2 * n_top[1];
				nz = lambda1 * n_bot[2] + lambda2 * n_top[2];
				for (j=0; j<Ny; j++) {
					for (i=0; i<Nx; i++) {
						id = (i + (j+k*Ny)*Nx);
						if (rand_init==0 || rand_init==1) {                     // 0: random random seed, 1: random with rand_seed
							nx = randvec();
							ny = randvec();
							nz = randvec();
						}
                        else if (rand_init==2 || rand_init==3) {                // 2d system; 2: random random seed; 3: random with rand_seed
                            if (q_init<0) {                                     // directors in yz plane
							    nx = 0;
							    ny = randvec();
							    nz = randvec();
                            } else if (q_init==0) {                             // directors in xz plane
							    nx = randvec();
							    ny = 0;
							    nz = randvec();
                            } else {                                            // directors in xy plane
							    nx = randvec();
							    ny = randvec();
							    nz = 0;
                            }
                        }
						else if (rand_init==-2 && npar>0) {                     // dipolar ansatz
                            nx = n_bot[0];
                            ny = n_bot[1];
                            nz = n_bot[2];
                            for (ipar=0; ipar<npar; ipar++) {
                                    dx = (double)i-p_pos[ipar*3+0];
                                    dy = (double)j-p_pos[ipar*3+1];
                                    dz = (double)k-p_pos[ipar*3+2];
                                    ir2= 1./(dx*dx+dy*dy+dz*dz);
                                    ir = sqrt(ir2);
                                    ir3= ir*ir2;
                                    nx += q_init*p_rad[ipar]*p_rad[ipar]*dx*ir3;
                                    ny += q_init*p_rad[ipar]*p_rad[ipar]*dy*ir3;
                                    nz += q_init*p_rad[ipar]*p_rad[ipar]*dz*ir3;
                            }
                        }
						else if (rand_init==-3 || rand_init==-4) {              // -3: BPI; -4: BPII
							double A_init = q_init;
							double cst;
							double isq2 = 1.0 / sqrt(2);
							double sq2 = sqrt(2);
							if (rand_init == -3){
								cst = 2 * q_ch * 0.68;
							}
							else{
								cst = 2 * q_ch * 0.86;
							}
							if (rand_init == -3){
								x = i * cst * isq2;
								y = j * cst * isq2;
								z = k * cst * isq2;
								q[0] = A_init * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								q[3] = A_init * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								q[5] = A_init * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								q[1] = A_init * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								q[2] = A_init * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								q[4] = A_init * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}
							else if (rand_init ==-4){
								x    = i - Nx * 0.5;
								y    = j - Ny * 0.5;
								z    = k - Nz * 0.5;
								q[0] = A_init * (cos(cst * z) - cos(cst * y));
								q[1] = A_init * sin(cst * z);
								q[2] = A_init * sin(cst * y);
								q[3] = A_init * (cos(cst * x) - cos(cst * z));
								q[4] = A_init * sin(cst * x);
								q[5] = A_init * (cos(cst * y) - cos(cst * x));
							}
						} else if (rand_init==-5 || rand_init==-6) {            // -5: rotated BPI; -6: rotated BPII
							double r, r2, rxy, omg, costhe, sinthe, cosphi, sinphi;
							nx = 0;
							ny = 0;
							nz = 0;
							for (ipar=0; ipar<npar; ipar++) {
                                dx = (double)i-p_pos[ipar*3+0];
                                dy = (double)j-p_pos[ipar*3+1];
                                dz = (double)k-p_pos[ipar*3+2];
                                r2 = dx*dx+dy*dy+dz*dz;
                                r  = sqrt(r2);
								ir = 1.0/r;
								omg= r*q_ch;
								if (rand_init==-6) omg += atan2(dy,dx);
								rxy = sqrt(dx*dx+dy*dy);
								if (rxy<1e-18) {
									nz += 1.;
								} else {
									costhe = dz*ir;
									sinthe = rxy*ir;
									cosphi = dx/rxy;
									sinphi = dy/rxy;
									nx    += cos(omg)*costhe*cosphi - sin(omg)*sinphi;
									ny    += cos(omg)*costhe*sinphi + sin(omg)*cosphi;
									nz    +=-cos(omg)*sinthe;
								}
							}
						} else if ( (rand_init==-7 || rand_init==-8) && npar>0) {   // escaped-radial or planar-radial, for cylindrical symmetry
							double cx, cy, cz, th, dr;
							cx = hLx - 0.5;
							cy = hLy - 0.5;
							cz = hLz ;
							dx = (double)i - cx;
							dy = (double)j - cy;
							dz = (double)k - cz;

							if (p_pos[0]!=0) {
								dr = sqrt((dy*dy+dz*dz));
								th = dr/p_rad[0]*1.5708;
								nx = cos(th);
								ny = 0;
								nz = 0;
								if (rand_init==-8) nx=1;
								if (dr>0) {
									ny = sin(th)*dy/dr;
									nz = sin(th)*dz/dr;
									if (rand_init==-8) {
										nx = 0;
										ny = dy/dr;
										nz = dz/dr;
									}
								}
							} else if (p_pos[1]!=0) {
								dr = sqrt((dx*dx+dz*dz));
								th = dr/p_rad[0]*1.5708;
								ny = cos(th);
								nx = 0;
								nz = 0;
								if (rand_init==-8) ny=1;
								if (dr>0) {
									nx = sin(th)*dx/dr;
									nz = sin(th)*dz/dr;
									if (rand_init==-8) {
										nx = dx/dr;
										ny = 0;
										nz = dz/dr;
									}
								}
							} else if (p_pos[2]!=0) {
								dr = sqrt((dx*dx+dy*dy));
								th = dr/p_rad[0]*1.5708;
								nz = cos(th);
								nx = 0;
								ny = 0;
								if (rand_init==-8) nz=1;
								if (dr>0) {
									nx = sin(th)*dx/dr;
									ny = sin(th)*dy/dr;
									if (rand_init==-8) {
										nx = dx/dr;
										ny = dy/dr;
										nz = 0;
									}
								}
							} else {
								nx = 0;
								ny = 0;
								nz = 0;
							}
						} else if (rand_init==-9) {                             // radial for spherical symmetry
							double cx, cy, cz, th, dr, idr;
							cx = p_pos[0]; 
							cy = p_pos[1];
							cz = p_pos[2];
							dx = (double)i - cx;
							dy = (double)j - cy;
							dz = (double)k - cz;
							dr = sqrt(dx*dx+dy*dy+dz*dz);
							if (dr>0) {
								idr = 1./dr;
								nx  = dx*idr;
								ny  = dy*idr;
								nz  = dz*idr;
							} else {
								nx  = 0;
								ny  = 0;
								nz  = 0;
							}
						} else if (rand_init==-10) {                            // twist-planar for planar cylindrical system
                            double cx, cy, cz, dr, idr, beta0, beta;
							cx = hLx - 0.5;
							cy = hLy - 0.5;
							cz = hLz ;
							dx = (double)i - cx;
							dy = (double)j - cy;
							dz = (double)k - cz;
                            if (K24>=2.*K2 && K2>0 && K3>0) {
                                beta0 = sqrt(K24*(K24-2.*K2)/(K2*K3));
                                beta0 = atan(beta0);
                            } else {
                                beta0 = 0;
                            }
                            if (p_pos[0]!=0) {
                                dr  = sqrt(dy*dy+dz*dz);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    beta= beta0*dr/p_rad[0];
                                    nx  = cos(beta);
                                    ny  =-sin(beta)*dz*idr;
                                    nz  = sin(beta)*dy*idr;
                                } else {
                                    nx  = 1.;
                                    ny  = 0.;
                                    nz  = 0.;
                                }
                            } else if (p_pos[1]!=0) {
                                dr  = sqrt(dx*dx+dz*dz);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    beta= beta0*dr/p_rad[0];
                                    nx  = sin(beta)*dz*idr;
                                    ny  = cos(beta);
                                    nz  =-sin(beta)*dx*idr;
                                } else {
                                    nx  = 0;
                                    ny  = 1.;
                                    nz  = 0.;
                                }
                            } else if (p_pos[2]!=0) {
                                dr  = sqrt(dx*dx+dy*dy);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    beta= beta0*dr/p_rad[0];
                                    nx  =-sin(beta)*dy*idr;
                                    ny  = sin(beta)*dx*idr;
                                    nz  = cos(beta);
                                } else {
                                    nx  = 0.;
                                    ny  = 0.;
                                    nz  = 1.;
                                }
                            } else {
								nx = 0;
								ny = 0;
								nz = 0;
                            }
						} else if (rand_init==-11) {                            // concentric for cylindrical symmetry
							double cx, cy, cz, th, dr, idr;
							cx = hLx - 0.5;
							cy = hLy - 0.5;
							cz = hLz ;
							dx = (double)i - cx;
							dy = (double)j - cy;
							dz = (double)k - cz;
                            if (p_pos[0]>0) {
                                dr = sqrt(dy*dy+dz*dz);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    nx  = 0;
                                    ny  =-dz*idr;
                                    nz  = dy*idr;
                                } else {
                                    nx  = 0;
                                    ny  = 0;
                                    nz  = 0;
                                }
                            } else if (p_pos[1]>0) {
                                dr = sqrt(dx*dx+dz*dz);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    nx  = dz*idr;
                                    ny  = 0;
                                    nz  =-dx*idr;
                                } else {
                                    nx  = 0;
                                    ny  = 0;
                                    nz  = 0;
                                }
                            } else if (p_pos[2]>0) {
                                dr = sqrt(dx*dx+dy*dy);
                                if (dr>1e-3) {
                                    idr = 1./dr;
                                    nx  =-dy*idr;
                                    ny  = dx*idr;
                                    nz  = 0;
                                } else {
                                    nx  = 0;
                                    ny  = 0;
                                    nz  = 0;
                                }
                            } else {
                                nx  = 0;
                                ny  = 0;
                                nz  = 0;
                            }
                        } else if (rand_init==-12 && npar>0) {      // escape radial with defects
							double cx, cy, cz, d, dr, idr, nr;
							cx = hLx - 0.5;
							cy = hLy - 0.5;
							cz = hLz ;

							dx = (double)i - cx;
							dy = (double)j - cy;
							dz = (double)k - cz;
                            if (p_pos[0]>0) {
                                d  = 0.75*Lx;
                                dx+= 0.5*d;
                                if (dx>d) dx-=Lx;
                                dr = sqrt(dy*dy+dz*dz);
                                if (dr>1e-3) {
                                    nx = dx*(d-dx)*(dx+Nx-d)*(fabs(p_rad[0])-dr);
                                    nr = fabs(p_rad[0])*d*(Nx-d)*dr;
                                    idr= 1./sqrt(nx*nx+nr*nr);
                                    nx = nx*idr;
                                    nr = nr*idr;
                                    idr= 1./dr;
                                    ny = nr*dy*idr;
                                    nz = nr*dz*idr;
                                } else {
                                    nx = 1.;
                                    ny = 0.;
                                    nz = 0.;
                                }
                            } else if (p_pos[1]>0) {
                                d  = 0.75*Ly;
                                dy+= 0.5*d;
                                if (dy>d) dy-=Ly;
                                dr = sqrt(dx*dx+dz*dz);
                                if (dr>1e-3) {
                                    ny = dy*(d-dy)*(dy+Ny-d)*(fabs(p_rad[0])-dr);
                                    nr = fabs(p_rad[0])*d*(Ny-d)*dr;
                                    idr= 1./sqrt(ny*ny+nr*nr);
                                    ny = ny*idr;
                                    nr = nr*idr;
                                    idr= 1./dr;
                                    nx = nr*dx*idr;
                                    nz = nr*dz*idr;
                                } else {
                                    nx = 0.;
                                    ny = 1.;
                                    nz = 0.;
                                }
                            } else if (p_pos[2]>0) {
                                d  = 0.75*Lz;
                                dz+= 0.5*d;
                                if (dz>d) dz-=Lz;
                                dr = sqrt(dx*dx+dy*dy);
                                if (dr>1e-3) {
                                    nz = dz*(d-dz)*(dz+Nz-d)*(fabs(p_rad[0])-dr);
                                    nr = fabs(p_rad[0])*d*(Nz-d)*dr;
                                    idr= 1./sqrt(nz*nz+nr*nr);
                                    nz = nz*idr;
                                    nr = nr*idr;
                                    idr= 1./dr;
                                    nx = nr*dx*idr;
                                    ny = nr*dy*idr;
                                } else {
                                    nx = 0.;
                                    ny = 0.;
                                    nz = 1.;
                                }
                            }
						} else if (rand_init==-13) {        // defect annihilation in channel
							if (i<Nx/5 || i>4*Nx/5) {
								nx = 1.;
								ny = 0.;
								nz = 0.;
							} else {
								nx = 0.;
								ny = 0.;
								nz = 1.;
							}
                        } else if (rand_init==-14) {        // defect annihilation in capillary
                            if (k<Nz/5 || k>4*Nz/5) {
                                nx = 0;
                                ny = 0;
                                nz = 1.;
                            } else {
                                dx = (double)i - hLx + 0.5;
                                dy = (double)j - hLy + 0.5;
                                ir = sqrt(dx*dx+dy*dy);
                                if (ir<1e-5) {
                                    nx = 0;
                                    ny = 0;
                                    nz = 1;
                                } else {
                                    ir = 1./ir;
                                    nx = dx*ir;
                                    ny = dy*ir;
                                    nz = 0;
                                }
                            }
                       } else if (rand_init==-31) {        // defect annihilation, defect line perpendicular to far-field director
                           double cx, cy, b, th;
                           cx = (double)Nx*0.5;
                           cy = (double)Ny*0.5;
                           b  = q_init;
                           if (b>0.5*cy) b = 0.5*cy;
                           if (i>cx-0.5*q_init && i<cx+0.5*q_init && j>cy-b && j<=cy) {
                                th = (cy-(double)j)/b;
                                nx = cos(th);
                                ny =-sin(th);
                                nz = 0;
                           } else if (i>cx-0.5*q_init && i<cx+0.5*q_init && j>cy && j<cy+b) {
                                th = ((double)j-cy)/b;
                                nx = cos(th);
                                ny = sin(th);
                                nz = 0;
                           } else {
                                nx = 0.;
                                ny = 1.;
                                nz = 0.;
                           }
                        } else if (rand_init==-32) {        // defect annihilation, defect line parallel to far-field director
                            double cx, cy, b, th;
                            cx = (double)Nx*0.5;
                            cy = (double)Ny*0.5;
                            b  = q_init;
                            if (b>0.5*cy) b = 0.5*cy;
                            if (i>cx-0.5*q_init && i<cx+0.5*q_init && j>cy-b && j<=cy) {
                                th = (cy-(double)j)/b;
                                nx = sin(th);
                                ny = cos(th);
                                nz = 0;
                            } else if (i>cx-0.5*q_init && i<cx+0.5*q_init && j>cy && j<cy+b) {
                                th = ((double)j-cy)/b;
                                nx =-sin(th);
                                ny = cos(th);
                                nz = 0;
                            } else {
                                nx = 1.;
                                ny = 0.;
                                nz = 0.;
                            }
                        } else if (rand_init==-40) {        // guess df from defect.in
                            double nx0, ny0, dx, dy, dr, th, sn, dot;
                            int ii;
                            for (ii=0; ii<ndefect; ii++) {
                                dx = i-d_x[ii];
                                dy = j-d_y[ii];
                                dr = sqrt(dx*dx+dy*dy);
                                if (dr==0) dr=0.001;
                                th = acos((dx*d_nx[ii]+dy*d_ny[ii])/dr);
                                sn = 1;
                                if (d_ch[ii]>0) sn = -1;
                                if (d_ny[ii]*dx-d_nx[ii]*dy<0) sn=-sn;
                                th=atan2(d_ny[ii],d_nx[ii])+0.5*(sn*th+PI);
                                
                                if (ii==0) {
                                    nx0 = cos(th)*exp(-dr*d_ixl);
                                    ny0 = sin(th)*exp(-dr*d_ixl);
                                    nx  = cos(th)*exp(-dr*d_ixl);
                                    ny  = sin(th)*exp(-dr*d_ixl);
                                } else {
                                    dot=nx0*cos(th)+ny0*sin(th);
                                    if (dot>=0) {
                                        nx += cos(th)*exp(-dr*d_ixl);
                                        ny += sin(th)*exp(-dr*d_ixl);
                                    } else {
                                        nx -= cos(th)*exp(-dr*d_ixl);
                                        ny -= sin(th)*exp(-dr*d_ixl);
                                    }
                                }
                            }
                            nz = 0;
                        }

						if (rand_init!=-3 && rand_init!=-4) ntoq(nx, ny, nz, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
						for (ii=0; ii<5; ii++) Q[id*5+ii] = q[ii];
						if (q[0]+q[3]+q[5]>1e-15 || q[0]+q[3]+q[5]<-1e-15 ) {
							printf("initial bulk Q not right\n");
							exit(-1);
						}
					}
				}
			}

			if (flow_on!=0) {
				for (k=0; k<Nz; k++) {
					for (j=0; j<Ny; j++) {
						for (i=0; i<Nx; i++) {
							id = i + (j+k*Ny)*Nx;
							Rho[id]  = rho;
							u[id*3]  = 0;
//							u[id*3+1]= uy_bottom + (uy_top-uy_bottom)*((double)k+0.5)/(double)Nz;
							u[id*3+1]= 0;
							u[id*3+2]= 0;
						}
					}
				}			
			}
		}		


        // assign surf Q to its nearest neighbor's Q
//       MPI_Barrier(MPI_COMM_WORLD);
//       MPI_Win_fence(0, winq);
//       if ( (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) && Q_on!=0) {
//           for (id=0; id<node; id++) {
//               if (myid*node+id<nsurf && surf[id*10]==2) {
//                   for (ii=0; ii<5; ii++) {
//                       ip = neighbsurf[id*6+ii];
//                       if (ip%5==0) {
//                           for (jj=0; jj<5; jj++) Qsurf[id*5+jj]=Q[ip+jj];
//                           break;
//                       }
//                   }
//               }
//           }
//       }
	} else {
		read_restart();
	}

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_fence(0, winq);
        MPI_Win_fence(0, winr);
        MPI_Win_fence(0, winu);
	if ( (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) && Q_on!=0 ) MPI_Win_fence(0, winQsurf);
}


void ntoq(double nx, double ny, double nz, double *q0, double *q1, double *q2, double *q3, double *q4, double *q5)
{
	double r2, r2i;
	
	r2 = nx*nx + ny*ny + nz*nz;
	if (debug_on!=0 && r2==0) {
		printf("warning: zero director encountered in function ntoq\n");
		*q0 = 0;
		*q1 = 0;
		*q2 = 0;
		*q3 = 0;
		*q4 = 0;
		*q5 = 0;
	} else {
		r2i= 1.0/r2;	
		*q0 = S_lc * (nx*nx*r2i - third);
		*q1 = S_lc * (nx*ny*r2i);
		*q2 = S_lc * (nx*nz*r2i);
		*q3 = S_lc * (ny*ny*r2i - third);
		*q4 = S_lc * (ny*nz*r2i);
		*q5 = -(*q0) - (*q3);
	}
}



real randvec()
{
	double r;
	
    r =  1.0 - 2.0 * (double)rand()/(double)RAND_MAX ;
    return r;
}
