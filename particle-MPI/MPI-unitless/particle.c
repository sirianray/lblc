/*
 *  particle.c
 *  
 *
 *  Created by Sirius on 3/3/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */

#include "particle.h"
#include "main.h"

real Lx, Ly, Lz, hLx, hLy, hLz;
double *p_m, *p_I, *p_rad, *p_pos, *p_vel, *p_acc, *p_qua, *p_angv, *p_torq;
double *p_Wtype, *p_W, *ubounce;
int *plink, *link_plist;

void p_allocate()
{
	if (npar>0) {
		p_m   = malloc(npar*sizeof(double));
		p_I   = malloc(npar*sizeof(double));
		p_rad = malloc(npar*sizeof(double));
		p_pos = malloc(3*npar*sizeof(double));
		p_vel = malloc(3*npar*sizeof(double));
		p_acc = malloc(3*npar*sizeof(double));
		p_qua = malloc(4*npar*sizeof(double));
		p_angv= malloc(3*npar*sizeof(double));
		p_torq= malloc(3*npar*sizeof(double));
		
		p_Wtype= malloc(  npar*sizeof(double));
		p_W    = malloc(  npar*sizeof(double));
		ubounce=malloc(15*point*sizeof(double));
	}
	
}

void p_deallocate()
{
	if (npar>0) {
		free(p_m);
		free(p_I);
		free(p_pos);
		free(p_vel);
		free(p_acc);
		free(p_qua);
		free(p_angv);
		free(p_torq);
		free(p_Wtype);
		free(p_W);
		free(ubounce);
	}
}

void p_init()
{
	int ipar;
	FILE *pos;
	
	Lx = (double)Nx;
	Ly = (double)Ny;
	Lz = (double)Nz;
	if (wall_on!=0) Lz = (double)(Nz-1);
	hLx= 0.5*Lx;
	hLy= 0.5*Ly;
	hLz= 0.5*Lz;

	if (npar>0) pos = fopen("ppos.in","r");
	
	for (ipar=0; ipar<npar; ipar++) {
		fscanf(pos,"%le %le %le %le %le %le\n",&p_pos[ipar*3],&p_pos[ipar*3+1],&p_pos[ipar*3+2],&p_rad[ipar],&p_Wtype[ipar],&p_W[ipar]);
		p_m[ipar]       = 100.;
		p_I[ipar]	= 100.;
		p_vel[ipar*3]   = 0;
		p_vel[ipar*3+1] = 0;
		p_vel[ipar*3+2] = 0;
		p_acc[ipar*3]   = 0;
		p_acc[ipar*3+1] = 0;
		p_acc[ipar*3+2] = 0;
		p_angv[ipar*3]  = 0;
		p_angv[ipar*3+1]= 0;
		p_angv[ipar*3+2]= 0;
		p_torq[ipar*3]  = 0;
		p_torq[ipar*3+1]= 0;
		p_torq[ipar*3+2]= 0;
		p_qua[ipar*4]   = 1.0;
		p_qua[ipar*4+1] = 0;
		p_qua[ipar*4+2] = 0;
		p_qua[ipar*4+3] = 0;
	}

	if (npar>0) fclose(pos);
}

int onsurface(double x, double y, double z, int ipar)
{
	double dx, dy, dz, r2p=p_rad[ipar]*p_rad[ipar], r2;
	
	dx = x - p_pos[ipar*3];
	dy = y - p_pos[ipar*3+1];
	dz = z - p_pos[ipar*3+2];
	if (dx<-hLx) dx+=Lx;
	if (dx> hLx) dx-=Lx;
	if (dy<-hLy) dy+=Ly;
	if (dy> hLy) dy-=Ly;
	if (wall_on==0) {
		if (dz<-hLz) dz+=Lz;
		if (dz> hLz) dz-=Lz;
	}
	r2 = dx*dx + dy*dy + dz*dz;
	if (r2>r2p) {
		return -1;
	} else if (r2<r2p) {
		return 1;
	} else {
		return 0;
	}
}

void p_iden()
{
	int ipar, xlo, xhi, ylo, yhi, zlo, zhi, x1, y1, z1, x2, y2, z2, i1, j1, k1, i2, j2, k2, i3, j3, k3;
	int i, j, k, ii, id, id1, id2, id3, di, p1, p2, nlink=0;
	int si, sj, sk, sid, sii, sjj;
	double x0, y0, z0, r, x, y, z, ubx, uby, ubz, omegax, omegay, omegaz, ub_dot_e, q[6];

	for (ipar=0; ipar<npar; ipar++) {
		x0  = p_pos[ipar*3];	// C.O.M of ipar'th particle
		y0  = p_pos[ipar*3+1];
		z0  = p_pos[ipar*3+2];
		r   = p_rad[ipar];
//		boundaries of search box for particle # ipar
		xlo = (int)(x0-r-1.01);
		xhi = (int)(x0+r+1.01);
		ylo = (int)(y0-r-1.01);
		yhi = (int)(y0+r+1.01);
		zlo = (int)(z0-r-1.01);
		zhi = (int)(z0+r+1.01);
		if (wall_on!=0) {
			if(zlo <  0)  zlo = 0;
			if(zhi >= Nz) zhi = Nz-1;
		}
		
		for (z1=zlo; z1<=zhi; z1++) {
			for (y1=ylo; y1<=yhi; y1++) {
				for (x1=xlo; x1<=xhi; x1++) {
					p1=onsurface(x1, y1, z1, ipar);
					i1 = x1;
					if(x1<0)   i1 += Nx;
					if(x1>=Nx) i1 -= Nx;
					j1 = y1;
					if(y1<0)   j1 += Ny;
					if(y1>=Ny) j1 -= Ny;
					k1 = z1;
					if (wall_on==0) {
						if(k1<0)   k1 += Nz;
						if(k1>=Nz) k1 -= Nz;
					}
					id1 = i1 + (j1+k1*Ny)*Nx;
					//	inside particle
					if (p1>0 && Q_on!=0 && id1/point==myid) info[id1%point]=-ipar-2;
					
					for (j=2; j<=10; j++) {
						z2 = z1 + e[j][2];
						if ( j!=3 && j!=5 && ( wall_on==0 || (z2>=0 && z2<Nz) ) ) {
							x2 = x1 + e[j][0];
							y2 = y1 + e[j][1];
							p2=onsurface(x2, y2, z2, ipar);

							if (p1 * p2 < 0) {
								i2 = x2;
								if(i2<0)   i2 += Nx;
								if(i2>=Nx) i2 -= Nx;
								j2 = y2;
								if(j2<0)   j2 += Ny;
								if(j2>=Ny) j2 -= Ny;
								k2 = z2;
								if (wall_on==0) {
									if(k2<0)   k2 += Nz;
									if(k2>=Nz) k2 -= Nz;
								}
								k   = bounce(j);
								id2 = i2 + (j2+k2*Ny)*Nx;
								
								x   = (double)x1 + 0.5*e[j][0];
								y   = (double)y1 + 0.5*e[j][1];
								z   = (double)z1 + 0.5*e[j][2];
								x  -= x0;
								y  -= y0;
								z  -= z0;
								
								if (flow_on!=0) {			//	modify streaming
									omegax = p_angv[ipar*3];
									omegay = p_angv[ipar*3+1];
									omegaz = p_angv[ipar*3+2];
									ubx    = p_vel[ipar*3]   + omegay*z - omegaz*y;
									uby    = p_vel[ipar*3+1] + omegaz*x - omegax*z;
									ubz    = p_vel[ipar*3+2] + omegax*y - omegay*x;
									ub_dot_e = ubx*e[j][0] + uby*e[j][1] + ubz*e[j][2];
									if(id1/point==myid) {
										nextf[(id1%point)*15+j]=(id1%point)*15+k;
										ubounce[(id1%point)*15+j] = ub_dot_e;
									}
									if(id2/point==myid) {
										nextf[(id2%point)*15+k]=(id2%point)*15+j;
										ubounce[(id2%point)*15+k] =-ub_dot_e;
									}
									nlink++;
									if (j<7) {
//										vsurf[nsurf*3]  = ubx;
//										vsurf[nsurf*3+1]= uby;
//										vsurf[nsurf*3+2]= ubz;
									}
								}
								
								if (j<7 && Q_on!=0) {			//	create surface point
									if (nsurf/node==myid) {
										id        = 10*(nsurf%node);
										surf[id]  = p_Wtype[ipar];
										surf[id+1]= p_W[ipar];
										normalize(&x,&y,&z);
										surf[id+2]= x;
										surf[id+3]= y;
										surf[id+4]= z;
										ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
										for (ii=0; ii<5; ii++) {
											surf[id+5+ii]  = q[ii];
										}
									}
									
									//	look for the nearest bulk point for the surface point
									if (p1<0) {
										si = i1;
										sj = j1;
										sk = k1;
										sid= id1;
										sii= j;
										sjj= k;
									} else {
										si = i2;
										sj = j2;
										sk = k2;
										sid= id2;
										sii= k;
										sjj= j;
									}
									if(sid/point==myid) neighb[(sid%point)*6+sii-1]=(nsurf-myid*node)*5-1;

									//	do the first pair neighbors
									i3 = si + e[sjj][0];
									j3 = sj + e[sjj][1];
									k3 = sk + e[sjj][2];
									if(i3<0)   i3 += Nx;
									if(i3>=Nx) i3 -= Nx;
									if(j3<0)   j3 += Ny;
									if(j3>=Ny) j3 -= Ny;
									if (wall_on==0) {
										if(k3<0)   k3 += Nz;
										if(k3>=Nz) k3 -= Nz;
									}
									id3 = i3 +(j3+k3*Ny)*Nx;
									if(nsurf/node==myid) {
										if (sii<sjj) {
											neighbsurf[(nsurf%node)*6+sii-1]=(sid-myid*point)*5;
											neighbsurf[(nsurf%node)*6+sjj-1]=(id3-myid*point)*5;
										} else {
											neighbsurf[(nsurf%node)*6+sjj-1]=(sid-myid*point)*5;
											neighbsurf[(nsurf%node)*6+sii-1]=(id3-myid*point)*5;
										}
									
										// find neighbors for surface point
										if (sii!=1 && sjj!=1) {
											if (x>0) di= 1;
											if (x<0) di=-1;
											i3 = si - di;
											j3 = sj;
											k3 = sk;
											if(i3<0)   i3 += Nx;
											if(i3>=Nx) i3 -= Nx;
											
											if (onsurface(i3,j3,k3,ipar)==-1) {
												neighbsurf[(nsurf%node)*6]  = 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-1;
												i3 = si + di;
												if(i3<0)   i3 += Nx;
												if(i3>=Nx) i3 -= Nx;
												neighbsurf[(nsurf%node)*6+1]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
											} else {
												i3 = si + di;
												if(i3<0)   i3 += Nx;
												if(i3>=Nx) i3 -= Nx;
												neighbsurf[(nsurf%node)*6]=   5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
												i3 = si + 2*di;
												if(i3<0)   i3 += Nx;
												if(i3>=Nx) i3 -= Nx;
												neighbsurf[(nsurf%node)*6+1]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
											}
										}
										
										if (sii!=3 && sjj!=3) {
											if (y>0) di= 1;
											if (y<0) di=-1;
											i3 = si;
											j3 = sj - di;
											k3 = sk;
											if(j3<0)   j3 += Ny;
											if(j3>=Ny) j3 -= Ny;
											
											if (onsurface(i3,j3,k3,ipar)==-1) {
												neighbsurf[(nsurf%node)*6+2]= 5*(i3 +(j3+k3*Ny)*Nx -myid*point)-1;
												j3 = sj + di;
												if(j3<0)   j3 += Ny;
												if(j3>=Nx) j3 -= Ny;
												neighbsurf[(nsurf%node)*6+3]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
											} else {
												j3 = sj + di;
												if(j3<0)   j3 += Ny;
												if(j3>=Nx) j3 -= Ny;
												neighbsurf[(nsurf%node)*6+2]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
												j3 = sj + 2*di;
												if(j3<0)   j3 += Ny;
												if(j3>=Nx) j3 -= Ny;
												neighbsurf[(nsurf%node)*6+3]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
											}
										}

										if (sii!=5 && sjj!=5) {
											if (z>0) di= 1;
											if (z<0) di=-1;
											i3 = si;
											j3 = sj;
											k3 = sk - di;
											if (wall_on==0) {
												if(k3<0)   k3 += Nz;
												if(k3>=Nz) k3 -= Nz;
											}
											
											if (onsurface(i3,j3,k3,ipar)==-1) {
												neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-1;
												k3 = sk + di;
												if (wall_on==0) {
													if(k3<0)   k3 += Nz;
													if(k3>=Nz) k3 -= Nz;
												}
												neighbsurf[(nsurf%node)*6+5]= 5*(i3 +(j3+k3*Ny)*Nx -myid*point);
											} else {
												k3 = sk + di;
												if (wall_on==0) {
													if(k3<0)   k3 += Nz;
													if(k3>=Nz) k3 -= Nz;
												}
												neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
												k3 = sk + 2*di;
												if (wall_on==0) {
													if(k3<0)   k3 += Nz;
													if(k3>=Nz) k3 -= Nz;
												}
												neighbsurf[(nsurf%node)*6+5]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
											}
										}
									}
									nsurf++;
								}
							}
						}
					}
				}
			}
		}
	}
	if (Q_on!=0 && (npar>0 || wall_on!=0) ) {
		MPI_Win_fence(0, winQsurf);
		MPI_Win_fence(0, winsurf);
		MPI_Win_fence(0, winneighbsurf);
	}
	if(flow_on!=0) MPI_Win_fence(0, winnf);
	MPI_Barrier(MPI_COMM_WORLD);
}
