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

int *p_shape;
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
		p_shape= malloc(  npar*sizeof(int));
		ubounce=malloc(15*point*sizeof(double));

		if (myid==root) printf("Particle related variables allocation done\n");
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
		free(p_shape);
		free(ubounce);
	}
}

void p_init()
{
	int ipar;
	FILE *pos;
	char ishape;
	
	Lx = (double)Nx;
	Ly = (double)Ny;
	Lz = (double)Nz;
	if (wall_z!=0) Lz = (double)(Nz-1);
	hLx= 0.5*Lx;
	hLy= 0.5*Ly;
	hLz= 0.5*Lz;

	if (npar>0) pos = fopen("ppos.in","r");
	
	for (ipar=0; ipar<npar; ipar++) {
		fscanf(pos,"%le %le %le %le %le %le %c\n",&p_pos[ipar*3],&p_pos[ipar*3+1],&p_pos[ipar*3+2],&p_rad[ipar],&p_Wtype[ipar],&p_W[ipar],&ishape);
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

		switch(ishape){
			case 'q':
				p_shape[ipar] = 2;
				break;
			case 's' :
				p_shape[ipar] = 0;
				break;
			case 'c' :
				p_shape[ipar] = 1;
				break;
			default:
				p_shape[ipar] = 0;
		}
	}

	if (npar>0) fclose(pos);
}

int onsurface_square(double x, double y, double z, int ipar, double *nx, double *ny, double *nz)
{
	double cx, cy, cz, dy, dz, ady, adz;

	cx = hLx - 0.5;
	cy = hLy - 0.5;
	cz = hLz - 0.5;

	dy = y - cy;
	dz = z - cz;
	ady= fabs(dy);
	adz= fabs(dz);

	*nx= 0;
	if (ady<adz) {
		*ny= 0;
		if (dz>0) {
			*nz =-1;
		} else {
			*nz = 1;
		}
	} else {
		*nz = 0;
		if (dy>0) {
			*ny =-1;
		} else {
			*ny = 1;
		}
	}

	if (ady<=p_rad[ipar] && adz<=p_rad[ipar]) {
		return -1;
	} else {
		return 1;
	}

}



int onsurface(double x, double y, double z, int ipar, double *nx, double *ny, double *nz)
{
	switch(p_shape[ipar]) {

   		case 1 :
		return onsurface_cylinder(x, y, z, ipar, nx, ny, nz);
      		break; 
	
		case 2: 
		return onsurface_square(x, y, z, ipar, nx, ny, nz);
      		break; 
  
   		default: 
		return onsurface_sphere(x, y, z, ipar, nx, ny, nz);
	}
}



int onsurface_cylinder(double x, double y, double z, int ipar, double *nx, double *ny, double *nz)
{
	double cx, cy, cz, dx, dy, dz, dr, as, r0;

	cx = hLx - 0.5;
	cy = hLy - 0.5;
	cz = hLz ;

	as = 1.;
	r0 = p_rad[ipar];
	if (p_rad[ipar]<0) {
		as = -1.;
		r0 = -r0;
	}

	if (p_pos[ipar*3]!=0) {
		dy = y - cy;
		dz = z - cz;
		dr = sqrt(dy*dy+dz*dz);
		*nx= 0;
		*ny= as*dy/dr;
		*nz= as*dz/dr;
		if ( (dr<=r0 && as<0) || (dr>r0 && as>0) ) {
			return -1;
		} else {
			return 1;
		}
	} else if (p_pos[ipar*3+1]!=0) {
                dx = x - cx;
                dz = z - cz;
                dr = sqrt(dx*dx+dz*dz);
		*nx= as*dx/dr;
		*ny= 0;
		*nz= as*dz/dr;
		if ( (dr<=r0 && as<0) || (dr>r0 && as>0) ) {
			return -1;
		} else {
			return 1;
		}
        } else if (p_pos[ipar*3+2]!=0) {
                dy = y - cy;
                dx = x - cx;
                dr = sqrt(dy*dy+dx*dx);
		*nx= as*dx/dr;
		*ny= as*dy/dr;
		*nz= 0;
		if ( (dr<=r0 && as<0) || (dr>r0 && as>0) ) {
			return -1;
		} else {
			return 1;
		}
        } else {
		*nx = 0;
		*ny = 0;
		*nz = 0;
		return -1;
	}
}

int onsurface_sphere(double x, double y, double z, int ipar, double *nx, double *ny, double *nz)
{
	double cx, cy, cz, dx, dy, dz, dr;

	cx = hLx - 0.5;
	cy = hLy - 0.5;
	cz = hLz ;

	dx = x - cx;
	dy = y - cy;
	dz = z - cz;
	dr = sqrt(dx*dx+dy*dy+dz*dz);
	if (p_rad[ipar]>0) {
		if (dr<p_rad[ipar]) {
			if (dr>0) {
				*nx = dx/dr;
				*ny = dy/dr;
				*nz = dz/dr;
			}
			return 1;
		} else {
			*nx = dx/dr;
			*ny = dy/dr;
			*nz = dz/dr;
			return -1;
		}
	} else if (p_rad[ipar]<0) {
		if (dr<-p_rad[ipar]) {
			if (dr>0) {
				*nx = -dx/dr;
				*ny = -dy/dr;
				*nz = -dz/dr;
			} else {
				*nx = 0;
				*ny = 0;
				*nz = 0;
			}
			return -1;
		} else {
			*nx = -dx/dr;
			*ny = -dy/dr;
			*nz = -dz/dr;
			return 1;
		}
	} else {
		*nx = 0;
		*ny = 0;
		*nz = 0;
		return 0;
	}

}

void p_iden()
{
	int ipar, xlo, xhi, ylo, yhi, zlo, zhi, x1, y1, z1, x2, y2, z2, i1, j1, k1, i2, j2, k2, i3, j3, k3;
	int i, j, k, ii, id, id1, id2, id3, di, p1, p2, nlink=0;
	int si, sj, sk, sid, sii, sjj, junk;
	double x, y, z, ubx, uby, ubz, omegax, omegay, omegaz, ub_dot_e, q[6], as, nx, ny, nz, junk1, junk2, junk3;

	for (ipar=0; ipar<npar; ipar++) {
		xlo = 0;
		xhi = Nx-1;
		ylo = 0;
		yhi = Ny-1;
		zlo = 0;
		zhi = Nz-1;
		
		for (z1=zlo; z1<=zhi; z1++) {
			for (y1=ylo; y1<=yhi; y1++) {
				for (x1=xlo; x1<=xhi; x1++) {
					p1=onsurface(x1, y1, z1, ipar, &junk1, &junk2, &junk3);
					i1 = x1;
					if(x1<0)   i1 += Nx;
					if(x1>=Nx) i1 -= Nx;
					j1 = y1;
					if(y1<0)   j1 += Ny;
					if(y1>=Ny) j1 -= Ny;
					k1 = z1;
					if (wall_z==0) {
						if(k1<0)   k1 += Nz;
						if(k1>=Nz) k1 -= Nz;
					}
					id1 = i1 + (j1+k1*Ny)*Nx;
					//	inside particle
					if (p1>0 && Q_on!=0 && id1/point==myid) info[id1%point]=-ipar-2;
					
					for (j=2; j<=10; j++) {
						z2 = z1 + e[j][2];
						if ( j!=3 && j!=5 && ( wall_z==0 || (z2>=0 && z2<Nz) ) ) {
							x2 = x1 + e[j][0];
							y2 = y1 + e[j][1];
							p2=onsurface(x2, y2, z2, ipar, &junk1, &junk2, &junk3);

							if (p1 * p2 < 0) {
								i2 = x2;
								if(i2<0)   i2 += Nx;
								if(i2>=Nx) i2 -= Nx;
								j2 = y2;
								if(j2<0)   j2 += Ny;
								if(j2>=Ny) j2 -= Ny;
								k2 = z2;
								if (wall_z==0) {
									if(k2<0)   k2 += Nz;
									if(k2>=Nz) k2 -= Nz;
								}
								k   = bounce(j);
								id2 = i2 + (j2+k2*Ny)*Nx;
								
								x   = (double)x1 + 0.5*e[j][0];
								y   = (double)y1 + 0.5*e[j][1];
								z   = (double)z1 + 0.5*e[j][2];
								
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
										if (rand_init==-5 || rand_init==-6) {
											double rr, rr2, irr, rrxy, omg, costhe, sinthe, cosphi, sinphi;
											nx = 0;
											ny = 0;
											nz = 0;
											rr2= x*x+y*y+z*z; 
											rr = sqrt(rr2);
											irr= 1.0/rr;
											omg= rr*q_ch;
											if (rand_init==-6) omg += atan2(-y,-x);
											rrxy = sqrt(x*x+y*y);
											if (rrxy<1e-18) {
												nz += 1.;
											} else {
												costhe =-z   *irr;
												sinthe = rrxy*irr;
												cosphi =-x/rrxy;
												sinphi =-y/rrxy;
												nx     = cos(omg)*costhe*cosphi - sin(omg)*sinphi;
												ny     = cos(omg)*costhe*sinphi + sin(omg)*cosphi;
												nz     =-cos(omg)*sinthe;
											}
										} else {
											junk=onsurface(x1+0.5*e[j][0],y1+0.5*e[j][1],z1+0.5*e[j][2],ipar,&nx,&ny,&nz);
										}
										normalize(&x,&y,&z);
										surf[id+2]= nx;
										surf[id+3]= ny;
										surf[id+4]= nz;
										ntoq(nx, ny, nz, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
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
									if (wall_z==0) {
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
											if (nx>0) di= 1;
											if (nx<0) di=-1;
											i3 = si - di;
											j3 = sj;
											k3 = sk;
											if(i3<0)   i3 += Nx;
											if(i3>=Nx) i3 -= Nx;
											
											if (onsurface(i3,j3,k3,ipar,&junk1,&junk2,&junk3)==-1) {
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
											if (ny>0) di= 1;
											if (ny<0) di=-1;
											i3 = si;
											j3 = sj - di;
											k3 = sk;
											if(j3<0)   j3 += Ny;
											if(j3>=Ny) j3 -= Ny;
											
											if (onsurface(i3,j3,k3,ipar,&junk1,&junk2,&junk3)==-1) {
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
											if (nz>0) di= 1;
											if (nz<0) di=-1;
											i3 = si;
											j3 = sj;
											k3 = sk - di;
											if (wall_z==0) {
												if(k3<0)   k3 += Nz;
												if(k3>=Nz) k3 -= Nz;
											}
											
											if (onsurface(i3,j3,k3,ipar,&junk1,&junk2,&junk3)==-1) {
												neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-1;
												k3 = sk + di;
												if (wall_z==0) {
													if(k3<0)   k3 += Nz;
													if(k3>=Nz) k3 -= Nz;
												}
												neighbsurf[(nsurf%node)*6+5]= 5*(i3 +(j3+k3*Ny)*Nx -myid*point);
											} else {
												k3 = sk + di;
												if (wall_z==0) {
													if(k3<0)   k3 += Nz;
													if(k3>=Nz) k3 -= Nz;
												}
												neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
												k3 = sk + 2*di;
												if (wall_z==0) {
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
	if (myid==root) printf("nsurf=%d\n",nsurf);
	if (Q_on!=0 && (npar>0 || wall_x!=0 || wall_y!=0 || wall_z!=0) ) {
		MPI_Win_fence(0, winQsurf);
		MPI_Win_fence(0, winsurf);
		MPI_Win_fence(0, winneighbsurf);
	}
	if(flow_on!=0) MPI_Win_fence(0, winnf);
	MPI_Barrier(MPI_COMM_WORLD);
}
