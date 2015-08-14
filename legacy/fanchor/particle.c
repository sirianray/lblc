/*
 *  particle.c
 *  
 *
 *  Created by Sirius on 4/4/14;
 *  Modified by Sirius on 6/8/2014;
 *  Copyright 2014 home. All rights reserved.
 *
 */

#include "particle.h"
int npar, nlink, nlist, p_move;
int link_point[Nx][Ny][Nz][NPMAX], link_listn[NLINKMAX][10], link_plist[NLINKMAX][4];
int link_pointo[Nx][Ny][Nz][NPMAX];
int link_point2[Nx][Ny][Nz][15];
char p_type[NPMAX];
real Lx, Ly, Lz;
//real **link_listf;
real p_m[NPMAX], p_I[NPMAX][2], p_rad[NPMAX], p_size[NPMAX], p_orien[NPMAX][3], p_orien2[NPMAX][3], p_norm[NPMAX][3];
real p_ncoil[NPMAX], p_b[NPMAX], p_thick[NPMAX], p_ph0[NPMAX], p_height[NPMAX], p_const1[NPMAX];
real p_pos[NPMAX][3], p_vel[NPMAX][3], p_acc[NPMAX][3], p_qua[NPMAX][4], p_angv[NPMAX][3], p_torq[NPMAX][3];

// for finite anchoring
int p_nsurf;
int p_surfi[Nx][Ny][Nz][3];
real p_surf[NLINKMAX][7], p_surfQ0[NLINKMAX][5], p_surfQ[NLINKMAX][5], p_surfdQ[NLINKMAX][5], p_surfW[NLINKMAX];
real p_anchor;

void p_init()
{
	int ipar, i, j, k, kp;
	real Lmin, temp;
	
	Lx = (real)Nx;
	Ly = (real)Ny;
	Lz = (real)Nz;
	if (Lx>Ly){
		Lmin = Ly;
	} else {
		Lmin = Lx;
	}
	if (Lmin>Lz)Lmin=Lz;
	
/*
		  p_type  p_m  p_ncoil  p_b  p_rad  p_thick  p_ph0  p_height  p_const1  p_size
 sphere        x    x        x           x                                           x
 helix         x    x        x    x      x        x      x         x         x       x
 wall          x
 (in developing)
 */
	npar=1;
	for (ipar=0; ipar<npar; ipar++) {
		p_type[ipar]  = 's';						// particle type: 's' - sphere; 'h' - helix
		p_m[ipar]     = 1.0;						// mass
		p_ncoil[ipar] = 5.0;						// number of helixs/coils
		p_b[ipar]     = 5.7;						// height of one circle of coil; +: right hand chiral; -: left hand chiral
		p_rad[ipar]   = 8.01;						// radius of the coil
		p_thick[ipar] = 0.99;						// radius of the string of the coil
		p_ph0[ipar]   = 2;							// initial phi in parameterized trajectory of the helix
		p_height[ipar]= p_ncoil[ipar] * p_b[ipar];	// height of the helix
		if (p_height[ipar]<0) p_height[ipar]=-p_height[ipar];
		
		temp = twopi * p_rad[ipar] * (1.0 / p_b[ipar]);
		p_const1[ipar]= 1.0 / (1.0 + temp*temp);	// constant for finding the nearest point to helix
		
		p_size[ipar] = p_rad[ipar];
		if (p_type[ipar]=='h') {					// helix size
			p_size[ipar] = sqrt(0.25 * p_height[ipar]* p_height[ipar] + p_rad[ipar]*p_rad[ipar]) + 4.0 * p_thick[ipar];
		}
		
		if(p_size[ipar]>=0.5*Lmin) {				// particle size should be suitable to the box
			printf("error: the size of helix #%d too large!",ipar);
			exit(-1);
		}
		
		p_move = 0;
		for (i=0; i<3; i++) {
			p_vel[ipar][i]  = 0;					// velocities			
			p_acc[ipar][i]  = 0;					// accelerations
			p_angv[ipar][i] = 0;					// angular velocities
			p_torq[ipar][i] = 0;					// torques
			p_move += (int)(p_vel[ipar][i]*p_vel[ipar][i] + p_angv[ipar][i]*p_angv[ipar][i]);
													// if particle is initially static, set p_move to 0
		}
	}
	
//	initialization of surface points for the finite anchoring case
	p_nsurf = 0;
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				p_surfi[i][j][k][0]=-1;
				p_surfi[i][j][k][1]=-1;
				p_surfi[i][j][k][2]=-1;
			}
		}
	}
	for (i=0; i<NLINKMAX; i++) {
		p_surf[i][3]=-1;		// default(inactive) surface point info list should have its direction of -1
	}
	
//	position & quternion
	p_pos[0][0] = 0.5 * (Lx);
	p_pos[0][1] = 0.5 * (Ly);
	p_pos[0][2] = 0.5 * (Lz-1.0);
//	p_qua[0][0] =-halfsqrt2;		// helix align in y direction
//	p_qua[0][1] = halfsqrt2;
	p_qua[0][0] = 1;			// helix align in z direction
	p_qua[0][1] = 0;
	p_qua[0][2] = 0;
	p_qua[0][3] = 0;
	
	for (ipar=0; ipar<npar; ipar++) {
		rotation( &p_orien[ipar][0], &p_orien[ipar][1], &p_orien[ipar][2],0,0,1,p_qua[ipar][0],p_qua[ipar][1],p_qua[ipar][2],p_qua[ipar][3]);
		rotation(&p_orien2[ipar][0],&p_orien2[ipar][1],&p_orien2[ipar][2],1,0,0,p_qua[ipar][0],p_qua[ipar][1],p_qua[ipar][2],p_qua[ipar][3]);
	}
	
// particle anchoring
	p_anchor = -1;
}

// identify boundary links and the associated points
void p_iden()
{
	int ipar, i, j, k, kk;
	int xlo, xhi, ylo, yhi, zlo, zhi;
	int x1, y1, z1, x2, y2, z2;
	int i1, i2, i3, j1, j2, j3, k1, k2, k3;
	int p1, p2;
	real x0, y0, z0, r, x3, y3, z3, nx, ny, nz;
	char itype;
	
	nlink = 0;
	nlist = 0;
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				for (ipar=0; ipar<npar; ipar++){
					link_point[i][j][k][ipar]  = 0;
				}
				for (ipar=0; ipar<15; ipar++){
					link_point2[i][j][k][ipar] = 0;
				}
			}
		}
	}
	
	for (ipar=0; ipar<npar; ipar++) {
		x0  = p_pos[ipar][0];	// C.O.M of ipar'th particle
		y0  = p_pos[ipar][1];
		z0  = p_pos[ipar][2];
		r   = p_size[ipar];
		itype  = p_type[ipar];
// boundaries of search box for particle # ipar
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
		for (x1=xlo; x1<=xhi; x1++) {
			for (y1=ylo; y1<=yhi; y1++) {
				for (z1=zlo; z1<=zhi; z1++) {
					p1=onsurface(x1, y1, z1, ipar, itype);
					i1 = x1;
					if(x1<0)   i1 += Nx;
					if(x1>=Nx) i1 -= Nx;
					j1 = y1;
					if(y1<0)   j1 += Ny;
					if(y1>=Ny) j1 -= Ny;
					k1 = z1;
					if (wall_on==0) {
						if(z1<0)   k1 += Nz;
						if(z1>=Nz) k1 -= Nz;
					}
					
					for (j=3; j<=10; j++) {
						z2 = z1 + e[j][2];
						if ( j!=6 && ( wall_on==0 || (z2>=0 && z2<Nz) ) ) {
							x2 = x1 + e[j][0];
							y2 = y1 + e[j][1];
							p2=onsurface(x2, y2, z2, ipar, itype);
							if (p1 * p2 < 0) {
								i2 = x2;
								if(x2<0)   i2 += Nx;
								if(x2>=Nx) i2 -= Nx;
								j2 = y2;
								if(y2<0)   j2 += Ny;
								if(y2>=Ny) j2 -= Ny;
								k2 = z2;
								if (wall_on==0) {
									if(z2<0)   k2 += Nz;
									if(z2>=Nz) k2 -= Nz;
								}
								
// link_listn: activity(effectiveness), ipar(which particle), i1, j1, k1,
//				f1(which velocity), i2, j2, k2, f2(which velocity)
								link_listn[nlink][0] = 1;
								link_listn[nlink][1] = ipar;
								k = bounce(j);
								
								// for finite anchoring, preparing for assignment
								if (j<=6) {
									x3 = 0.5 * (double)(x1 + x2);
									y3 = 0.5 * (double)(y1 + y2);
									z3 = 0.5 * (double)(z1 + z2);
									i3 = floor(x3);
									j3 = floor(y3);
									k3 = floor(z3);
									if (i3<0) i3+=Nx;
									if (j3<0) j3+=Ny;
									if (k3<0) k3+=Nz;
									// the following debug part will be removed in the future
									if (debug_on && !( (int)x3==i3 && (int)y3==j3 && (int)z3==k3 )) {
										printf("error in p_iden()!\n");
									}
									find_normal(&nx, &ny, &nz, x3, y3, z3, ipar, itype);
									p_surf[p_nsurf][4]   = nx;
									p_surf[p_nsurf][5]   = ny;
									p_surf[p_nsurf][6]   = nz;
									p_surfQ0[p_nsurf][0] = S * (nx*nx-third);
									p_surfQ0[p_nsurf][1] = S * (nx*ny);
									p_surfQ0[p_nsurf][2] = S * (nx*nz);
									p_surfQ0[p_nsurf][3] = S * (ny*ny-third);
									p_surfQ0[p_nsurf][4] = S * (ny*nz);
									p_surfW[p_nsurf]     = p_anchor;
									kk = axis_proj(j);
									p_surfi[i3][j3][k3][kk] = p_nsurf;
								}								
								
// if boundary "side in" direction is specified, the following if condition can be modified to account for that
								if (p1<0) {
									link_listn[nlink][2] = i1;
									link_listn[nlink][3] = j1;
									link_listn[nlink][4] = k1;
									link_listn[nlink][5] = j;
									link_listn[nlink][6] = i2;
									link_listn[nlink][7] = j2;
									link_listn[nlink][8] = k2;
									link_listn[nlink][9] = k;
									link_point2[i1][j1][k1][j] = -(ipar+1);
									link_point2[i2][j2][k2][k] = ipar+1;
									// for finite anchoring
									if (j<=6) {
										p_surf[p_nsurf][0]   = i2;
										p_surf[p_nsurf][1]   = j2;
										p_surf[p_nsurf][2]   = k2;
										p_surf[p_nsurf][3]   = k;
										p_nsurf++;
									}									
								}
								else {
									link_listn[nlink][6] = i1;
									link_listn[nlink][7] = j1;
									link_listn[nlink][8] = k1;
									link_listn[nlink][9] = j;
									link_listn[nlink][2] = i2;
									link_listn[nlink][3] = j2;
									link_listn[nlink][4] = k2;
									link_listn[nlink][5] = k;
									link_point2[i1][j1][k1][j] = ipar+1;
									link_point2[i2][j2][k2][k] = -(ipar+1);
									// for finite anchoring
									if (j<=6) {
										p_surf[p_nsurf][0]   = i1;
										p_surf[p_nsurf][1]   = j1;
										p_surf[p_nsurf][2]   = k1;
										p_surf[p_nsurf][3]   = j;
										p_nsurf++;
									}
								}
								nlink++;
								
								// build link point list
								if (link_point[i1][j1][k1][ipar]==0) {
									link_plist[nlist][0] = i1;
									link_plist[nlist][1] = j1;
									link_plist[nlist][2] = k1;
									link_plist[nlist][3] = ipar;
									nlist++;
									link_point[i1][j1][k1][ipar]  = p1;
								} else {
									link_point[i1][j1][k1][ipar] += p1;
								}
								
								if (link_point[i2][j2][k2][ipar]==0) {
									link_plist[nlist][0] = i2;
									link_plist[nlist][1] = j2;
									link_plist[nlist][2] = k2;
									link_plist[nlist][3] = ipar;
									nlist++;
									link_point[i2][j2][k2][ipar]  = p2;
								} else {
									link_point[i2][j2][k2][ipar] += p2;
								}
								
							}
						}
					}
				}
			}
		}
	}
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				for (ipar=0; ipar<npar; ipar++){
					link_pointo[i][j][k][ipar] = link_point[i][j][k][ipar];
				}
			}
		}
	}
	
	for (i=0; i<p_nsurf; i++) {
		for (j=0; j<5; j++) {
			p_surfQ[i][j] = p_surfQ0[i][j];
		}
	}
}

int onsurface(real x, real y, real z, int ipar, char surf)
{
	real dx, dy, dz, dr2, r2;
	real phi, th, x1, y1, z1, hx1, hx2, hy;
	real psi, b, r, phi0, x2, y2 ,z2, tp, ti, temp, tt, ni, tir, c1, c2, hn;
	
	if (ipar<0 || ipar>=npar) {
		printf("error in subroutine onsurface: particle number wrong\n");
		exit(-1);
	}
	
	if (surf=='s') {				// sphere
		r2 = p_rad[ipar] * p_rad[ipar];
		dx = x - p_pos[ipar][0];
		dy = y - p_pos[ipar][1];
		dz = z - p_pos[ipar][2];
		if(dx<-0.5*Lx)dx+=Lx;
		if(dx> 0.5*Lx)dx-=Lx;
		if(dy<-0.5*Ly)dy+=Ly;
		if(dy> 0.5*Ly)dy-=Ly;
		if (wall_on==0) {
			if(dz<-0.5*Lz)dz+=Lz;
			if(dz> 0.5*Lz)dz-=Lz;
		}
		dr2= dx*dx + dy*dy +dz*dz;
		if (dr2<r2) {
			return 1;
		} else if (dr2 >r2) {
			return -1;
		} else {
			return 0;
			printf("warning: lattice point on particle surface\n");
		}
	}
	
	if (surf=='h') {				// helix; For more details, see notes 6/9/2014
		th  = p_thick[ipar];
		r   = p_rad[ipar];
		hy = 0.5*p_height[ipar]+th;	// half of helix real length
		hx1= r - th;				// bounds of the helix real radius (hx1: inner; hx2: outer)
		hx2= r + th;
		x1 = x - p_pos[ipar][0];	// vector relative to helix C.O.M
		y1 = y - p_pos[ipar][1];
		z1 = z - p_pos[ipar][2];
		if(x1<-0.5*Lx)x1+=Lx;
		if(x1> 0.5*Lx)x1-=Lx;
		if(y1<-0.5*Ly)y1+=Ly;
		if(y1> 0.5*Ly)y1-=Ly;
		if (wall_on==0) {
			if(z1<-0.5*Lz)z1+=Lz;
			if(z1> 0.5*Lz)z1-=Lz;
		}
		rotationback(&x2, &y2, &z2, x1, y1, z1, p_qua[ipar][0], p_qua[ipar][1], p_qua[ipar][2], p_qua[ipar][3]);
		// rotate back to helix reference
		
		if ( z2 > hy || z2 < -hy ) {
			return -1;
		}
		
		dr2 = x2*x2 + y2*y2;
		if (dr2<hx1*hx1 || dr2>hx2*hx2) {
			return -1;
		}
		
//		if the point is inside the cylindrical shell, do the following
		b   = p_b[ipar];
		phi0= p_ph0[ipar];
		hn  = 0.5 * p_ncoil[ipar];
		c2  = p_const1[ipar];		// coefficients for finding nearest point
		c1  = 1.0 - c2;
		tp   = z2/b;
		ni   = floor(tp);			// which coil the particle is within
		psi  = atan2(y2, x2) - phi0;// relative angle
		tt   = psi / twopi;
		tt   = tt - floor(tt);
		temp = tp - ni - tt;		// to tell which coil the particle is nearer to
		if ( temp > 0.5 ) {
			ni += 1;
		} else if ( temp < -0.5 ) {
			ni += -1;
		}
		
		ti = c1 * tt + c2 * (tp-ni);
		tir= ti + ni;
		if (tir < - hn) {			// ensure the nearest point is on the helix
			tir = - hn;
		} else if (tir > hn ) {
			tir = hn;
		}

		phi= twopi * tir + phi0;
		dx = x2 - r * cos(phi);
		dy = y2 - r * sin(phi);
		dz = z2 - b * tir;
		dr2= dx * dx + dy * dy + dz * dz;
		r2 = th*th;
		
		if (dr2 > r2) {
			return -1;
		} else if (dr2 < r2) {
			return 1;
		} else {
			return 0;
			printf("warning: lattice point on particle surface\n");
		}
		
	}
	
	return -9;
}

void find_normal(real *nx, real *ny, real *nz, real x, real y, real z, int ipar, char surf)
{
	real dx, dy, dz, dr2, r2, ir;
	real phi, th, x1, y1, z1;
	real psi, b, r, phi0, x2, y2 ,z2, tp, ti, temp, tt, ni, tir, c1, c2, hn;
	
	if (ipar<0 || ipar>=npar) {
		printf("error in subroutine find_normal: particle number wrong\n");
		exit(-1);
	}
	
	if (surf=='s') {				// sphere
//		r2 = p_rad[ipar] * p_rad[ipar];
		dx = x - p_pos[ipar][0];
		dy = y - p_pos[ipar][1];
		dz = z - p_pos[ipar][2];
		if(dx<-0.5*Lx)dx+=Lx;
		if(dx> 0.5*Lx)dx-=Lx;
		if(dy<-0.5*Ly)dy+=Ly;
		if(dy> 0.5*Ly)dy-=Ly;
		if (wall_on==0) {
			if(dz<-0.5*Lz)dz+=Lz;
			if(dz> 0.5*Lz)dz-=Lz;
		}
		dr2 = dx*dx + dy*dy +dz*dz;
		ir  = sqrt(1.0/dr2);
		*nx = dx * ir;
		*ny = dy * ir;
		*nz = dz * ir;
	}
	
	if (surf=='h') {				// helix; For more details, see notes 6/9/2014
//		th  = p_thick[ipar];
		r   = p_rad[ipar];
//		hy = 0.5*p_height[ipar]+th;	// half of helix real length
//		hx1= r - th;				// bounds of the helix real radius (hx1: inner; hx2: outer)
//		hx2= r + th;
		x1 = x - p_pos[ipar][0];	// vector relative to helix C.O.M
		y1 = y - p_pos[ipar][1];
		z1 = z - p_pos[ipar][2];
		if(x1<-0.5*Lx)x1+=Lx;
		if(x1> 0.5*Lx)x1-=Lx;
		if(y1<-0.5*Ly)y1+=Ly;
		if(y1> 0.5*Ly)y1-=Ly;
		if (wall_on==0) {
			if(z1<-0.5*Lz)z1+=Lz;
			if(z1> 0.5*Lz)z1-=Lz;
		}
		rotationback(&x2, &y2, &z2, x1, y1, z1, p_qua[ipar][0], p_qua[ipar][1], p_qua[ipar][2], p_qua[ipar][3]);
		// rotate back to helix reference

		
		b   = p_b[ipar];
		phi0= p_ph0[ipar];
		hn  = 0.5 * p_ncoil[ipar];
		c2  = p_const1[ipar];		// coefficients for finding nearest point
		c1  = 1.0 - c2;
		tp   = z2/b;
		ni   = floor(tp);			// which coil the particle is within
		psi  = atan2(y2, x2) - phi0;// relative angle
		tt   = psi / twopi;
		tt   = tt - floor(tt);
		temp = tp - ni - tt;		// to tell which coil the particle is nearer to
		if ( temp > 0.5 ) {
			ni += 1;
		} else if ( temp < -0.5 ) {
			ni += -1;
		}
		
		ti = c1 * tt + c2 * (tp-ni);
		tir= ti + ni;
		if (tir < - hn) {			// ensure the nearest point is on the helix
			tir = - hn;
		} else if (tir > hn ) {
			tir = hn;
		}
		
		phi= twopi * tir + phi0;
		dx = x2 - r * cos(phi);
		dy = y2 - r * sin(phi);
		dz = z2 - b * tir;
		
		rotation(nx,ny,nz,dx,dy,dz,p_qua[ipar][0], p_qua[ipar][1], p_qua[ipar][2], p_qua[ipar][3]);
		dr2 =  dx * dx + dy * dy + dz * dz;
		ir  = sqrt(1.0/dr2);
		*nx = (*nx) * ir;
		*ny = (*ny) * ir;
		*nz = (*nz) * ir;
		
	}
	
}


//	rotate (nx, ny, nz) to (mx, my, mz)
void rotation(real *mx, real *my, real *mz, real nx, real ny, real nz, real q0, real q1, real q2, real q3)
{
	real r2, r2i=1.0;
	real q0q0, q1q1, q2q2, q3q3;
	
	q0q0 = q0 * q0;
	q1q1 = q1 * q1;
	q2q2 = q2 * q2;
	q3q3 = q3 * q3;
	
	r2 = q0q0 + q1q1 + q2q2 + q3q3;
	r2i= 1.0/r2;
	
	*mx = ( nx * (q0q0+q1q1-q2q2-q3q3) + ny * 2.0 * (q1*q2-q0*q3  ) + nz * 2.0 * (q1*q3+q0*q2  ) ) * r2i;
	*my = ( nx * (q1*q2+q0*q3        ) + ny * (q0q0-q1q1+q2q2-q3q3) + nz * 2.0 * (q2*q3-q0*q1  ) ) * r2i;
	*mz = ( nx * 2.0 * (q1*q3-q0*q2  ) + ny * 2.0 * (q2*q3+q0*q1  ) + nz * (q0q0-q1q1-q2q2+q3q3) ) * r2i;
}

// rotate (mx, my, mz) to (nx, ny, nz)
void rotationback(real *mx, real *my, real *mz, real nx, real ny, real nz, real q0, real q1, real q2, real q3)
{
	real r2, r2i=1.0;
	real q0q0, q1q1, q2q2, q3q3;
	
	q0q0 = q0 * q0;
	q1q1 = q1 * q1;
	q2q2 = q2 * q2;
	q3q3 = q3 * q3;
	
	r2 = q0q0 + q1q1 + q2q2 + q3q3;
	r2i= 1.0/r2;
	
	*mx = ( nx * (q0q0+q1q1-q2q2-q3q3) + ny * 2.0 * (q1*q2+q0*q3  ) + nz * 2.0 * (q1*q3-q0*q2  ) ) * r2i;
	*my = ( nx * (q1*q2-q0*q3        ) + ny * (q0q0-q1q1+q2q2-q3q3) + nz * 2.0 * (q2*q3+q0*q1  ) ) * r2i;
	*mz = ( nx * 2.0 * (q1*q3+q0*q2  ) + ny * 2.0 * (q2*q3-q0*q1  ) + nz * (q0q0-q1q1-q2q2+q3q3) ) * r2i;
}


void p_iden_up(real Q0[Nx][Ny][Nz][3][3])
{
	int ilist, ipar, i, j, k, ii, jj, kk;
	int i1, i2, i3, j1, j2, j3, k1, k2, k3;
	int p1, p2, p1o, p2o;
	int ilink, newlink, nlinkold=nlink;
	int isurf, p_nsurf_new;
	int temp;
	real x3, y3, z3, nx, ny, nz;
	
	for (ilist=0; ilist<nlist; ilist++) {
		ipar = link_plist[ilist][3];
		if (ipar > -1) {
			i1  = link_plist[ilist][0];
			j1  = link_plist[ilist][1];
			k1  = link_plist[ilist][2];
			p1o = link_pointo[i1][j1][k1][ipar];
			p1=onsurface(i1, j1, k1, ipar, p_type[ipar]);
			if(p1*p1o<0) {			//	only consider points that crossed the surface
				newlink = 0;
				for (ii=1; ii<15; ii++) {
					i2 = i1 + e[ii][0];
					j2 = j1 + e[ii][1];
					k2 = k1 + e[ii][2];
					
					if ( wall_on==0 || (k2>=0 && k2<Nz) ) {
						x3 = (double)i1 + 0.5*e[ii][0];
						y3 = (double)j1 + 0.5*e[ii][1];
						z3 = (double)k1 + 0.5*e[ii][2];
						if (x3<0) x3+=Nx;
						if (y3<0) y3+=Ny;
						if (z3<0) z3+=Nz;
						i3 = x3;
						j3 = y3;
						k3 = z3;
						if(i2 <  0)  i2 += Nx;
						if(i2 >= Nx) i2 -= Nx;
						if(j2 <  0)  j2 += Ny;
						if(j2 >= Ny) j2 -= Ny;
						if(k2 <  0)  k2 += Nz;
						if(k2 >= Nz) k2 -= Nz;
						p2o = link_pointo[i2][j2][k2][ipar];
						p2=onsurface(i2, j2, k2, ipar, p_type[ipar]);
						jj = bounce(ii);
						
						if (p1o*p2o>=0) {
							if (p1*p2<0) {		//	new link
								newlink             += p1;
								link_listn[nlink][0] = 1;
								link_listn[nlink][1] = ipar;
								
								if (p1<0) {		// point #1 is outside the particle
									link_listn[nlink][2] = i1;
									link_listn[nlink][3] = j1;
									link_listn[nlink][4] = k1;
									link_listn[nlink][5] = ii;
									link_listn[nlink][6] = i2;
									link_listn[nlink][7] = j2;
									link_listn[nlink][8] = k2;
									link_listn[nlink][9] = jj;
									link_point2[i1][j1][k1][ii] = -(ipar+1);
									link_point2[i2][j2][k2][jj] =  ipar+1;
									// for finite anchoring
									if (jj<=6) {
										p_surf[p_nsurf][0]   = i2;
										p_surf[p_nsurf][1]   = j2;
										p_surf[p_nsurf][2]   = k2;
										p_surf[p_nsurf][3]   = jj;
										if (p_anchor>=0) {
											p_surfQ[p_nsurf][0]  = 0.5 * (Q0[i1][j1][k1][0][0] + Q0[i2][j2][k2][0][0]);
											p_surfQ[p_nsurf][1]  = 0.5 * (Q0[i1][j1][k1][0][1] + Q0[i2][j2][k2][0][1]);
											p_surfQ[p_nsurf][2]  = 0.5 * (Q0[i1][j1][k1][0][2] + Q0[i2][j2][k2][0][2]);
											p_surfQ[p_nsurf][3]  = 0.5 * (Q0[i1][j1][k1][1][1] + Q0[i2][j2][k2][1][1]);
											p_surfQ[p_nsurf][4]  = 0.5 * (Q0[i1][j1][k1][1][2] + Q0[i2][j2][k2][1][2]);
											p_surfW[p_nsurf]     = p_anchor;
										}
										
										kk = axis_proj(ii);
										if (debug_on==1 && p_surfi[i3][j3][k3][kk]!=-1) {
											printf("error in p_surfi of 'new link'\n");
										}
										p_surfi[i3][j3][k3][kk]=p_nsurf;
										p_nsurf++;
									}
								}
								else {			// point #2 is outside the particle
									link_listn[nlink][2] = i2;
									link_listn[nlink][3] = j2;
									link_listn[nlink][4] = k2;
									link_listn[nlink][5] = jj;
									link_listn[nlink][6] = i1;
									link_listn[nlink][7] = j1;
									link_listn[nlink][8] = k1;
									link_listn[nlink][9] = ii;
									link_point2[i1][j1][k1][ii] = ipar+1;
									link_point2[i2][j2][k2][jj] = -(ipar+1);
									// for finite anchoring
									if (ii<=6) {
										p_surf[p_nsurf][0]   = i1;
										p_surf[p_nsurf][1]   = j1;
										p_surf[p_nsurf][2]   = k1;
										p_surf[p_nsurf][3]   = ii;
										// to assign the Q tensor to the new surface point, we do average over Q's
										// ( second order accuracy)
										if (p_anchor>=0) {
											p_surfQ[p_nsurf][0]  = 0.5 * (Q0[i1][j1][k1][0][0] + Q0[i2][j2][k2][0][0]);
											p_surfQ[p_nsurf][1]  = 0.5 * (Q0[i1][j1][k1][0][1] + Q0[i2][j2][k2][0][1]);
											p_surfQ[p_nsurf][2]  = 0.5 * (Q0[i1][j1][k1][0][2] + Q0[i2][j2][k2][0][2]);
											p_surfQ[p_nsurf][3]  = 0.5 * (Q0[i1][j1][k1][1][1] + Q0[i2][j2][k2][1][1]);
											p_surfQ[p_nsurf][4]  = 0.5 * (Q0[i1][j1][k1][1][2] + Q0[i2][j2][k2][1][2]);
											p_surfW[p_nsurf]     = p_anchor;
										}
										kk = axis_proj(jj);
										if (debug_on==1 && p_surfi[i3][j3][k3][kk]!=-1) {
											printf("error in p_surfi of 'new link'\n");										
										}
										p_surfi[i3][j3][k3][kk]=p_nsurf;
										p_nsurf++;
									}
								}
								nlink++;
								
								if (link_point[i2][j2][k2][ipar]==0) {
									link_plist[nlist][0] = i2;
									link_plist[nlist][1] = j2;
									link_plist[nlist][2] = k2;
									link_plist[nlist][3] = ipar;
									nlist++;
								}
								link_point[i2][j2][k2][ipar] += p2;
							}

						} else {					// old link
							if (p1*p2>0) {			// no longer a link
								link_point[i2][j2][k2][ipar] -= p2;
								link_point2[i1][j1][k1][ii] = 0;
								link_point2[i2][j2][k2][jj] = 0;
								// for finite anchoring
								if (ii<=6) {
									kk = axis_proj(ii);
									if (debug_on==1 && p_surfi[i3][j3][k3][kk]==-1) {
										printf("error in p_surfi of 'old link'\n");										
									}
									p_surfi[i3][j3][k3][kk]=-1;
								}

							} else {				// still a link
								newlink += p1;
								if (p1<0) {
									link_listn[nlink][0] = 1;
									link_listn[nlink][1] = ipar;
									link_listn[nlink][2] = i1;
									link_listn[nlink][3] = j1;
									link_listn[nlink][4] = k1;
									link_listn[nlink][5] = ii;
									link_listn[nlink][6] = i2;
									link_listn[nlink][7] = j2;
									link_listn[nlink][8] = k2;
									link_listn[nlink][9] = jj;
									nlink++;
									link_point2[i1][j1][k1][ii] = -(ipar+1);
									link_point2[i2][j2][k2][jj] = ipar+1;
									// for finite anchoring
									kk = p_surfi[i3][j3][k3][axis_proj(ii)];
									if (debug_on==1 && p_surf[kk][3]!=ii) {
										printf("error in p_surfi of 'still a link'\n");
									}
									p_surf[kk][3] = jj;
								}
							}
						}
					}
				}
				if(newlink!=0) {
					link_point[i1][j1][k1][ipar] = newlink;
				} else {
					link_point[i1][j1][k1][ipar] = 0;
					link_plist[ilist][3] = -1;
				}
			}
		}
	}
	
	for (ilist=0; ilist<nlist; ilist++) {
		ipar = link_plist[ilist][3];
		if (ipar > -1) {
			i1  = link_plist[ilist][0];
			j1  = link_plist[ilist][1];
			k1  = link_plist[ilist][2];
			p1  = link_point[i1][j1][k1][ipar];
			if (p1==0) link_plist[ilist][3]=-1;
		}
	}
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				for (ipar=0; ipar<npar; ipar++){
					link_pointo[i][j][k][ipar] = link_point[i][j][k][ipar];
					
					if(debug_on!=0){
						for (ii=0, temp=0; ii<15; ii++) {
							if(link_point2[i][j][k][ii]==ipar+1){
								temp += 1;
							} else if(link_point2[i][j][k][ii]==-ipar-1){
								temp -= 1;
							}
						}
						if (link_point[i][j][k][ipar] != temp) {
							printf("error in link_point: [%d %d %d]: %d, %d\n", i, j, k, link_point[i][j][k][ipar],temp);
							exit(-1);
						}
					}
					
				}
			}
		}
	}
	
	for (ilink=0; ilink<nlinkold; ilink++) {
		if (link_listn[ilink][0]==1) {
			ipar = link_listn[ilink][1];
			i1   = link_listn[ilink][2];
			j1   = link_listn[ilink][3];
			k1   = link_listn[ilink][4];
//			ii   = link_listn[ilink][5];
			i2   = link_listn[ilink][6];
			j2   = link_listn[ilink][7];
			k2   = link_listn[ilink][8];
//			jj   = link_listn[ilist][9];
			if (!(link_point[i1][j1][k1][ipar]<0 && link_point[i2][j2][k2][ipar]>0)) {
				link_listn[ilink][0]=0;
			}
		}
	}
	// for finite anchoring
	p_nsurf_new=0;
	for (isurf=0; isurf<p_nsurf; isurf++) {
		i  = p_surf[isurf][0];
		j  = p_surf[isurf][1];
		k  = p_surf[isurf][2];
		ii = p_surf[isurf][3];
		x3 = (real)i + 0.5 * e[ii][0];
		y3 = (real)j + 0.5 * e[ii][1];
		z3 = (real)k + 0.5 * e[ii][2];
		i1 = floor(x3);
		j1 = floor(y3);
		k1 = floor(z3);
		if (i1<0) i1 += Nx;
		if (j1<0) j1 += Ny;
		if (k1<0) k1 += Nz;
		
		jj = axis_proj(ii);

		if (p_surfi[i1][j1][k1][jj]>=0) {
			p_surfi[i1][j1][k1][jj] = p_nsurf_new;
			if (isurf>p_nsurf_new) {
				p_surf[p_nsurf_new][0] = i;
				p_surf[p_nsurf_new][1] = j;
				p_surf[p_nsurf_new][2] = k;
				p_surf[p_nsurf_new][3] = ii;				
				// need to find a way to identify ipar
				ipar = 0;
				find_normal(&nx, &ny, &nz, x3, y3, z3, ipar, p_type[ipar]);
				p_surf[p_nsurf_new][4]  = nx;
				p_surf[p_nsurf_new][5]  = ny;
				p_surf[p_nsurf_new][6]  = nz;
				for (jj=0; jj<5; jj++) {
					p_surfQ[p_nsurf_new][jj]  = p_surfQ[isurf][jj];
//					p_surfQ0[p_nsurf_new][jj] = p_surfQ0[isurf][jj];
				}
				p_surfQ0[p_nsurf_new][0] = S * (nx*nx-third);
				p_surfQ0[p_nsurf_new][1] = S * (nx*ny);
				p_surfQ0[p_nsurf_new][2] = S * (nx*nz);
				p_surfQ0[p_nsurf_new][3] = S * (ny*ny-third);
				p_surfQ0[p_nsurf_new][4] = S * (ny*nz);
				if (p_anchor<0) {
					p_surfQ[p_nsurf_new][0] = p_surfQ0[p_nsurf_new][0];
					p_surfQ[p_nsurf_new][1] = p_surfQ0[p_nsurf_new][1];
					p_surfQ[p_nsurf_new][2] = p_surfQ0[p_nsurf_new][2];
					p_surfQ[p_nsurf_new][3] = p_surfQ0[p_nsurf_new][3];
					p_surfQ[p_nsurf_new][4] = p_surfQ0[p_nsurf_new][4];
				}
			}
			p_nsurf_new++;
		}
	}
	p_nsurf = p_nsurf_new;
}


void p_link_up()
{
	// update list of links - remove inactive links
	int i, j, nlinknew=0, nlistnew=0, ipar;
	
	for (i=0; i<nlink; i++) {
		if (link_listn[i][0]==1) {
			if (i!=nlinknew) {
				for (j=0; j<10; j++) link_listn[nlinknew][j] = link_listn[i][j];
			}
			nlinknew++;
		}
	}

	nlink = nlinknew;
	
	for (i=0; i<nlist; i++) {
		ipar = link_plist[i][3];
		if (ipar>-1) {
			if (i!=nlistnew) {
				for (j=0; j<4; j++) link_plist[nlistnew][j] = link_plist[i][j];
			}
			nlistnew++;
		}
		
	}
	nlist = nlistnew;

}

void p_bc(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15])
{
	//	perform bounce-back bc on all active links
	int ilink, ipar, i1, j1, k1, i2, j2, k2, ii, jj;
	real dx, dy, dz, omegax, omegay, omegaz, ubx, uby, ubz, ub_dot_e, rho_b;
	real B, f1o, f2o, f_exch, m_exch, ex, ey, ez;
	int ilist, iig, jjg;
	real dr2i, rrad;
	
	//	cal_rho(Rho,f1);
	
	for (ilink=0; ilink<nlink; ilink++) {
		if (link_listn[ilink][0]==1) {
			ipar = link_listn[ilink][1]; 
			i1   = link_listn[ilink][2];
			j1   = link_listn[ilink][3];
			k1   = link_listn[ilink][4];
			ii   = link_listn[ilink][5];
			i2   = link_listn[ilink][6];
			j2   = link_listn[ilink][7];
			k2   = link_listn[ilink][8];
			jj   = link_listn[ilink][9];
			if(ii<7) {
				B=third;
			} else {
				B=one24th;
			}
			ex = e[ii][0];
			ey = e[ii][1];
			ez = e[ii][2];
			//	need more elaborate details
			dx = i1 - p_pos[ipar][0] + 0.5 * ex;
			dy = j1 - p_pos[ipar][1] + 0.5 * ey;
			dz = k1 - p_pos[ipar][2] + 0.5 * ez;
			if(dx<-0.5*Lx)dx+=Lx;
			if(dx> 0.5*Lx)dx-=Lx;
			if(dy<-0.5*Ly)dy+=Ly;
			if(dy> 0.5*Ly)dy-=Ly;
			omegax = p_angv[ipar][0];
			omegay = p_angv[ipar][1];
			omegaz = p_angv[ipar][2];
			ubx    = p_vel[ipar][0] + omegay*dz-omegaz*dy;
			uby    = p_vel[ipar][1] + omegaz*dx-omegax*dz;
			ubz    = p_vel[ipar][2] + omegax*dy-omegay*dx;
			ub_dot_e = ubx*ex + uby*ey + ubz*ez;
//			rho_b  = 0.5*(Rho[i1][j1][k1]+Rho[i2][j2][k2]);
//			f_exch = 2.0*B*rho_b*ub_dot_e;
			f_exch = 2.0*B*rho*ub_dot_e;
			f1o    = f1[i1][j1][k1][ii];
			f2o    = f1[i2][j2][k2][jj];
			f2[i1][j1][k1][jj] = f1o - f_exch;
			f2[i2][j2][k2][ii] = f2o + f_exch;
			
			m_exch          = f1o -f2o -2.0*f_exch;
			p_acc[ipar][0] += 2.0 * m_exch * ex;
			p_acc[ipar][1] += 2.0 * m_exch * ey;
			p_acc[ipar][2] += 2.0 * m_exch * ez;
			p_torq[ipar][0]+= 2.0 * m_exch * (dy*ez - dz*ey);
			p_torq[ipar][1]+= 2.0 * m_exch * (dz*ex - dx*ez);
			p_torq[ipar][2]+= 2.0 * m_exch * (dx*ey - dy*ex);
		}
	}
	
	/*
	 ipar=0;
	 printf("force: %f, %f, %f\n",p_acc[ipar][0],p_acc[ipar][1],p_acc[ipar][2]);
	 printf("torque: %f, %f, %f\n",p_torq[ipar][0],p_torq[ipar][1],p_torq[ipar][2]);
	 */
}

void p_update()
{
	//	update particle position, velocity, quaternion, etc.
	int ipar;
	real x, y;
	
	p_move = 0;
	for (ipar=0; ipar<npar; ipar++) {
		x = p_pos[ipar][0] + 0.;
		if(x>Nx)x-=Nx;
		y = p_pos[ipar][1] + 0.0;
		if(y>Ny)y-=Ny;
		
		if (p_pos[ipar][0]-x>1e-9 || p_pos[ipar][0]-x<-1e-9 || p_pos[ipar][1]-y>1e-9 || p_pos[ipar][1]-y<-1e-9) {
			// if the particle position has changed, set p_move=1, otherwise p_move=0
			p_move = 1;
		}
		
		p_pos[ipar][0] = x;
		p_pos[ipar][1] = y;		
		
		p_vel[ipar][1] = 0.0;
		
		// reset force and torque terms
		p_acc[ipar][0]  = 0;
		p_acc[ipar][1]  = 0;
		p_acc[ipar][2]  = 0;
		p_torq[ipar][0] = 0;
		p_torq[ipar][1] = 0;
		p_torq[ipar][2] = 0;
	}

}


int axis_proj(int ii)
{
	if (ii==1 || ii==3) {			// x direction
		return 0;
	} else if (ii==2 || ii==4) {	// y direction
		return 1;
	} else if (ii==5 || ii==6) {	// z direction
		return 2;
	}
	printf("error in axis_proj!\n");// error
	return -1;
}
