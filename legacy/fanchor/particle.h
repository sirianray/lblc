/*
 *  particle.h
 *  
 *
 *  Created by Sirius on 4/4/14;
 *  Modified by Sirius on 6/8/2014;
 *  Copyright 2014 home. All rights reserved.
 *
 */

/*
 LINK VARAIBLES:
 nlink:   number of bb links
 nlist:   number of registered points (meaning it is/was on a bb link)
 link_point:  Nx*Ny*Nz*NPMAX matrix that stores # of links of lattice point [i,j,k] to ipar'th particle
 link_pointo: link_point at previous time step
 link_point2: Nx*Ny*Nz*15 matrix that stores particle number that the link has crossed:
			-(ipar+1): inside ipar's; ipar+1: outside ipar's; 0: not on surface
 link_listn:  integer info of the links - [effectiveness, ipar, i1, j1, k1, ii, i2, j2, k2, jj]; 1; outside; 2: inside
 link_listf:  float number info of the links - including the position, velocity & interception of the surface point
 associated to the link
 link_plist:  info of registered points - [i, j, k, ipar]
 
 PARTICLE RELATED VARIABLES:
 npar:                     number of particles
 p_type:                   particle shape type: 's' - spherical; 'h' - helical
 p_m, p_I, p_rad & p_size: masses, moments of inertia, radii and sizes of particles
 p_pos, p_vel, p_acc:      (x/y/z) coordinates, velocities & accelerations of particles
 p_qua:                    quaternions of particles
 p_orien & p_orien2:       biaxil orientations of the particle
 p_angv & p_torq:          angular velocities & torques of particles
 p_move:					detect if particle has moved (if not, some subroutines can be skipped)
 
 HELIX SPECIAL VARIABLES:
 p_ncoil:  number of helical coils (can be uninteger)
 p_b:      height of each coil
 p_thick:  radius of the helical thread (cylindrical like)
 p_ph0:    phase of the trajectory equation for helix
 p_height: helix length
 p_const1: constants for calculating the nearest point on the thread to a given point = b^2/((2*pi*r)^2+b^2)
 
 IMPORTANT NOTICE:
 - lattice points may appear multiple times in link_plist as they may simultaneously be the surface points to multiple
	particles; link_plist[i][3] must be within [0,npar-1]; -1 means the point is no longer a surface point to ipar'th
	particle; points with -1 will be removed in p_list_up
 - link_point[i][j][k][ipar] is the number of inner/outer links of lattice point [i,j,k] to ipar'th particle;
	+: inside; -: outside; 0: it is not the surface point to particle ipar but it may still be inside the particle
 - the definition of p_rad depends on the particle type. If the particle is spherical, p_rad is the radius; If the
	particle is helix, p_rad is the radius of the coil. p_size = p_rad for sphere; p_size = half diagonal length for helix
 - parametric equation for the helix thread trajectory: x = r*cos(twopi * t + phi0), x = r*sin(twopi * t + phi0) , z = bt.
 */

#include "lb.h"
#define NPMAX 10
#define NLINKMAX  22*Nx*Nx
#define halfsqrt2 0.7071067811865475
#define twopi     6.2831853071795862

extern int npar, nlink, nlist, p_move;
extern int link_point[Nx][Ny][Nz][NPMAX], link_listn[NLINKMAX][10], link_plist[NLINKMAX][4];
extern int link_pointo[Nx][Ny][Nz][NPMAX];
extern int link_point2[Nx][Ny][Nz][15];
extern char p_type[NPMAX];
extern real Lx, Ly, Lz;
//extern real **link_listf;
extern real p_m[NPMAX], p_I[NPMAX][2], p_rad[NPMAX], p_size[NPMAX], p_orien[NPMAX][3], p_orien2[NPMAX][3], p_norm[NPMAX][3];
extern real p_ncoil[NPMAX], p_b[NPMAX], p_thick[NPMAX], p_ph0[NPMAX], p_height[NPMAX], p_const1[NPMAX];
extern real p_pos[NPMAX][3], p_vel[NPMAX][3], p_acc[NPMAX][3], p_qua[NPMAX][4], p_angv[NPMAX][3], p_torq[NPMAX][3];

// for finite anchoring
// p_nsurf:  number of surface points
// p_surfi:  label of the surface points; default value=-1 (not a surface point)
// p_surf:   info of surface points: inner end coordinate, direction & normal;
//				default direction=-1 (not an active point); available directio=1-6 (active point)
// p_surfQ & p_surfQ0: Q tensor associated with the surface point
// p_surfdQ: gradQ dot n - derivative of Q tensor along the surface normal
// p_anchor & p_w: anchoring strength
extern int p_nsurf;
extern int p_surfi[Nx][Ny][Nz][3];
extern real p_surf[NLINKMAX][7], p_surfQ0[NLINKMAX][5], p_surfQ[NLINKMAX][5], p_surfdQ[NLINKMAX][5], p_surfW[NLINKMAX];
extern real p_anchor;

void p_init();		// allocation & initialization of p_* variables
void p_iden();		// identify all links and associated points

void p_iden_up(real Q0[Nx][Ny][Nz][3][3]);						// update the identifications of links & points
void p_link_up();												// update link/point lists
void p_bc(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15]);	// apply bounce-back bc to the links
void p_update();												// update particle positions, etc.

// geometry related functions
int onsurface(real x, real y, real z, int ipar, char surf);		// -1: outside; 0: on; 1: inside; -9: error
void find_normal(real *nx, real *ny, real *nz, real x, real y, real z, int ipar, char surf);
// rotate (nx, ny, nz) to (mx, my, mz)
void rotation(real *mx, real *my, real *mz, real nx, real ny, real nz, real q0, real q1, real q2, real q3);
// rotate (mx, my, mz) to (nx, ny, nz)
void rotationback(real *mx, real *my, real *mz, real nx, real ny, real nz, real q0, real q1, real q2, real q3);
// axis projection: from 1:6 (six directions) to 0:2 (3 dimensions)
int axis_proj(int ii);