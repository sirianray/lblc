/*
 *  geometry.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */

#include "main.h"

void init_surf()
{
	int i, j, k, id, ii, id0, iv;
	double x, y, z, q[6];
	
	id = -10;
	if (wall_x!=0 && Q_on!=0 && myid==0) {
		x  = n_xlo[0];
		y  = n_xlo[1];
		z  = n_xlo[2];
		for (k=0; k<Nz; k++) {
			for (j=0; j<Ny; j++) {
				id        += 10;
				surf[id]   = type_xlo;
				surf[id+1] = W_xlo;
				surf[id+2] = 1.0;
				surf[id+3] = 0;
				surf[id+4] = 0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
		
		x  = n_xhi[0];
		y  = n_xhi[1];
		z  = n_xhi[2];
		for (k=0; k<Nz; k++) {
			for (j=0; j<Ny; j++) {
				id        += 10;
				surf[id]   = type_xhi;
				surf[id+1] = W_xhi;
				surf[id+2] =-1.0;
				surf[id+3] = 0;
				surf[id+4] = 0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
	}
	
	if (wall_y!=0 && Q_on!=0 && myid==0) {
		x  = n_ylo[0];
		y  = n_ylo[1];
		z  = n_ylo[2];
		for (k=0; k<Nz; k++) {
			for (i=0; i<Nx; i++) {
				id        += 10;
				surf[id]   = type_ylo;
				surf[id+1] = W_ylo;
				surf[id+2] = 0;
				surf[id+3] = 1.0;
				surf[id+4] = 0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
		
		x  = n_yhi[0];
		y  = n_yhi[1];
		z  = n_yhi[2];
		for (k=0; k<Nz; k++) {
			for (i=0; i<Nx; i++) {
				id        += 10;
				surf[id]   = type_yhi;
				surf[id+1] = W_yhi;
				surf[id+2] = 0;
				surf[id+3] =-1.0;
				surf[id+4] = 0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
	}
	
	if (wall_z!=0 && Q_on!=0 && myid==0) {
		x  = n_bot[0];
		y  = n_bot[1];
		z  = n_bot[2];
		for (j=0; j<Ny; j++) {
			for (i=0; i<Nx; i++) {
				id        += 10;
				surf[id]   = type_bot;
				surf[id+1] = W_bot;
				surf[id+2] = 0;
				surf[id+3] = 0;
				surf[id+4] = 1.0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
		
		x  = n_top[0];
		y  = n_top[1];
		z  = n_top[2];
		for (j=0; j<Ny; j++) {
			for (i=0; i<Nx; i++) {
				id        += 10;
				surf[id]   = type_top;
				surf[id+1] = W_top;
				surf[id+2] = 0;
				surf[id+3] = 0;
				surf[id+4] =-1.0;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
			}
		}
	}
	
	if ((wall_x!=0 || wall_y!=0 || wall_z!=0) && Q_on!=0) MPI_Win_fence(0, winsurf);
	MPI_Barrier(MPI_COMM_WORLD);
}


void build_neighbor()
{
	int i, j, k, ii, id0, id[6], im, ip, jm, jp, km, kmm, kp, kpp, isurf=0;
	int nwally=0, nwallz=0;

	if (myid == root) printf("Building mesh neighbors..\n");
	
	if (wall_x!=0) {
		nwally += 2*Ny*Nz;
		nwallz += 2*Ny*Nz;
	}
	if (wall_y!=0) {
		nwallz += 2*Nx*Nz;
	}
	
	for (k=0; k<Nz; k++) {
		km = k - 1;
		kp = k + 1;
		if (wall_z==0) {
			if (km<0) km+=Nz;
			if (kp>=Nz) kp-=Nz;
		}
		for (j=0; j<Ny; j++) {
			jm = j - 1;
			jp = j + 1;
			if (wall_y==0) {
				if (jm<0) jm+=Ny;
				if (jp>=Ny) jp-=Ny;
			}
			
			for (i=0; i<Nx; i++) {
				im = i - 1;
				ip = i + 1;
				if (wall_x==0) {
					if (im<0) im+=Nx;
					if (ip>=Nx) ip-=Nx;
				}				
				
				id0  = i + (j + k *Ny)*Nx;

				if (id0/point==myid) {
					id0  = id0%point;
					if (im>=0) {
						id[0]= (im+ (j + k *Ny)*Nx - myid * point) * 5;
					} else {
						id[0]= ( j + k *Ny  - myid * node) * 5 - 1;
					}
					if (ip<Nx) {
						id[1]= (ip+ (j + k *Ny)*Nx - myid * point) * 5;
					} else {
						id[1]= ( j+(k+Nz)*Ny- myid * node) * 5 - 1;
					}
					if (jm>=0) {
						id[2]= (i + (jm + k *Ny)*Nx - myid * point) * 5;
					} else {
						id[2]= ( nwally + i + k *Nx - myid * node) * 5 - 1;
					}
					if (jp<Ny) {
						id[3]= (i + (jp + k *Ny)*Nx - myid * point) * 5;
					} else {
						id[3]= ( nwally + i+(k+Nz)*Nx- myid * node) * 5 - 1;
					}
					if (km>=0) {
						id[4] = (i + (j + km*Ny)*Nx - myid * point) * 5;
					} else {
						id[4] = (nwallz + i +  j*Nx - myid * node ) * 5 - 1;
					}
					if (kp<Nz) {
						id[5] = (i + (j + kp*Ny)*Nx - myid * point) * 5;
					} else {
						id[5] = (nwallz+i+(j+Ny)*Nx - myid * node ) * 5 - 1;
					}
					info[id0]=-1;
					for (ii=0; ii<6; ii++) {
						neighb[6*id0+ii] = id[ii];
					}
				}
			}
		}
	}
	
	if (wall_x!=0) {
		for (k=0; k<Nz; k++) {
			km = k - 1;
			kp = k + 1;
			if (km<0)  km+=Nz;
			if (kp>=Nz)kp-=Nz;
			for (j=0; j<Ny; j++) {
				jm = j - 1;
				jp = j + 1;
				if (jm<0)  jm+=Ny;
				if (jp>=Ny)jp-=Ny;
				// x low wall
				id0  =  j + k* Ny;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (    (j+k*Ny)*Nx - myid * point) * 5;
					id[1]= (1 + (j+k*Ny)*Nx - myid * point) * 5;
					id[2]= (jm + k *Ny - myid * node) * 5;
					id[3]= (jp + k *Ny - myid * node) * 5;
					id[4]= (j  + km*Ny - myid * node) * 5;
					id[5]= (j  + kp*Ny - myid * node) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
				// x high wall
				id0  =  j + (k +Nz)*Ny;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (Nx-1+(j+k*Ny)*Nx - myid * point) * 5;
					id[1]= (Nx-2+(j+k*Ny)*Nx - myid * point) * 5;
					id[2]= (jm  + k *Ny - myid * node) * 5;
					id[3]= (jp  + k *Ny - myid * node) * 5;
					id[4]= (j   + km*Ny - myid * node) * 5;
					id[5]= (j   + kp*Ny - myid * node) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
			}
		}
	}
	
	if (wall_y!=0) {
		for (k=0; k<Nz; k++) {
			km = k - 1;
			kp = k + 1;
			if (km<0)  km+=Nz;
			if (kp>=Nz)kp-=Nz;
			for (i=0; i<Nx; i++) {
				im = i - 1;
				ip = i + 1;
				if (im<0)  im+=Nx;
				if (ip>=Nx)ip-=Nx;
				// y low wall
				id0  = nwally + i + k* Nx;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (nwally + im + k *Nx - myid*node) * 5;
					id[1]= (nwally + ip + k *Nx - myid*node) * 5;
					id[2]= (i + (     k*Ny)*Nx - myid*point) * 5;
					id[3]= (i + (1 +  k*Ny)*Nx - myid*point) * 5;
					id[4]= (nwally + i +  km*Nx - myid*node) * 5;
					id[5]= (nwally + i +  kp*Nx - myid*node) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
				// y high wall
				id0  = nwally + i + (k+Nz)* Nx;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (nwally + im + k *Nx - myid*node) * 5;
					id[1]= (nwally + ip + k *Nx - myid*node) * 5;
					id[2]= (i + (Ny-1+k*Ny)*Nx - myid*point) * 5;
					id[3]= (i + (Ny-2+k*Ny)*Nx - myid*point) * 5;
					id[4]= (nwally + i +  km*Nx - myid*node) * 5;
					id[5]= (nwally + i +  kp*Nx - myid*node) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
			}
		}
	}

	if (wall_z!=0) {
                for (j=0; j<Ny; j++) {
                        jm = j - 1;
                        jp = j + 1;
                        if (jm<0) jm+=Ny;
                        if (jp>=Ny) jp-=Ny;
                        for (i=0; i<Nx; i++) {
                                im = i - 1;
                                ip = i + 1;
                                if (im<0) im+=Nx;
                                if (ip>=Nx) ip-=Nx;
			// bottom wall	
				id0  = nwallz + i + j* Nx;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (nwallz+im+ j *Nx - myid * node) * 5;
					id[1]= (nwallz+ip+ j *Nx - myid * node) * 5;
					id[2]= (nwallz+i + jm*Nx - myid * node) * 5;
					id[3]= (nwallz+i + jp*Nx - myid * node) * 5;
					id[4]= (i +     j    *Nx - myid *point) * 5;
					id[5]= (i +    (j+Ny)*Nx - myid *point) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
				// top wall
				id0  =  nwallz + i + (j +Ny)*Nx;
				if (id0/node==myid) {
					id0  = id0%node;
					id[0]= (nwallz+im+ (j +Ny)*Nx - myid * node) * 5;
					id[1]= (nwallz+ip+ (j +Ny)*Nx - myid * node) * 5;
					id[2]= (nwallz+i + (jm+Ny)*Nx - myid * node) * 5;
					id[3]= (nwallz+i + (jp+Ny)*Nx - myid * node) * 5;
					id[4]= (i + (j+(Nz-1)*Ny)*Nx - myid * point) * 5;
					id[5]= (i + (j+(Nz-2)*Ny)*Nx - myid * point) * 5;
					for (ii=0; ii<6; ii++) {
						neighbsurf[6*id0+ii] = id[ii];
					}
				}
			}
		}
	}
}

void build_stream()
{
	int i, j, k, in, jn, kn, ii, iip, id, idn, iproc;

	if (myid == root) printf("Building streaming functions..\n");
	
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				id = 15 * (i + (j+k*Ny)*Nx);
				iproc = id/(point*15);
				if (iproc==myid) {
					id = id%(point*15);
					for(ii=0;ii<15;ii++) {
						in = i + e[ii][0];
						jn = j + e[ii][1];
						kn = k + e[ii][2];
						if (wall_x==0) {
							if(in==-1)in=Nx-1;
							if(in==Nx)in=0;
						}
						if (wall_y==0) {
							if(jn==-1)jn=Ny-1;
							if(jn==Ny)jn=0;
						}
						if (wall_z==0) {
							if(kn==-1)kn=Nz-1;
							if(kn==Nz)kn=0;
						}
						
						iip = bounce(ii);
						
						if (in>=0 && in<Nx && jn>=0 && jn<Ny && kn>=0 && kn<Nz) {
							idn = 15 * (in + (jn+kn*Ny)*Nx - myid * point);
							nextf[id+ii] = idn+ii;
						}
						else {
							nextf[id+ii] = id+iip;
						}
					}
				}
			}
		}
	}
}

void lattice_vec()
{
	int i;
	
	for(i=0;i<3;i++){
		e[0][i] = 0;
	}
	
	e[1][0] =-1;
	e[1][1] = 0;
	e[1][2] = 0;
	
	e[2][0] = 1;
	e[2][1] = 0;
	e[2][2] = 0;
	
	e[3][0] = 0;
	e[3][1] =-1;
	e[3][2] = 0;
	
	e[4][0] = 0;
	e[4][1] = 1;
	e[4][2] = 0;
	
	e[5][0] =  0;
	e[5][1] =  0;
	e[5][2] = -1;
	
	e[6][0] =  0;
	e[6][1] =  0;
	e[6][2] =  1;
	
	e[7][0] = +1;
	e[7][1] = +1;
	e[7][2] = +1;
	
	e[8][0] = -1;
	e[8][1] = 1;
	e[8][2] = 1;
	
	e[9][0] = -1;
	e[9][1] = -1;
	e[9][2] = +1;
	
	e[10][0] = +1;
	e[10][1] = -1;
	e[10][2] = +1;
	
	e[11][0] = +1;
	e[11][1] = +1;
	e[11][2] = -1;
	
	e[12][0] = -1;
	e[12][1] = +1;
	e[12][2] = -1;
	
	e[13][0] = -1;
	e[13][1] = -1;
	e[13][2] = -1;
	
	e[14][0] = +1;
	e[14][1] = -1;
	e[14][2] = -1;	
}

int bounce(int i)
{	
	switch( i ) 
	{
		case 0:
			return 0;
			break;
		case 1 :
			return 2;
			break;
		case 2 :
			return 1;
			break;
		case 3 :
			return 4;
			break;
		case 4 :
			return 3;
			break;
		case 5 :
			return 6;
			break;
		case 6 :
			return 5;
			break;
		case 7 :
			return 13;
			break;
		case 8 :
			return 14;
			break;
		case 9 :
			return 11;
			break;
		case 10 :
			return 12;
			break;
		case 11 :
			return 9;
			break;
		case 12 :
			return 10;
			break;
		case 13 :
			return 7;
			break;
		case 14 :
			return 8;
			break;
		default:
			printf("error: index of bounce() is out of range!");
			return -1;
	}
}
