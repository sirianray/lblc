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
	int i, j, id, ii, id0, iv;
	double x, y, z, q[6];
	
	if (wall_on!=0) {
		x  = n_bot[0];
		y  = n_bot[1];
		z  = n_bot[2];
		for (j=0; j<Ny; j++) {
			for (i=0; i<Nx; i++) {
				id = 10 * (i + j*Nx);
				if (Q_on!=0) {
					surf[id]   = type_bot;
					surf[id+1] = W_bot;
					surf[id+2] = 0;
					surf[id+3] = 0;
					surf[id+4] = 1.0;
					ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
					for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
				}
				
				if (flow_on!=0) {
					iv = id/10*3;
					vsurf[iv]  = 0;
					vsurf[iv+1]= uy_bottom;
					vsurf[iv+2]= 0;
				}				
			}
		}
		id0 = id + 10;
		
		x  = n_top[0];
		y  = n_top[1];
		z  = n_top[2];
		for (j=0; j<Ny; j++) {
			for (i=0; i<Nx; i++) {
				id = id0 + 10 * (i + j*Nx);
				if (Q_on!=0) {
					surf[id]   = type_top;
					surf[id+1] = W_top;
					surf[id+2] = 0;
					surf[id+3] = 0;
					surf[id+4] =-1.0;
					ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
					for (ii=0; ii<5; ii++) surf[id+5+ii]=q[ii];
				}
				
				if (flow_on!=0) {
					iv = id/10*3;
					vsurf[iv]  = 0;
					vsurf[iv+1]= uy_top;
					vsurf[iv+2]= 0;
				}				
			}
		}
	}
}


void build_neighbor()
{
	int i, j, k, ii, id0, id[6], im, ip, jm, jp, km, kmm, kp, kpp, isurf=0;
	
//	bottom wall
	if (wall_on!=0) {
		k  = 0;
		kp = 1;
		kpp= 2;
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
				
				id0  = i +  j          *Nx;
				id[0]= im+ (j + 0 *Ny) *Nx;
				id[1]= ip+ (j + 0 *Ny) *Nx;
				id[2]= i + (jm+ 0 *Ny) *Nx;
				id[3]= i + (jp+ 0 *Ny) *Nx;
				id[4]= i + (j + kp *Ny)*Nx;
				id[5]= i + (j + kpp*Ny)*Nx;
				
				for (ii=0; ii<6; ii++) {
					neighbor[6*id0+ii] = 5*id[ii];
				}				
				info[id0] = isurf*10;
				isurf++;
			}
		}
	}

//	bulk point
	for (k=0; k<Nz; k++) {
		km = k - 1;
		kp = k + 1;
		if (wall_on==0) {
			if (km<0) km+=Nz;
			if (kp>=Nz) kp-=Nz;
		}		
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
				
				id0  = (i + (j + k *Ny)*Nx) + bulk0;
				id[0]= (im+ (j + k *Ny)*Nx) + bulk0;
				id[1]= (ip+ (j + k *Ny)*Nx) + bulk0;
				id[2]= (i + (jm+ k *Ny)*Nx) + bulk0;
				id[3]= (i + (jp+ k *Ny)*Nx) + bulk0;
				id[4]= (i + (j + km*Ny)*Nx) + bulk0;
				id[5]= (i + (j + kp*Ny)*Nx) + bulk0;
				
				info[id0]=-1;
				for (ii=0; ii<6; ii++) {
					neighbor[6*id0+ii]= 5*id[ii];					
				}
			}
		}
	}
	
//	top wall
	if (wall_on!=0) {
		k  = Nz;
		km = Nz-1;
		kmm= Nz-2;
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
				
				id0  = (i + (j + k  *Ny)*Nx) + bulk0;
				id[0]= (im+ (j + k  *Ny)*Nx) + bulk0;
				id[1]= (ip+ (j + k  *Ny)*Nx) + bulk0;
				id[2]= (i + (jm+ k  *Ny)*Nx) + bulk0;
				id[3]= (i + (jp+ k  *Ny)*Nx) + bulk0;
				id[4]= (i + (j + km *Ny)*Nx) + bulk0;
				id[5]= (i + (j + kmm*Ny)*Nx) + bulk0;
				
				for (ii=0; ii<6; ii++) {
					neighbor[6*id0+ii] = 5*id[ii];
				}				
				info[id0] = isurf*10;
				isurf++;
			}
		}
	}
	
	//	before assigning particle surface point, change sign of wall neighbor point
	if (wall_on!=0) {
		for (i=0; i<qpoints*6; i++) {
			if (info[neighbor[i]/5]>-1) {
				neighbor[i] =-neighbor[i]-1;
			}
		}
	}
}

void build_stream()
{
	int i, j, k, in, jn, kn, ii, iip, id, idn;
	
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny;j++){
			for(i=0;i<Nx;i++){
				id = 15 * (i + (j+k*Ny)*Nx);
				for(ii=0;ii<15;ii++) {
					in = i + e[ii][0];
					jn = j + e[ii][1];
					kn = k + e[ii][2];
					if(in==-1)in=Nx-1;
					if(in==Nx)in=0;
					if(jn==-1)jn=Ny-1;
					if(jn==Ny)jn=0;
					
					iip = bounce(ii);
					
					if (kn>=0 && kn<Nz) {
						idn = 15 * (in + (jn+kn*Ny)*Nx);
						nextf[id+ii] = idn+ii;
					}
					else {
						if (wall_on==0) {
							if(kn==-1)kn=Nz-1;
							if(kn==Nz)kn=0;
							idn = 15 * (in + (jn+kn*Ny)*Nx);
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
