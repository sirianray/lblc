/*
 *  memory.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */
#include "main.h"

void allocate()
{
	int mpoint=points;
	
	if (wall_on!=0) {
		mpoint=qpoints;
	}
	if (npar>0) {
		if (Nx>Ny) {
			mpoint+=Nx*Nz*7;
		} else {
			mpoint+=Ny*Nz*7;
		}
	}

	Q   = malloc(mpoint*5*sizeof(real));
	info= malloc(mpoint*sizeof(int));
	Rho = malloc(points  *sizeof(real));
	u   = malloc(points*3*sizeof(real));
	
	if (Q_on!=0) {
		neighbor= malloc(mpoint*sizeof(int)*6);
		H       = malloc(mpoint*sizeof(real)*5);		
		if (wall_on!=0 && npar==0) {
			surf = malloc(nsurf*10*sizeof(real));
		} else if (npar>0) {
			surf = malloc((mpoint-points)*10*sizeof(real));
		}

		sigma_q=malloc(points*9*sizeof(real));
		sigma_p=malloc(points*3*sizeof(real));
	}
	
	if (flow_on!=0) {		
		W    = malloc(points*9 *sizeof(real));
		f    = malloc(points*15*sizeof(real));
		f2   = malloc(points*15*sizeof(real));
		Cf   = malloc(points*15*sizeof(real));
		p    = malloc(points*15*sizeof(real));
//		sigma= malloc(points*9 *sizeof(real));
		
		nextf= malloc(points*15*sizeof(int));
		if (wall_on!=0 && npar==0) {
			vsurf = malloc(nsurf*3*sizeof(real));
		} else if (npar>0) {
			vsurf = malloc((mpoint-points)*3*sizeof(real));
		}
	}
}


void deallocate()
{	
	free(Q);
	free(Rho);
	free(u);
	
	if (Q_on!=0) {		
		free(H);
		free(neighbor);
		free(info);
		if(wall_on!=0 || npar>0) free(surf);
		free(sigma_p);
		free(sigma_q);
	}
	
	if (flow_on!=0) {
		
		free(W);
		free(f);
		free(f2);
		free(Cf);
		free(p);
//		free(sigma);
		free(nextf);
		free(vsurf);
	}
	
}