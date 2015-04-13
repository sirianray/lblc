#include "main.h"
#include <time.h>

void init()
{
	int i, j, k, ii, ip, id;
	time_t t;
	double nx, ny, nz, q[6], sita0, sita, lambda1, lambda2;
	
	if (newrun_on!=0) {
		sita0 = n_top[0]*n_bot[0]+n_top[1]*n_bot[1]+n_top[2]*n_bot[2];
		sita0 = sita0/sqrt(n_top[0]*n_top[0]+n_top[1]*n_top[1]+n_top[2]*n_top[2]);        
		sita0 = sita0/sqrt(n_bot[0]*n_bot[0]+n_bot[1]*n_bot[1]+n_bot[2]*n_bot[2]);
		sita0 = acos(sita0);
		
		if (rand_init==0) {				// random random seed
			srand((unsigned)time(&t));
		} else if (rand_init>0) {		// random seed from input
			srand((unsigned)rand_seed);
		}
		
		if (wall_on!=0 && Q_on!=0) {
			for (i=0; i<Nx*Ny; i++) {
				ip = info[i];
				for (j=0; j<5; j++) Q[i*5+j] = surf[ip+5+j];
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
					id = (i + (j+k*Ny)*Nx) + bulk0;
					if (rand_init>=0) {
						nx = randvec();
						ny = randvec();
						nz = randvec();
					}
					ntoq(nx, ny, nz, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
					for (ii=0; ii<5; ii++) Q[id*5+ii] = q[ii];
					if (q[0]+q[3]+q[5]>1e-15 || q[0]+q[3]+q[5]<-1e-15 ) {
						printf("initial bulk Q not right\n");
						exit(-1);
					}
				}
			}
		}

		if (wall_on!=0 && Q_on!=0) {
			for (i=Nx*Ny*(Nz+1); i<Nx*Ny*(Nz+2); i++) {
				ip = info[i];
				for (j=0; j<5; j++) Q[i*5+j] = surf[ip+5+j];
			}
		}
		
		if (flow_on!=0) {
			for (k=0; k<Nz; k++) {
				for (j=0; j<Ny; j++) {
					for (i=0; i<Nx; i++) {
						id = i + (j+k*Ny)*Nx;
						Rho[id]  = rho;
						u[id*3]  = 0;
						u[id*3+1]= uy_bottom + (uy_top-uy_bottom)*((double)k+0.5)/(double)Nz;
						u[id*3+2]= 0;
					}
				}
			}			
		}
		
	} else {
		read_restart();
	}
}


void ntoq(double nx, double ny, double nz, double *q0, double *q1, double *q2, double *q3, double *q4, double *q5)
{
	double r2, r2i;
	
	r2 = nx*nx + ny*ny + nz*nz;
	if (debug_on!=0 && r2<1e-8) {
		printf("error: zero director encountered in function ntoq\n");
		exit(-1);
	}
	r2i= 1.0/r2;	
	*q0 = S_lc * (nx*nx*r2i - third);
	*q1 = S_lc * (nx*ny*r2i);
	*q2 = S_lc * (nx*nz*r2i);
	*q3 = S_lc * (ny*ny*r2i - third);
	*q4 = S_lc * (ny*nz*r2i);
	*q5 = -(*q0) - (*q3);
}



real randvec()
{
	double r;
	
    r =  1.0 - 2.0 * (double)rand()/(double)RAND_MAX ;
    return r;
}