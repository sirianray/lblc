#include "main.h"
#include "particle.h"
#include <time.h>

void init()
{
	int i, j, k, ii, ip, id, ipar;
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
			
			if (rand_init==0) {			// random random seed
				srand((unsigned)time(&t));
			} else if (rand_init==1) {		// random seed from input
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
						if (rand_init==0 || rand_init==1) {
							nx = randvec();
							ny = randvec();
							nz = randvec();
						}
						else if (rand_init==-2 && npar>0) {
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
						else if (rand_init==-3 || rand_init==-4) {
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
						} else if (rand_init==-32 || rand_init==-42) {
							double A_init = q_init;
							double cst;
							double isq2 = 1.0 / sqrt(2);
							double sq2 = sqrt(2);
							double xi, yi, zi;
							double cosquatpi=cos(0.25*PI);
							if (rand_init == -32){
								cst = 2 * q_ch * 0.68;
							}
							else{
								cst = 2 * q_ch * 0.86;
							}
							if (rand_init == -32){
								xi = (i-Nx*0.5+0.5) * cst * isq2;
								yi = (j-Ny*0.5+0.5) * cst * isq2;
								zi = (k-Nz*0.5+0.5) * cst * isq2;
								x  = xi;
								y  = cosquatpi*yi + cosquatpi*zi;
								z  =-cosquatpi*yi + cosquatpi*zi;
								q[0] = A_init * (- sin(y) * cos(x) - sin(x) * cos(z) + 2 * sin(z) * cos(y));
								q[3] = A_init * (- sin(z) * cos(y) - sin(y) * cos(x) + 2 * sin(x) * cos(z));
								q[5] = A_init * (- sin(x) * cos(z) - sin(z) * cos(y) + 2 * sin(y) * cos(x));
								q[1] = A_init * (- sq2 * sin(x) * sin(z) - sq2 * cos(y) * cos(z) + sin(x) * cos(y));
								q[2] = A_init * (- sq2 * sin(z) * sin(y) - sq2 * cos(x) * cos(y) + sin(z) * cos(x));
								q[4] = A_init * (- sq2 * sin(y) * sin(x) - sq2 * cos(z) * cos(x) + sin(y) * cos(z));
							}
							else if (rand_init ==-42){
								xi   = i - Nx * 0.5 + 0.5;
								yi   = j - Ny * 0.5 + 0.5;
								zi   = k - Nz * 0.5 + 0.5;
								x    = xi;
								y    = cosquatpi*yi + cosquatpi*zi;
								z    =-cosquatpi*yi + cosquatpi*zi;
								q[0] = A_init * (cos(cst * z) - cos(cst * y));
								q[1] = A_init * sin(cst * z);
								q[2] = A_init * sin(cst * y);
								q[3] = A_init * (cos(cst * x) - cos(cst * z));
								q[4] = A_init * sin(cst * x);
								q[5] = A_init * (cos(cst * y) - cos(cst * x));
							}
						} else if (rand_init==-5 || rand_init==-6) {
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
						}

						if (rand_init!=-3 && rand_init!=-4 && rand_init!=-32 && rand_init!=-42) {
							ntoq(nx, ny, nz, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
						}
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
