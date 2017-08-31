#include "main.h"

//Read in parameters from input FILE
void read_param()
{
	int ipar, i, flag;
	double x, y, z, r, Wt, W, r2, ir;
	char junk[256], ishape;
	FILE *param;

	param = fopen("param.in", "r");

	if (param==NULL) {	// make sure param file exists
		printf("input file param.in not exist\n");
		exit(-1);
	} else if (myid==root) {
		printf("Reading param.in ...\n");
	}

	flag = 1;		// skip info lines which start with #
	while ( flag==1 && fgets (junk , 256 , param) != NULL ) {
		if(junk[0]=='#') {
			flag=1;
		} else {
			flag=0;
		}
	}

	sscanf(junk,  "newrun_on %d\n", &newrun_on);
//	fscanf(param, "newrun_on %d\n", &newrun_on);
	fscanf(param, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
	fscanf(param, "npar %d\n", &npar);
	fscanf(param, "patch_on %d\n", &patch_on);
	fscanf(param, "wall_x wall_y wall_z %d %d %d\n", &wall_x, &wall_y, &wall_z);
	fscanf(param, "flow_on %d\n", &flow_on);
	fscanf(param, "Q_on %d\n", &Q_on);

	fscanf(param, "rand_init rand_seed q_init %d %d %le\n", &rand_init, &rand_seed, &q_init);
	
	fscanf(param, "t_max t_print t_write %d %d %d\n", &t_max, &t_print, &t_write);
	
	fscanf(param, "n_evol_Q %d\n", &n_evol_Q);
	
	fscanf(param, "type_xlo W_xlo n_xlo %d %lf %lf %lf %lf\n", &type_xlo, &W_xlo, &n_xlo[0], &n_xlo[1], &n_xlo[2]);
	fscanf(param, "type_xhi W_xhi n_xhi %d %lf %lf %lf %lf\n", &type_xhi, &W_xhi, &n_xhi[0], &n_xhi[1], &n_xhi[2]);
	fscanf(param, "type_ylo W_ylo n_ylo %d %lf %lf %lf %lf\n", &type_ylo, &W_ylo, &n_ylo[0], &n_ylo[1], &n_ylo[2]);
	fscanf(param, "type_yhi W_yhi n_yhi %d %lf %lf %lf %lf\n", &type_yhi, &W_yhi, &n_yhi[0], &n_yhi[1], &n_yhi[2]);
	fscanf(param, "type_bot W_bot n_bot %d %lf %lf %lf %lf\n", &type_bot, &W_bot, &n_bot[0], &n_bot[1], &n_bot[2]);
	fscanf(param, "type_top W_top n_top %d %lf %lf %lf %lf\n", &type_top, &W_top, &n_top[0], &n_top[1], &n_top[2]);
	
	fscanf(param, "tau_f %le\n", &tau_f);
	fscanf(param, "ux_lo uy_lo uz_lo %lf %lf %lf\n", &ux_lo, &uy_lo, &uz_lo);
	fscanf(param, "ux_hi uy_hi uz_hi %lf %lf %lf\n", &ux_hi, &uy_hi, &uz_hi);
	fscanf(param, "xforce yforce zforce %le %le %le\n", &xforce, &yforce, &zforce);
	
	fscanf(param, "L1 L2 L3 L4 %le %le %le %le\n", &L1, &L2, &L3, &L4);
	fscanf(param, "rho %le\n", &rho);
	fscanf(param, "A_ldg %le\n", &A_ldg);
	fscanf(param, "q_ch %le\n", &q_ch);
	fscanf(param, "U %le\n", &U_lc);
	fscanf(param, "xi %le\n", &xi);
	fscanf(param, "Gamma_rot %le\n", &Gamma_rot);
	
	fscanf(param, "debug_on %d\n",&debug_on);

	fscanf(param, "Q_tol u_tol %le %le\n", &Q_tol, &u_tol);
	
	fclose(param);

	if (npar>0) param = fopen("ppos.in","r");

	points = Nx*Ny*Nz;
	point=ceil((double)points/(double)numprocs);
	lpoint=point;
	if(myid==numprocs-1)lpoint=points-(numprocs-1)*point;

	bulk0  = Nx*Ny;
	nsurf  = 0;
	nodes  = 0;	
	if (wall_x!=0) {
		nsurf += Ny*Nz*2;
		nodes += Ny*Nz*2;
	}
	if (wall_y!=0) {
		nsurf += Nx*Nz*2;
		nodes += Nx*Nz*2;
	}
	if (wall_z!=0) {
		nsurf += Nx*Ny*2;
		nodes += Nx*Ny*2;
	}

	for (ipar=0; ipar<npar; ipar++) {
		fscanf(param, "%le %le %le %le %le %le %c\n", &x, &y, &z, &r, &Wt, &W, &ishape);
		switch (ishape) {
			case 'c':
			if (x!=0) {
				nodes += 2*(Ny+Nz)*Nx;
			} else if (y!=0) {
				nodes += 2*(Nx+Nz)*Ny;
			} else if (z!=0) {
				nodes += 2*(Nx+Ny)*Nz;
			}
			break;

			case 'q':
			if (x!=0) {
				nodes += 2*(Ny+Nz)*Nx;
			} else if (y!=0) {
				nodes += 2*(Nx+Nz)*Ny;
			} else if (z!=0) {
				nodes += 2*(Nx+Ny)*Nz;
			}
			break;

			default:
			if (r<4.0) {
				nodes += ceil(r*r*31.0);
			} else if (r<8.0) {
				nodes += ceil(r*r*22.0);
			} else if (r<14.0) {
				nodes += ceil(r*r*20.0);
			} else {
				nodes += ceil(r*r*19.5);
			}
		}
	}
	if (npar>0) fclose(param);

	node  =ceil((double)nodes/(double)numprocs);
//	qpoint=ceil((double)qpoints/(double)numprocs);

	if(myid==0) printf("nodes, node=%d %d\n",nodes,node);

	itau_f = 1.0/tau_f;
	qdt    = 1/(double)n_evol_Q;
	S_lc   = getS(U_lc);	
	Fld0   = getF0(U_lc);
	xi1    = 0.5*(xi+1.0);
	xi2    = 0.5*(xi-1.0);
	twqL_ch= 2.*q_ch*L1;
	
	if (flow_on!=0) uconverge=0;
	if (Q_on!=0)    qconverge=0;

//  calculate frank elastic modulii
    K1 = (2.*L1 + L2 - two3rd*S_lc*L3 + L4)*S_lc*S_lc;
    K2 = (   L1      -  third*S_lc*L3     )*S_lc*S_lc*2.;
    K3 = (2.*L1 + L2 +four3rd*S_lc*L3 + L4)*S_lc*S_lc;
    K4 = L4*S_lc*S_lc;
    K24= K2+K4;
    if (myid==0) printf("Elastic constants: K1=%f, K2=%f, K3=%f, K24=%f\n",K1,K2,K3,K24);

	r2 = n_top[0]*n_top[0]+n_top[1]*n_top[1]+n_top[2]*n_top[2];
	if (r2>0 && r2!=1) {
		ir = 1.0/sqrt(r2);
		for (i=0; i<3; i++) {
			n_top[i] *= ir;
		}
	} else if (r2==0) {
		printf("Top wall preferred orientation wrong\n");
		exit(-1);
	}

        r2 = n_bot[0]*n_bot[0]+n_bot[1]*n_bot[1]+n_bot[2]*n_bot[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_bot[i] *= ir;
                }
        } else if (r2==0) {
                printf("Bottom wall preferred orientation wrong\n");
                exit(-1);
        }

        r2 = n_xlo[0]*n_xlo[0]+n_xlo[1]*n_xlo[1]+n_xlo[2]*n_xlo[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_xlo[i] *= ir;
                }
        } else if (r2==0) {
                printf("x bottom wall preferred orientation wrong\n");
                exit(-1);
        }

        r2 = n_xhi[0]*n_xhi[0]+n_xhi[1]*n_xhi[1]+n_xhi[2]*n_xhi[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_xhi[i] *= ir;
                }
        } else if (r2==0) {
                printf("x top wall preferred orientation wrong\n");
                exit(-1);
        }

        r2 = n_ylo[0]*n_ylo[0]+n_ylo[1]*n_ylo[1]+n_ylo[2]*n_ylo[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_ylo[i] *= ir;
                }
        } else if (r2==0) {
                printf("y bottom wall preferred orientation wrong\n");
                exit(-1);
        }

        r2 = n_yhi[0]*n_yhi[0]+n_yhi[1]*n_yhi[1]+n_yhi[2]*n_yhi[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_yhi[i] *= ir;
                }
        } else if (r2==0) {
                printf("y top wall preferred orientation wrong\n");
                exit(-1);
        }
	MPI_Barrier(MPI_COMM_WORLD);	
}

//Print title information
void title()
{
	if (myid == 0){
		printf("-----------------------------------------\n");
		printf("------3D LATTICE BOLTZMANN MPI CODE------\n");
		printf("--BASED ON ALGORITHM FROM YEOMANS ET AL--\n");
		printf("---AUTHORS: RUI ZHANG & TYLER ROBERTS----\n");
		printf("-----UNIVERSITY OF CHICAGO, JAN 2014-----\n");
		printf("-----------------------------------------\n");
	}
}
